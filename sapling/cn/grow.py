"""Sapling-CN on bulk input: grow a clone tree one clone at a time (beam search).

    python -m sapling.cn.grow <instance dir> [--clones k|true] [--beam 10] [--workers 8] [--json out.json]

The input is snvs.tsv, cn.tsv and clusters.tsv (from clonesim's generatebulk), with any
number of samples. Events are SNV clusters (all SNVs of a cluster are introduced on the
same edge) and copy-number events (segment s acquires a called state B). Every edge
carries at least one event.

A tree is scored by choosing, for every SNV, the multiplicities that fit its reads best
among those allowed by its edge and the states of its segment, then solving for the
usages of every sample with fastppm-CN, and alternating the two (at most `rounds` times).
The loss is the sum over the samples. Scoring runs in the C++ core (sapling.cn._core, built
from cpp/), on several threads; score() below is the same scoring in Python.

The search starts from the normal clone and one tumor clone that carries every cluster.
Each step adds one clone, as a new leaf or by splitting an edge, carrying one event moved
from another edge or a copy-number event not yet in the tree, and keeps the best `beam`
trees. Each kept tree is then polished: the best move of whole clones or of single events
(at the same number of clones) is applied while one lowers the loss. A move lowers the loss
when its loss is below the lower bound of the current tree from fastppm-CN. The search stops
at the requested number of clones, and the best tree is polished once more.

Candidates are screened with one round of alternation and at most 10 iterations of
fastppm-CN per sample; the best of them are then scored fully.
"""
import argparse
import json
import time

import numpy as np

import fastppm

from . import _core
from .data import NORMAL, load_bulk, load_truth
from .run import ancestor_agreement, fmt, relations
from .tree import Tree, snv_option_kinds, snv_options


# ---------------------------------------------------------------------- scoring

def _clip(p):
    return np.clip(p, 1e-12, 1 - 1e-12)


def solve(bulk, tree, M, kappa=-1.0, eps=1.0, u_init=None):
    """fastppm-CN on every sample, started from the usages u_init (P x k) if given. Returns
    (loss, lower bound, u of shape P x k, converged)."""
    S = bulk.num_segments
    st = tree.states(S)
    M = np.ascontiguousarray(M, dtype=np.int32)
    seg = np.ascontiguousarray(bulk.seg, dtype=np.int32)
    L = lb = 0.0
    ok = True
    u = np.zeros((bulk.num_samples, tree.k))
    for p in range(bulk.num_samples):
        res = fastppm.regress_copy_number(
            parent=tree.parent, states=st, multiplicity=M, segment=seg,
            variant=np.ascontiguousarray(bulk.v[:, p], dtype=np.int32),
            total=np.ascontiguousarray(bulk.d[:, p], dtype=np.int32),
            state_fractions=[[(x, y, float(m)) for (x, y), m in bulk.mu[p][s].items()] for s in range(S)],
            mean_copy_number=bulk.chat[p].tolist(), segment_weight=bulk.weight[p].tolist(),
            kappa=kappa, eps=eps, initial_usages=None if u_init is None else np.asarray(u_init[p], dtype=float).tolist())
        # calls of states that no clone has: fastppm-CN leaves them out, their error is mu^2
        missing = sum(bulk.weight[p, s] * m * m for s in range(S) for x, m in bulk.mu[p][s].items() if x not in set(st[s]))
        L += res["loss"] + missing
        lb += res["lower_bound"] + missing
        u[p] = res["usages"]
        ok &= res["converged"]
    return L, lb, u, ok


def choose_multiplicities(bulk, tree, u, options):
    """For every SNV, the allowed multiplicities (over the clones) that fit its reads best."""
    M = np.zeros((len(bulk.seg), tree.k), dtype=int)
    for (c, s), (idx, O) in options.items():
        prob = _clip((O @ u.T) / bulk.chat[:, s][None, :])                        # options x samples
        V, D = bulk.v[idx], bulk.d[idx]
        nll = -(V @ np.log(prob).T + (D - V) @ np.log1p(-prob).T)                    # SNVs x options
        M[idx] = O[np.argmin(nll, axis=1)]
    return M


def score(bulk, tree, u0, rounds=3, kappa=-1.0, eps=1.0, M=None):
    """Loss of a tree, alternating the multiplicities and the usages from usages u0
    (P x k). With M given, the multiplicities are fixed. Returns dict(loss, lower, u, M)."""
    S = bulk.num_segments
    if M is not None:
        L, lb, u, ok = solve(bulk, tree, M, kappa, eps, u_init=u0)
        return dict(loss=L, lower=lb, u=u, M=M, converged=ok)
    st, ch = tree.states(S), tree.children()
    options = {}
    for w in range(1, tree.k):
        for c in tree.clusters[w]:
            for s, idx in bulk.groups[c].items():
                O = snv_options(st[s], ch, tree.parent, w, s in tree.cn[w])
                if len(O) == 0:
                    return dict(loss=np.inf, lower=np.inf, u=u0, M=None, converged=False)
                options[(c, s)] = (idx, O)
    u, prev = np.asarray(u0, dtype=float), None
    for _ in range(rounds):
        M = choose_multiplicities(bulk, tree, u, options)
        if prev is not None and np.array_equal(M, prev):
            break
        L, lb, u, ok = solve(bulk, tree, M, kappa, eps, u_init=u)  # warm start from the last usages
        prev = M
    return dict(loss=L, lower=lb, u=u, M=M, converged=ok)


# ----------------------------------------------------------------------- search

def make_scorer(bulk, kappa=-1.0, eps=1.0, threads=8):
    """The C++ scoring core (sapling.cn._core) for this input."""
    return _core.Scorer(
        segment=bulk.seg, variant=bulk.v, total=bulk.d, cluster=bulk.cluster,
        calls=[[[(x, y, float(m)) for (x, y), m in bulk.mu[p][s].items()] for s in range(bulk.num_segments)]
               for p in range(bulk.num_samples)],
        mean_copy_number=bulk.chat.tolist(), segment_weight=bulk.weight.tolist(), kappa=kappa, eps=eps, threads=threads)


def _encode(tree):
    return (tree.parent, [sorted(c) for c in tree.clusters], [[(s, b[0], b[1]) for s, b in x.items()] for x in tree.cn])


def all_events(bulk):
    return ([("c", c) for c in range(bulk.cluster.max() + 1)]
            + [("s", s, b) for s in range(bulk.num_segments) for b in bulk.cn_events(s)])


def _dedupe(items, tree_of=lambda x: x[0]):
    seen, out = set(), []
    for x in items:
        key = tree_of(x).key()
        if key not in seen:
            seen.add(key)
            out.append(x)
    return out


class Search:
    def __init__(self, bulk, beam=10, workers=8, kappa=-1.0, eps=1.0, rounds=3, screen=3, log=print):
        """Candidates are screened with one round of alternation; the best screen * beam of a
        growth step (screen * 5 when polishing) are scored fully."""
        self.bulk, self.beam, self.screen, self.log = bulk, beam, screen, log
        bulk.groups = [{} for _ in range(bulk.cluster.max() + 1)]
        for c in range(len(bulk.groups)):
            idx = np.flatnonzero(bulk.cluster == c)
            for s in np.unique(bulk.seg[idx]):
                bulk.groups[c][int(s)] = idx[bulk.seg[idx] == s]
        self.events = all_events(bulk)
        self.rounds = rounds
        self.scorer = make_scorer(bulk, kappa, eps, threads=workers)
        self.solves = 0

    def _map(self, cands, rounds=None, max_iterations=100000):
        """(loss, lower bound, usages) of every (tree, initial usages)."""
        self.solves += len(cands)
        if not cands:
            return []
        loss, lower, u = self.scorer.score([_encode(t) for t, _ in cands],
                                           [np.asarray(uu, dtype=float) for _, uu in cands], rounds or self.rounds,
                                           max_iterations)
        return list(zip(loss.tolist(), lower.tolist(), u))

    def _screened(self, cands, keep):
        """Score every candidate with one round of alternation and at most 10 iterations of
        fastppm-CN per sample (its loss is that of feasible usages), then the best `keep` of
        them fully. Returns the full results, the others marked infinite."""
        cheap = self._map(cands, rounds=1, max_iterations=10)
        top = np.argsort([r[0] for r in cheap])[:keep]
        full = self._map([(cands[j][0], cheap[j][2]) for j in top])  # started from the screen's usages
        out = [(np.inf, np.inf, None)] * len(cands)
        for j, r in zip(top, full):
            out[j] = r
        return out

    def _new_clone_candidates(self, tree, u):
        """Trees with one more clone, which carries one event moved or added to it."""
        S, out = self.bulk.num_segments, []
        for e in self.events:
            t_minus, u_minus, contracted = tree.remove(e, u)
            if contracted is not None:
                continue  # e was the only event of its clone: moving it does not add a clone
            out += t_minus.expand(e, u_minus, S, new_clone=True, existing=False)
        return out

    def grow(self, k):
        """Beam search from 2 to k clones (the normal clone included). Returns the final beam
        as (tree, usages, loss, lower bound), best first. self.history gets the best tree of
        every step, as (tree, usages, loss, lower bound)."""
        P = self.bulk.num_samples
        t0 = Tree.trunk(self.bulk.cluster.max() + 1)
        r = self._map([(t0, np.tile([0.2, 0.8], (P, 1)))])[0]
        beam = [(t0, r[2], r[0], r[1])]
        self.history = [beam[0]]
        start = time.time()
        for kk in range(3, k + 1):
            cands = []
            for t, u, *_ in beam:
                cands += self._new_clone_candidates(t, u)
            cands = _dedupe(cands)
            res = self._screened(cands, self.screen * self.beam)
            ranked = sorted(zip(cands, res), key=lambda x: x[1][0])
            kept = [(t, r[2], r[0], r[1]) for (t, _), r in ranked[: self.beam] if np.isfinite(r[0])]
            if not kept:
                raise RuntimeError(f"no valid tree with {kk} clones")
            beam = [self.polish(t, u, loss, lower, log=False)[:4] for t, u, loss, lower in kept]
            beam = sorted(_dedupe(beam), key=lambda b: b[2])
            self.history.append(beam[0])
            if self.log:
                self.log(f"  {kk} clones: {len(cands)} candidates, best {beam[0][2]:.1f}, {time.time() - start:.1f}s")
        return beam

    # ------------------------------------------------------------- polishing

    def _ok(self, tree):
        """Lost alleles stay lost, and every copy-number event changes the state of its clone."""
        S = self.bulk.num_segments
        if not tree.valid(S):
            return False
        st = tree.states(S)
        return all(st[s][tree.parent[v]] != b for v in range(1, tree.k) for s, b in tree.cn[v].items())

    def _clone_moves(self, tree, u):
        """Trees with the same number of clones, from moves of whole clones:
        (a) reattach a clone, with its subtree and everything on its edge, under another parent;
        (b) remove a clone (its children attach to its parent) and put everything on its edge on
            a new clone, as a leaf or splitting an edge;
        (c) swap everything on the edges of two clones."""
        k, out = tree.k, []
        below = [None] * k
        for v in reversed(tree.preorder()):
            below[v] = {v}.union(*(below[c] for c in tree.children()[v])) if tree.children()[v] else {v}
        for v in range(1, k):                                                    # (a)
            for p in range(k):
                if p != tree.parent[v] and p not in below[v]:
                    par = list(tree.parent)
                    par[v] = p
                    out.append((Tree(par, tree.clusters, tree.cn), np.array(u)))
        for v in range(1, k):                                                    # (b)
            keep = [w for w in range(k) if w != v]
            new = {w: j for j, w in enumerate(keep)}
            pv = tree.parent[v]
            par = [-1 if tree.parent[w] < 0 else new[pv if tree.parent[w] == v else tree.parent[w]] for w in keep]
            uu = np.array(u)
            uu[..., pv] += uu[..., v]
            uu = uu[..., keep]
            base = Tree(par, [tree.clusters[w] for w in keep], [tree.cn[w] for w in keep])
            ch = base.children()
            for p in range(base.k):
                for c in [None] + ch[p]:
                    if p == new[pv] and (c is None and not tree.children()[v]):
                        continue  # puts v back where it was
                    t = Tree(base.parent + [p], base.clusters + [tree.clusters[v]], base.cn + [tree.cn[v]])
                    if c is not None:
                        t.parent[c] = base.k
                    out.append((t, np.append(uu, np.zeros_like(uu[..., :1]), axis=-1)))
        for v in range(1, k):                                                    # (c)
            for w in range(v + 1, k):
                cl, cn = list(tree.clusters), list(tree.cn)
                cl[v], cl[w], cn[v], cn[w] = cl[w], cl[v], cn[w], cn[v]
                out.append((Tree(tree.parent, cl, cn), np.array(u)))
        return out

    def _event_moves(self, tree, u):
        """Trees with the same number of clones, from moving one event (or dropping a
        copy-number event)."""
        S, out = self.bulk.num_segments, []
        node = tree.node_of()
        for e in self.events:
            t_minus, u_minus, _ = tree.remove(e, u)
            out += [(t2, u2) for t2, u2 in t_minus.expand(e, u_minus, S) if t2.k == tree.k]
            if e[0] == "s" and e in node and t_minus.k == tree.k:
                out.append((t_minus, u_minus))
        return out

    def polish(self, tree, u, loss, lower, max_sweeps=50, log=True):
        """Apply the best of all clone moves and single-event moves while its loss is below
        the lower bound of the current tree (so the move certainly lowers the minimum of
        fastppm-CN, for the chosen multiplicities). Returns (tree, u, loss, lower, moves)."""
        moves = []
        for sweep in range(max_sweeps):
            t0 = time.time()
            cands = [(t, uu, "clone") for t, uu in self._clone_moves(tree, u)]
            cands += [(t, uu, "event") for t, uu in self._event_moves(tree, u)]
            cands = [c for c in _dedupe(cands) if c[0].key() != tree.key() and self._ok(c[0])]
            res = self._screened([(t, uu) for t, uu, _ in cands], self.screen * 5)
            j = int(np.argmin([r[0] for r in res]))
            if self.log and log:
                self.log(f"  polish sweep {sweep + 1}: {len(cands)} trees, best {res[j][0]:.1f} "
                         f"({res[j][0] - loss:+.1f}, {cands[j][2]} move), {time.time() - t0:.1f}s")
            if res[j][0] >= lower:
                break
            tree, u, loss, lower = cands[j][0], res[j][2], res[j][0], res[j][1]
            moves.append(cands[j][2])
        return tree, u, loss, lower, moves

    def close(self):
        pass


# ------------------------------------------------------------------- evaluation

def true_tree(bulk, truth):
    """The true tree as a Tree: each cluster on the edge where its SNVs are introduced, and
    every change of state on its edge (called or not)."""
    parent, k = truth["parent"], len(truth["parent"])
    clusters = [set() for _ in range(k)]
    for c in range(bulk.cluster.max() + 1):
        clusters[int(truth["arises"][bulk.cluster == c][0])].add(c)
    st = truth["states"]
    cn = [{s: st[s][v] for s in range(bulk.num_segments) if v > 0 and st[s][v] != st[s][parent[v]]} for v in range(k)]
    return Tree(parent, clusters, cn)


def agreement(t, truth, events):
    """Fraction of ordered pairs of events (placed in both trees) with the same ancestor
    relation: over all pairs, and by type (cc: clusters, ss: copy-number events, cs: mixed)."""
    r1, r2 = relations(t, events), relations(truth, events)
    out = ancestor_agreement(t, truth, events)
    common = [p for p in r1 if p in r2]
    out["all"] = float(np.mean([r1[p] == r2[p] for p in common])) if common else float("nan")
    return out


def order_counts(bulk, tree, M):
    """For every copy-number event (s, B) on an edge (p -> w) that also carries clusters:
    {(w, s): {c: [before, undetermined, after, unmatched]}}, counting the SNVs of each cluster
    c on the edge that lie on s. An SNV arose before the change if its multiplicities M[i]
    (over the clones) come only from the options of an SNV introduced before it, after the
    change if only from those introduced after it, and is undetermined if from both. Unmatched
    counts SNVs whose multiplicities are not among the options (none, for an inferred tree).
    Clusters without SNVs on s are left out."""
    S = bulk.num_segments
    st, ch = tree.states(S), tree.children()
    out = {}
    for w in range(1, tree.k):
        for s in tree.cn[w]:
            kinds = None
            for c in sorted(tree.clusters[w]):
                idx = np.flatnonzero((bulk.cluster == c) & (bulk.seg == s))
                if len(idx) == 0:
                    continue
                if kinds is None:
                    kinds = snv_option_kinds(st[s], ch, tree.parent, w, True)
                n = [0, 0, 0, 0]
                for i in idx:
                    b, a = kinds.get(tuple(int(x) for x in M[i]), (None, None))
                    n[3 if b is None else 1 if (b and a) else 0 if b else 2] += 1
                out.setdefault((w, s), {})[c] = n
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("instance")
    ap.add_argument("--clones", default="true", help="number of clones, the normal clone included ('true': as in the true tree)")
    ap.add_argument("--beam", type=int, default=10)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--kappa", type=float, default=-1.0)
    ap.add_argument("--eps", type=float, default=1.0)
    ap.add_argument("--screen", type=int, default=3)
    ap.add_argument("--json", type=str)
    a = ap.parse_args()

    bulk = load_bulk(a.instance)
    truth = load_truth(a.instance, bulk)
    srch = Search(bulk, beam=a.beam, workers=a.workers, kappa=a.kappa, eps=a.eps, screen=a.screen)
    tt = true_tree(bulk, truth)
    gt = score(bulk, tt, truth["u"], M=truth["M"], kappa=a.kappa, eps=a.eps)
    k = tt.k if a.clones == "true" else int(a.clones)
    events = srch.events
    print(f"{a.instance}: {tt.k} clones, {bulk.cluster.max() + 1} clusters, {len(events) - bulk.cluster.max() - 1} "
          f"copy-number events, {bulk.num_samples} samples; ground truth loss {gt['loss']:.1f}; searching for {k} clones, beam {a.beam}")
    t0 = time.time()
    beam = srch.grow(k)
    tree, u, loss, lower = beam[0]
    t_grow = time.time() - t0
    grown = agreement(tree, tt, events)
    print(f"  grown: loss {loss:.1f} ({loss - gt['loss']:+.1f} from the ground truth), agreement {fmt(grown)}, {t_grow:.0f}s")
    tree, u, loss, lower, moves = srch.polish(tree, u, loss, lower)
    # final scoring, for the multiplicities that go with the usages
    lf, _, uf, Mf = srch.scorer.score([_encode(tree)], [np.asarray(u, dtype=float)], srch.rounds, 100000, True)
    u, M = uf[0], Mf[0]
    if abs(lf[0] - loss) > 0.05:  # the scoring starts from the final usages, so it can differ slightly
        print(f"  final scoring: loss {lf[0]:.1f}")
        loss = float(lf[0])
    dt = time.time() - t0
    srch.close()
    final = agreement(tree, tt, events)
    print(f"  polished: loss {loss:.1f} ({loss - gt['loss']:+.1f} from the ground truth), {len(moves)} moves, "
          f"agreement {fmt(final)}, {srch.solves} trees scored, {dt:.0f}s")
    if a.json:
        json.dump(dict(instance=a.instance, clones=k, gt_loss=gt["loss"], loss_grown=beam[0][2], loss=loss,
                       agree_grown=grown, agree=final, solves=srch.solves, time_grow=t_grow, time=dt, moves=moves,
                       tree=dict(parent=tree.parent, clusters=[sorted(x) for x in tree.clusters],
                                 cn=[{str(s): list(b) for s, b in x.items()} for x in tree.cn], u=u.tolist(),
                                 M=M.tolist()),
                       history=[dict(parent=t.parent, clusters=[sorted(x) for x in t.clusters],
                                     cn=[{str(s): list(b) for s, b in x.items()} for x in t.cn], u=np.asarray(uu).tolist(),
                                     loss=float(ls)) for t, uu, ls, _ in srch.history]),
                  open(a.json, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
