"""Sapling-CN: grow a clone tree one clone at a time (beam search), then refine it by moving events.

Events are SNV clusters and copy-number events (segment s acquires state B). The search
starts from the normal clone and one tumor clone that carries every SNV cluster. Each
step adds one clone, as a new leaf or by splitting an edge, and gives it one event: an
event moved from another edge, or a copy-number event not yet in the tree. Events then
move onto the new clone while that lowers the loss. Every tree scores all SNVs and all
copy-number calls, and trees in a step have the same number of clones, so their losses
are compared directly. The search stops at the requested number of clones.
"""
import multiprocessing as mp
import time

import numpy as np

from .score import score
from .tree import Tree

_DATA, _KW = None, None


def _init(data, kw):
    global _DATA, _KW
    _DATA, _KW = data, kw


def _score(args):
    tree, u0 = args
    r = score(_DATA, tree, u0, **_KW)
    return r["loss"], r["lower"], r["u"]


def tied(beam):
    """Trees of the beam that the certificate cannot rank below the best: lower bound <= best loss."""
    best = beam[0][2]
    return [b for b in beam if b[3] <= best]


def all_events(data):
    return [("c", c) for c in range(len(data.clusters))] + [("s", s, b) for s, b in data.cn_events()]


def _dedupe(items, tree_of=lambda x: x[0]):
    seen, out = set(), []
    for x in items:
        key = tree_of(x).key()
        if key not in seen:
            seen.add(key)
            out.append(x)
    return out


class Search:
    def __init__(self, data, beam=25, workers=8, kappa=1.0, eps=1.0, rounds=3, log=print):
        self.data, self.beam, self.log = data, beam, log
        self.events = all_events(data)
        kw = dict(kappa=kappa, eps=eps, rounds=rounds)
        self.tol = 2 * eps
        _init(data, kw)
        self.pool = mp.get_context("fork").Pool(workers, _init, (data, kw)) if workers > 1 else None
        self.solves = 0

    def _map(self, cands):
        self.solves += len(cands)
        if self.pool is None:
            return [_score(c) for c in cands]
        return self.pool.map(_score, cands, chunksize=max(1, len(cands) // (8 * self.pool._processes)))

    def _new_clone_candidates(self, tree, u):
        """Trees with one more clone, which carries one event moved or added to it."""
        S, out = self.data.num_segments, []
        for e in self.events:
            t_minus, u_minus, contracted = tree.remove(e, u)
            if contracted is not None:
                continue  # e was the only event of its clone: moving it does not add a clone
            out += t_minus.expand(e, u_minus, S, new_clone=True, existing=False)
        return out

    def _absorb(self, tree, u, loss, lower, w, max_rounds=20):
        """Move events onto the edge into w, and add missing copy-number events to any edge.
        Each round scores every single move on the current tree, then applies the improving
        ones greedily in order of gain, each checked again on the tree as it is by then.
        Rounds repeat until no move lowers the loss by more than tol."""
        S = self.data.num_segments

        def apply(t, uu, e, v):
            tm, um, contracted = t.remove(e, uu)
            if contracted is not None:
                return None
            t2 = tm.place(e, v)
            return None if t2 is None or not t2.valid(S) else (t2, um)

        for _ in range(max_rounds):
            node = tree.node_of()
            moves = [(e, v) for e in self.events if node.get(e) != w
                     for v in ([w] if e in node else range(1, tree.k))]
            trial = [(m, apply(tree, u, *m)) for m in moves]
            trial = [(m, c) for m, c in trial if c is not None]
            if not trial:
                break
            res = self._map([c for _, c in trial])
            base, accepted = loss, 0
            for j in sorted(range(len(trial)), key=lambda j: res[j][0]):
                if res[j][0] >= base - self.tol:
                    break  # no gain on its own
                cand = apply(tree, u, *trial[j][0])  # on the current tree, which earlier moves changed
                if cand is None:
                    continue
                r = self._map([cand])[0]
                if r[0] < loss - self.tol:
                    tree, u, loss, lower = cand[0], r[2], r[0], r[1]
                    accepted += 1
            if accepted == 0:
                break
        return tree, u, loss, lower

    def grow(self, k):
        """Beam search from 2 to k clones (the normal clone included). Returns the final beam
        as (tree, usages, loss, lower bound), best first."""
        t0 = Tree.trunk(len(self.data.clusters))
        r = self._map([(t0, np.array([0.2, 0.8]))])[0]
        beam = [(t0, r[2], r[0], r[1])]
        start = time.time()
        for kk in range(3, k + 1):
            cands = []
            for t, u, *_ in beam:
                cands += self._new_clone_candidates(t, u)
            cands = _dedupe(cands)
            res = self._map(cands)
            ranked = sorted(zip(cands, res), key=lambda x: x[1][0])
            kept = [(t, r[2], r[0], r[1]) for (t, _), r in ranked[: self.beam] if np.isfinite(r[0])]
            if not kept:
                raise RuntimeError(f"no valid tree with {kk} clones")
            beam = [self._absorb(t, u, loss, lower, kk - 1) for t, u, loss, lower in kept]
            beam = sorted(_dedupe(beam), key=lambda b: b[2])
            if self.log:
                self.log(f"  {kk} clones: {len(cands)} candidates, best {beam[0][2]:.1f}, {time.time() - start:.1f}s")
        return beam

    def refine(self, tree, u, loss, max_sweeps=3):
        """Move one event at a time to its best position (or drop a copy-number event), keeping
        at most as many clones as the tree has; accept moves that lower the loss by > tol."""
        S, k = self.data.num_segments, tree.k
        for sweep in range(max_sweeps):
            moved = 0
            for e in self.events:
                placed = e in tree.node_of()
                t_minus, u_minus, _ = tree.remove(e, u)
                cands = [(t2, u2) for t2, u2 in t_minus.expand(e, u_minus, S) if t2.k <= k]
                if e[0] == "s" and placed:
                    cands.append((t_minus, u_minus))
                cands = _dedupe(cands)
                if not cands:
                    continue
                res = self._map(cands)
                j = int(np.argmin([r[0] for r in res]))
                if res[j][0] < loss - self.tol:
                    tree, u, loss = cands[j][0], res[j][2], res[j][0]
                    moved += 1
            if self.log:
                self.log(f"  refine sweep {sweep + 1}: {moved} moves, loss {loss:.1f}")
            if moved == 0:
                break
        return tree, u, loss

    def close(self):
        if self.pool is not None:
            self.pool.close()
