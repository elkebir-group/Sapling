"""Place SNVs and copy-number events on a fixed topology, by alternating as k-means does.

The topology is a parent array with the normal clone at node 0. A placement puts every
called copy-number state of a segment (an event) on an edge, at most one per segment and
edge, and every SNV i on an edge e_i, with multiplicities m_{i,v} (the number of mutated
copies in clone v) among those allowed by e_i and the states of its segment (see
tree.snv_options). Each SNV takes the allowed multiplicities that fit its reads best.

Without clusters, the loss given the usages of every sample separates over segments, and
each segment moves one event at a time to its best edge (or swaps two events), with every
SNV of the segment at its best edge for each trial.
With clusters (all SNVs of a cluster share their edge), a round first moves every cluster
to its best edge given the usages and events, then optimizes the events of each segment
with the edges of the SNVs fixed.
Given the placement, fastppm-CN then solves for the usages of each sample. The loop stops
when a round no longer lowers the loss by more than the tolerance.
"""
import numpy as np

import fastppm

from .data import NORMAL
from .tree import snv_options


def _clip(p):
    return np.clip(p, 1e-12, 1 - 1e-12)


class Placer:
    def __init__(self, bulk, parent, kappa=-1.0, eps=1.0, cap=256):
        self.bulk, self.parent, self.kappa, self.eps, self.cap = bulk, list(parent), kappa, eps, cap
        self.k = len(parent)
        self.ch = [[] for _ in range(self.k)]
        for v in range(1, self.k):
            self.ch[parent[v]].append(v)
        order, stack = [], [0]
        while stack:
            x = stack.pop()
            order.append(x)
            stack.extend(reversed(self.ch[x]))
        self.preorder = order
        self.subtree = [None] * self.k
        for v in reversed(order):
            self.subtree[v] = [v] + [y for c in self.ch[v] for y in self.subtree[c]]
        S = bulk.num_segments
        self.events = [bulk.cn_events(s) for s in range(S)]
        self.V = [bulk.v[bulk.snvs_of[s]] for s in range(S)]
        self.D = [bulk.d[bulk.snvs_of[s]] for s in range(S)]
        self._options = {}

    # ---------------------------------------------------------------- states

    def states(self, nodes):
        """States of every node on a segment whose events are {state: node}."""
        at = {v: b for b, v in nodes.items()}
        st = [NORMAL] * self.k
        for v in self.preorder[1:]:
            st[v] = at.get(v, st[self.parent[v]])
        return st

    def valid(self, nodes, st):
        """At most one event per edge, each event changes the state, and lost alleles stay lost."""
        if len(set(nodes.values())) < len(nodes):
            return False
        for b, v in nodes.items():
            a0, b0 = st[self.parent[v]]
            if (a0, b0) == b or (a0 == 0 and b[0] > 0) or (b0 == 0 and b[1] > 0):
                return False
        return True

    # ------------------------------------------------------------ SNV options

    def options(self, st, nodes):
        """Allowed multiplicities of an SNV of the segment, one row per option, stacked over
        the edge where the SNV is introduced. Returns (O, edge of each row). Cached on the
        states of the subtree of the edge."""
        at = set(nodes.values())
        blocks, owner = [], []
        for w in range(1, self.k):
            key = (w, w in at, st[self.parent[w]], tuple(st[y] for y in self.subtree[w]))
            O = self._options.get(key)
            if O is None:
                O = snv_options(st, self.ch, self.parent, w, w in at, self.cap)
                self._options[key] = O
            blocks.append(O)
            owner.append(np.full(len(O), w))
        return np.vstack(blocks), np.concatenate(owner)

    # ------------------------------------------------------------------- loss

    def snv_table(self, s, nodes, u):
        """For the SNVs of segment s and every edge w: the binomial loss of the SNV introduced
        on w with its best allowed multiplicities (inf if none), and those multiplicities.
        Returns (A, best) with A of shape (SNVs, k) and best[i][w] a row of length k."""
        st = self.states(nodes)
        O, owner = self.options(st, nodes)
        prob = _clip((O @ u.T) / self.bulk.chat[:, s][None, :])     # options x samples
        V, D = self.V[s], self.D[s]
        nll = -(V @ np.log(prob).T + (D - V) @ np.log1p(-prob).T)     # SNVs x options
        A = np.full((len(V), self.k), np.inf)
        arg = np.zeros((len(V), self.k), dtype=int)
        start = 0
        for w in range(1, self.k):
            stop = start + int((owner == w).sum())
            if stop > start:
                j = start + np.argmin(nll[:, start:stop], axis=1)
                A[:, w] = nll[np.arange(len(V)), j]
                arg[:, w] = j
            start = stop
        return A, arg, O

    def segment_loss(self, s, nodes, u, snvs=True, edges=None):
        """Loss of segment s given usages u (P x k): squared error of the state fractions plus
        the binomial loss of its SNVs, each with its best allowed multiplicities, on its best
        edge or on the edge given by `edges`. Returns (loss, multiplicities, edges)."""
        b = self.bulk
        st = self.states(nodes)
        cn = 0.0
        for p in range(b.num_samples):
            g = {}
            for v in range(self.k):
                g[st[v]] = g.get(st[v], 0.0) + u[p, v]
            mu = b.mu[p][s]
            cn += b.weight[p, s] * sum((g.get(x, 0.0) - mu.get(x, 0.0)) ** 2 for x in set(g) | set(mu))
        if not snvs or len(b.snvs_of[s]) == 0:
            return cn, np.zeros((0, self.k), dtype=int), np.zeros(0, dtype=int)
        A, arg, O = self.snv_table(s, nodes, u)
        if edges is None:
            edges = np.argmin(A, axis=1)
        idx = np.arange(len(edges))
        return cn + float(A[idx, edges].sum()), O[arg[idx, edges]], edges

    def optimize_segment(self, s, nodes, u, snvs=True, edges=None, need=()):
        """Coordinate descent over the events of segment s from `nodes`: move one event to
        its best edge, or swap two events, while that lowers the loss. Every node in `need`
        (edges with no other event) must keep an event of this segment."""
        best = self.segment_loss(s, nodes, u, snvs, edges)
        improved = True
        while improved:
            improved = False
            trials = []
            for e in nodes:
                for x in range(1, self.k):
                    if x != nodes[e]:
                        trials.append({**nodes, e: x})
            ev = list(nodes)
            for i in range(len(ev)):
                for j in range(i + 1, len(ev)):
                    trials.append({**nodes, ev[i]: nodes[ev[j]], ev[j]: nodes[ev[i]]})
            for t in trials:
                if any(x not in t.values() for x in need) or not self.valid(t, self.states(t)):
                    continue
                r = self.segment_loss(s, t, u, snvs, edges)
                if r[0] < best[0] - 1e-9:
                    nodes, best, improved = t, r, True
        return nodes, best

    def random_nodes(self, rng, tries=1000):
        """A random valid placement of the events of every segment."""
        out = []
        for s in range(self.bulk.num_segments):
            for _ in range(tries):
                nodes = {e: int(x) for e, x in zip(self.events[s], rng.choice(np.arange(1, self.k), len(self.events[s]), replace=False))}
                if self.valid(nodes, self.states(nodes)):
                    break
            else:
                raise RuntimeError(f"no valid placement of the events of segment {s}")
            out.append(nodes)
        return out

    def occupancy(self, nodes, cl_edge, skip_segment=None):
        """Number of clusters and copy-number events on each edge."""
        occ = np.zeros(self.k, dtype=int)
        if cl_edge is not None:
            np.add.at(occ, cl_edge, 1)
        for s, x in enumerate(nodes):
            if s != skip_segment:
                np.add.at(occ, list(x.values()), 1)
        return occ

    def fill_empty(self, nodes, cl_edge, rng):
        """Move clusters or events from edges with several onto empty edges, until every edge
        carries at least one (one random valid move per empty edge)."""
        nodes, cl_edge = [dict(x) for x in nodes], np.array(cl_edge)
        for w in rng.permutation(np.arange(1, self.k)):
            occ = self.occupancy(nodes, cl_edge)
            if occ[w] > 0:
                continue
            items = [("c", c) for c in range(len(cl_edge)) if occ[cl_edge[c]] > 1]
            items += [("s", s, B) for s, x in enumerate(nodes) for B, v in x.items() if occ[v] > 1]
            for j in rng.permutation(len(items)):
                it = items[j]
                if it[0] == "c":
                    cl_edge[it[1]] = w
                    break
                trial = {**nodes[it[1]], it[2]: w}
                if self.valid(trial, self.states(trial)):
                    nodes[it[1]] = trial
                    break
            else:
                raise RuntimeError(f"cannot put an event on edge {w}")
        return nodes, cl_edge

    def random_start(self, rng):
        """Random usages, and a random placement with an event on every edge."""
        u0 = rng.dirichlet(np.ones(self.k), size=self.bulk.num_samples)
        cl = rng.integers(1, self.k, self.bulk.cluster.max() + 1) if self.bulk.cluster is not None else None
        nodes = self.random_nodes(rng)
        if cl is not None:
            nodes, cl = self.fill_empty(nodes, cl, rng)
        return u0, nodes, cl

    # ----------------------------------------------------------------- usages

    def solve(self, nodes, M, edge, snvs=True):
        """fastppm-CN on every sample, given the placement. Returns (loss, lower bound, u).
        With snvs=False, only the copy-number calls are fitted.

        A node whose edge carries no event has the genotype of its parent, so it is merged
        into its parent for the solve and gets usage 0."""
        b = self.bulk
        S = b.num_segments
        st = [self.states(nodes[s]) for s in range(S)]
        busy = np.zeros(self.k, dtype=bool)
        busy[0] = True
        busy[edge] = True
        for x in nodes:
            busy[list(x.values())] = True
        keep = [v for v in self.preorder if busy[v]]
        new = {v: j for j, v in enumerate(keep)}

        def kept_ancestor(v):
            while not busy[v]:
                v = self.parent[v]
            return v
        par = [-1] + [new[kept_ancestor(self.parent[v])] for v in keep[1:]]
        st_kept = [[x[v] for v in keep] for x in st]
        rows = np.arange(len(b.seg)) if snvs else np.zeros(0, dtype=int)
        if snvs:
            M_kept = np.ascontiguousarray(M[:, keep], dtype=np.int32)
        else:
            M_kept = np.zeros((0, len(keep)), dtype=np.int32)  # M is not needed (it may be None)
        seg = np.ascontiguousarray(b.seg[rows], dtype=np.int32)
        L, lb, u = 0.0, 0.0, np.zeros((b.num_samples, self.k))
        for p in range(b.num_samples):
            res = fastppm.regress_copy_number(
                parent=par,
                states=st_kept,
                multiplicity=M_kept,
                segment=seg,
                variant=np.ascontiguousarray(b.v[rows, p], dtype=np.int32),
                total=np.ascontiguousarray(b.d[rows, p], dtype=np.int32),
                state_fractions=[[(x, y, float(m)) for (x, y), m in b.mu[p][s].items()] for s in range(S)],
                mean_copy_number=b.chat[p].tolist(), segment_weight=b.weight[p].tolist(),
                kappa=self.kappa, eps=self.eps)
            # calls of states that no node has: fastppm-CN leaves them out, their error is mu^2
            missing = sum(b.weight[p, s] * m * m for s in range(S) for x, m in b.mu[p][s].items() if x not in set(st[s]))
            L += res["loss"] + missing
            lb += res["lower_bound"] + missing
            u[p, keep] = res["usages"]
        return L, lb, u

    # ------------------------------------------------------------- alternation

    def cluster_step(self, nodes, cl_edge, u):
        """Move each cluster in turn to its best edge given the usages and the events,
        never emptying the edge it leaves."""
        b = self.bulk
        cost = np.zeros((b.cluster.max() + 1, self.k))
        for s in range(b.num_segments):
            if len(b.snvs_of[s]):
                A, _, _ = self.snv_table(s, nodes[s], u)
                np.add.at(cost, b.cluster[b.snvs_of[s]], A)
        cost[:, 0] = np.inf
        cl_edge = np.array(cl_edge)
        occ = self.occupancy(nodes, cl_edge)
        for c in range(len(cl_edge)):
            if cl_edge[c] != 0 and occ[cl_edge[c]] == 1:
                continue  # the only event of its edge
            w = int(np.argmin(cost[c]))
            occ[cl_edge[c]] -= 1
            occ[w] += 1
            cl_edge[c] = w
        return cl_edge

    def run(self, u, nodes, cl_edge=None, snvs=True, max_rounds=50, tol=None, log=None):
        """Alternate from usages u (P x k), a placement of the events, and (with clusters) the
        edge of every cluster, until a round lowers the loss by at most tol (default: the gap
        tolerance of fastppm-CN, 2 eps per sample). With clusters, every edge keeps at least
        one cluster or copy-number event. Returns the best round, as dict(loss, lower, u,
        nodes, cl_edge, M, edge, rounds, history). With snvs=False, only the copy-number
        calls are fitted."""
        b = self.bulk
        n = len(b.seg)
        tol = 2 * self.eps * b.num_samples if tol is None else tol
        clustered = snvs and b.cluster is not None
        if clustered:
            if cl_edge is None:
                raise ValueError("with clusters, run needs the starting edge of every cluster")
            if (self.occupancy(nodes, cl_edge)[1:] == 0).any():
                raise ValueError("every edge needs a cluster or a copy-number event at the start")
        best, history = None, []
        for r in range(1, max_rounds + 1):
            M = np.zeros((n, self.k), dtype=int)
            edge = np.zeros(n, dtype=int)
            fixed = None
            if clustered:
                cl_edge = self.cluster_step(nodes, cl_edge, u)
                fixed = cl_edge[b.cluster]
            new_nodes, step = list(nodes), 0.0
            for s in range(b.num_segments):
                es = None if fixed is None else fixed[b.snvs_of[s]]
                need = ()
                if clustered:
                    occ = self.occupancy(new_nodes, cl_edge, skip_segment=s)
                    need = [v for v in range(1, self.k) if occ[v] == 0]
                ns, (Ls, rows, owner) = self.optimize_segment(s, new_nodes[s], u, snvs, es, need)
                new_nodes[s] = ns
                if snvs:
                    M[b.snvs_of[s]] = rows
                    edge[b.snvs_of[s]] = owner
                step += Ls
            nodes = new_nodes
            L, lb, u = self.solve(nodes, M, edge if snvs else np.zeros(0, dtype=int), snvs)
            history.append((step, L, lb))
            if log:
                log(f"    round {r}: placement step {step:.1f}, after solving {L:.1f} (lower bound {lb:.1f})")
            improved = best is None or L < best["loss"] - tol
            if best is None or L < best["loss"]:
                best = dict(loss=L, lower=lb, u=u, nodes=nodes, cl_edge=None if cl_edge is None else np.array(cl_edge),
                            M=M, edge=edge)
            if not improved:
                break
        best.update(rounds=r, history=history)
        return best

    def score(self, nodes, cl_edge, u):
        """Loss after choosing the multiplicities of every SNV given usages u (the SNV edges
        follow the clusters) and solving for the usages. Returns (loss, lower bound, u, M, edge)."""
        b = self.bulk
        edge = cl_edge[b.cluster]
        M = np.zeros((len(b.seg), self.k), dtype=int)
        for s in range(b.num_segments):
            if len(b.snvs_of[s]):
                _, rows, _ = self.segment_loss(s, nodes[s], u, True, edge[b.snvs_of[s]])
                M[b.snvs_of[s]] = rows
        L, lb, uu = self.solve(nodes, M, edge)
        return L, lb, uu, M, edge

    def swap(self, nodes, cl_edge, v, w):
        """The placement with the bundles of edges v and w exchanged (all clusters and
        copy-number events on them), or None if a copy-number event becomes invalid."""
        cl = np.array(cl_edge)
        cl[cl_edge == v], cl[cl_edge == w] = w, v
        out = []
        for x in nodes:
            y = {B: (w if e == v else v if e == w else e) for B, e in x.items()}
            if y != x and not self.valid(y, self.states(y)):
                return None
            out.append(y)
        return out, cl

    def search(self, u, nodes, cl_edge, max_sweeps=50, log=None):
        """Alternate until the loss stops falling, then try the bundle swap of every pair of
        edges and keep the best, if it lowers the loss by more than the tolerance; repeat.
        Returns the result of the last alternation, with the number of swaps kept."""
        tol = 2 * self.eps * self.bulk.num_samples
        res = self.run(u, nodes, cl_edge, log=log)
        swaps = 0
        for _ in range(max_sweeps):
            best = None
            for v in range(1, self.k):
                for w in range(v + 1, self.k):
                    t = self.swap(res["nodes"], res["cl_edge"], v, w)
                    if t is None:
                        continue
                    L = self.score(t[0], t[1], res["u"])
                    if best is None or L[0] < best[0][0]:
                        best = (L, t, (v, w))
            if best is None or best[0][0] >= res["loss"] - tol:
                break
            swaps += 1
            if log:
                log(f"    swap the bundles of edges {best[2][0]} and {best[2][1]}: {res['loss']:.1f} -> {best[0][0]:.1f}")
            res = self.run(best[0][2], best[1][0], best[1][1], log=log)
        res["swaps"] = swaps
        return res

    def init_from_calls(self, rng, restarts=20, log=None):
        """Usages and event placement fitted to the copy-number calls alone: the best of
        `restarts` alternations from random usages and placements."""
        best = None
        for _ in range(restarts):
            u0 = rng.dirichlet(np.ones(self.k), size=self.bulk.num_samples)
            res = self.run(u0, self.random_nodes(rng), snvs=False)
            if best is None or res["loss"] < best["loss"]:
                best = res
        if log:
            log(f"    fit to the copy-number calls: loss {best['loss']:.2f}")
        nodes, cl = best["nodes"], None
        if self.bulk.cluster is not None:
            # clusters to their best edges given these usages, then an event on every edge
            cl = self.cluster_step(nodes, np.zeros(self.bulk.cluster.max() + 1, dtype=int), best["u"])
            nodes, cl = self.fill_empty(nodes, cl, rng)
        return best["u"], nodes, cl
