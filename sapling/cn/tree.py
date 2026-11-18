"""Clone trees whose edges carry events.

Node 0 is the normal clone, in state (1, 1) on every segment. Every other edge carries
SNV clusters (any number) and copy-number events (at most one per segment). A
copy-number event (s, B) sets segment s to state B, and a node's state on s is that of
the nearest event on s at or above it. Each state is introduced by at most one event
(no homoplasy), and an allele with zero copies stays lost. A copy-number event may be
absent from the tree; its call then counts as unexplained in the loss.
"""
from itertools import product

import numpy as np

from .data import NORMAL


class Tree:
    __slots__ = ("parent", "clusters", "cn", "_states")

    def __init__(self, parent, clusters, cn):
        self.parent = list(parent)
        self.clusters = [set(x) for x in clusters]
        self.cn = [dict(x) for x in cn]
        self._states = None

    @classmethod
    def trunk(cls, num_clusters):
        """The normal clone and one tumor clone that carries every SNV cluster."""
        return cls([-1, 0], [set(), set(range(num_clusters))], [{}, {}])

    @property
    def k(self):
        return len(self.parent)

    def copy(self):
        return Tree(self.parent, self.clusters, self.cn)

    def children(self):
        ch = [[] for _ in range(self.k)]
        for v in range(1, self.k):
            ch[self.parent[v]].append(v)
        return ch

    def preorder(self):
        ch, out, stack = self.children(), [], [0]
        while stack:
            x = stack.pop()
            out.append(x)
            stack.extend(reversed(ch[x]))
        return out

    def states(self, num_segments):
        """states[s][v] for every segment and node."""
        if self._states is None:
            st = [[NORMAL] * self.k for _ in range(num_segments)]
            for v in self.preorder()[1:]:
                p = self.parent[v]
                for s in range(num_segments):
                    st[s][v] = self.cn[v].get(s, st[s][p])
            self._states = st
        return self._states

    def valid(self, num_segments):
        """An allele with zero copies cannot be regained."""
        st = self.states(num_segments)
        for v in range(1, self.k):
            p = self.parent[v]
            for s in self.cn[v]:
                (a0, b0), (a1, b1) = st[s][p], st[s][v]
                if (a0 == 0 and a1 > 0) or (b0 == 0 and b1 > 0):
                    return False
        return True

    def events_at(self, v):
        return {("c", c) for c in self.clusters[v]} | {("s", s, b) for s, b in self.cn[v].items()}

    def node_of(self):
        """event -> node, for the events in the tree, as ('c', c) and ('s', s, state)."""
        return {e: v for v in range(1, self.k) for e in self.events_at(v)}

    def key(self):
        """Canonical form: the set of (events of a node, events of its parent)."""
        lab = [frozenset(self.events_at(v)) for v in range(self.k)]
        lab[0] = frozenset([("root",)])
        return frozenset((lab[v], lab[self.parent[v]]) for v in range(1, self.k))

    # ------------------------------------------------------------------ edits

    def has_slot(self, v, event):
        return v > 0 and (event[0] == "c" or event[1] not in self.cn[v])

    def place(self, event, v):
        """Copy of the tree with `event` added on the edge into v (None if v has no slot for it)."""
        if not self.has_slot(v, event):
            return None
        t = self.copy()
        if event[0] == "c":
            t.clusters[v].add(event[1])
        else:
            t.cn[v][event[1]] = event[2]
        return t

    def expand(self, event, u, num_segments, new_clone=True, existing=True):
        """Trees with `event` added: on a new clone (a leaf, or splitting an edge) and/or on an
        existing edge. Returns (tree, initial usages) pairs; a new clone gets index k."""
        out = []
        if new_clone:
            ch = self.children()
            for p in range(self.k):
                t = Tree(self.parent + [p], self.clusters + [set()], self.cn + [{}]).place(event, self.k)
                uu = np.append(u, u[p] / 2)
                uu[p] /= 2
                out.append((t, uu))
                for v in ch[p]:
                    t = Tree(self.parent + [p], self.clusters + [set()], self.cn + [{}])
                    t.parent[v] = self.k
                    out.append((t.place(event, self.k), np.append(u, 0.0)))
        if existing:
            for v in range(1, self.k):
                t = self.place(event, v)
                if t is not None:
                    out.append((t, np.array(u)))
        return [(t, uu) for t, uu in out if t is not None and t.valid(num_segments)]

    def remove(self, event, u):
        """(tree without `event`, usages, index of the contracted node or None). A node left
        without events is contracted into its parent, which takes its usage."""
        node = self.node_of()
        u = np.array(u)
        if event not in node:
            return self.copy(), u, None
        v = node[event]
        t = self.copy()
        if event[0] == "c":
            t.clusters[v].discard(event[1])
        else:
            del t.cn[v][event[1]]
        if t.clusters[v] or t.cn[v]:
            return t, u, None
        p = t.parent[v]
        u[p] += u[v]
        keep = [w for w in range(t.k) if w != v]
        new = {w: j for j, w in enumerate(keep)}
        par = [-1 if t.parent[w] < 0 else new[p if t.parent[w] == v else t.parent[w]] for w in keep]
        return Tree(par, [t.clusters[w] for w in keep], [t.cn[w] for w in keep]), u[keep], v

    # ------------------------------------------------------- SNV copy counts

    def snv_options(self, w, s, num_segments, cap=256):
        """Copy-count vectors (over nodes) of an SNV that arises on the edge into w, on segment s."""
        return snv_options(self.states(num_segments)[s], self.children(), self.parent, w, s in self.cn[w], cap)


def _counts(xp, xc, c):
    """Mutated copies in a child with xc copies of the allele, given c of xp in its parent."""
    if c == 0 or xp == xc:
        return [c]
    if c == xp:
        return [xc]
    if xc > xp:
        return list(range(c, c + xc - xp + 1))
    return list(range(max(0, xc - (xp - c)), min(c, xc) + 1))


def snv_options(st, ch, parent, w, event_at_w, cap=256):
    """Copy-count vectors (over nodes) of an SNV that arises on the edge into w.

    st[v] is the state of every node on the segment of the SNV, ch the children lists, and
    event_at_w whether the edge into w changes the segment. The SNV hits one copy of one
    allele. If the edge into w also changes the segment, the SNV may arise before or after
    that change. Below w, an amplification of the mutated allele adds 0..delta mutated copies
    (all of them if every copy is mutated), and a deletion keeps any feasible number of
    mutated copies. A count that reaches 0 stays 0.
    """
    k = len(parent)

    def below(h, x, c):
        """Assignments [(node, count), ...] for the strict subtree of x, given count c at x."""
        per_child = []
        for y in ch[x]:
            opts = []
            for cy in _counts(st[x][h], st[y][h], c):
                for rest in below(h, y, cy):
                    opts.append([(y, cy)] + rest)
                    if len(opts) >= cap:
                        break
            per_child.append(opts)
        res = []
        for combo in product(*per_child):
            res.append([a for part in combo for a in part])
            if len(res) >= cap:
                break
        return res or [[]]

    vecs = set()
    p = parent[w]
    for h in (0, 1):
        starts = []
        if st[w][h] >= 1:
            starts.append(1)                                                   # after the change
        if event_at_w and st[p][h] >= 1:
            starts += [c for c in _counts(st[p][h], st[w][h], 1) if c >= 1]   # before it
        for cw in sorted(set(starts)):
            for rest in below(h, w, cw):
                vec = [0] * k
                vec[w] = cw
                for y, cy in rest:
                    vec[y] = cy
                vecs.add(tuple(vec))
    return np.array(sorted(vecs), dtype=int).reshape(-1, k)


def from_truth(truth):
    return Tree(truth["parent"], [{c} if c >= 0 else set() for c in truth["cluster"]], truth["cn"])
