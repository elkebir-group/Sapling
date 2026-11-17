"""Input data for Sapling-CN: SNV read counts, copy-number calls, and SNV clusters.

`load_clonesim` reads a clonesim instance (tree.tsv, node.tsv, proportions.tsv),
simulates reads, and keeps the ground truth for evaluation.
"""
import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

NORMAL = (1, 1)


@dataclass
class Data:
    v: np.ndarray            # variant reads per SNV
    d: np.ndarray            # total reads per SNV
    seg: np.ndarray          # segment index of each SNV, in 0..S-1
    mu: list                 # mu[s] = {state: fraction}, the copy-number calls
    chat: np.ndarray         # mean copy number per segment
    clusters: list           # clusters[c] = array of SNV indices
    groups: dict = field(default_factory=dict)  # (c, s) -> SNV indices of cluster c on segment s
    truth: dict = None
    weight: np.ndarray = None  # w_s, weight of the squared loss on the calls of segment s

    def __post_init__(self):
        if self.weight is None:
            self.weight = np.ones(len(self.mu))
        self.cluster_segments = []
        for c, idx in enumerate(self.clusters):
            segs = [int(s) for s in np.unique(self.seg[idx])]
            self.cluster_segments.append(segs)
            for s in segs:
                self.groups[(c, s)] = idx[self.seg[idx] == s]

    @property
    def num_segments(self):
        return len(self.mu)

    def cn_events(self, min_mu=1e-9):
        """Copy-number events (s, state): every non-normal state with a positive call."""
        return [(s, st) for s in range(self.num_segments) for st, m in sorted(self.mu[s].items())
                if st != NORMAL and m > min_mu]

    def cluster_prevalence(self, c):
        """Median cellular prevalence of cluster c, assuming one mutated copy."""
        idx = self.clusters[c]
        return float(np.median(np.minimum(1.0, self.v[idx] / np.maximum(self.d[idx], 1) * self.chat[self.seg[idx]])))


def load_clonesim(path, depth=100, seed=0):
    path = Path(path)
    rng = np.random.default_rng(seed)
    edges = [tuple(map(int, l.split())) for l in open(path / "tree.tsv") if l.strip()]
    par = {c: p for p, c in edges}
    nodes = sorted({x for e in edges for x in e})
    root = next(x for x in nodes if x not in par)
    ch = {x: [] for x in nodes}
    for p, c in edges:
        ch[p].append(c)
    order, stack = [], [root]
    while stack:
        x = stack.pop()
        order.append(x)
        stack.extend(sorted(ch[x], reverse=True))
    idx = {x: j for j, x in enumerate(order)}
    parent = [-1] + [idx[par[x]] for x in order[1:]]
    k = len(order)

    states, mult, snv_seg = {}, {}, {}
    with open(path / "node.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            j, s, m = idx[int(r["node"])], int(r["segment"]), int(r["m"])
            states.setdefault(s, [None] * k)[j] = (int(r["x"]), int(r["y"]))
            if m >= 0:
                mult.setdefault(m, np.zeros(k, dtype=int))[j] = int(r["xbar"]) + int(r["ybar"])
                snv_seg[m] = s
    seg_ids = sorted(states)
    sid = {s: j for j, s in enumerate(seg_ids)}
    states = [states[s] for s in seg_ids]
    snvs = [i for i in sorted(mult) if mult[i].any()]
    M = np.array([mult[i] for i in snvs])
    seg = np.array([sid[snv_seg[i]] for i in snvs])

    u = np.zeros(k)
    with open(path / "proportions.tsv") as fh:
        next(fh)
        for line in fh:
            a, b = line.split()
            u[idx[int(a)]] = float(b)
    u /= u.sum()

    mu = []
    for st in states:
        frac = {}
        for j, ab in enumerate(st):
            frac[ab] = frac.get(ab, 0.0) + u[j]
        mu.append(frac)
    chat = np.array([sum((a + b) * x for (a, b), x in frac.items()) for frac in mu])

    y = M @ u
    d = rng.poisson(depth * chat[seg] / 2)
    v = rng.binomial(d, np.clip(y / chat[seg], 0, 1))

    # clusters: SNVs grouped by the node where they arise
    intro = np.array([next(w for w in range(k) if M[i, w] > 0 and (parent[w] < 0 or M[i, parent[w]] == 0))
                      for i in range(len(snvs))])
    cl_nodes = sorted(set(intro.tolist()))
    clusters = [np.flatnonzero(intro == w) for w in cl_nodes]
    cn = [{s: states[s][w] for s in range(len(states)) if w > 0 and states[s][w] != states[s][parent[w]]}
          for w in range(k)]
    truth = dict(parent=parent, states=states, mults=M, u=u, labels=order,
                 cluster=[cl_nodes.index(w) if w in cl_nodes else -1 for w in range(k)], cn=cn)
    return Data(v=v, d=d, seg=seg, mu=mu, chat=chat, clusters=clusters, truth=truth)
