"""Input data for Sapling-CN: SNV read counts and copy-number calls, in one or more samples.

`load_bulk` reads the input: snvs.tsv (snv, sample, segment, variant, total), cn.tsv
(segment, sample, x, y, fraction), and optionally clusters.tsv (snv, cluster), as written
by clonesim's generatebulk. The SNVs of a cluster are introduced on the same edge.
`load_truth` reads the true tree of a clonesim instance (tree.tsv, node.tsv,
proportions.tsv), for evaluation only.
`load_clonesim` is the input of the single-sample prototype (run.py, search.py), which
also takes SNV clusters from the truth.
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


@dataclass
class Bulk:
    """Read counts and copy-number calls of n SNVs on S segments in P samples."""
    v: np.ndarray        # (n, P) variant reads
    d: np.ndarray        # (n, P) total reads
    seg: np.ndarray      # (n,) segment index of each SNV, in 0..S-1
    mu: list             # mu[p][s] = {state: fraction}, the copy-number calls
    snv_ids: np.ndarray  # (n,) SNV ids of the input
    segment_ids: np.ndarray  # (S,) segment ids of the input
    cluster: np.ndarray = None  # (n,) cluster index of each SNV, in 0..C-1 (None without clusters.tsv)

    def __post_init__(self):
        S = len(self.segment_ids)
        self.chat = np.array([[sum((a + b) * x for (a, b), x in self.mu[p][s].items()) for s in range(S)]
                              for p in range(self.num_samples)])
        self.snvs_of = [np.flatnonzero(self.seg == s) for s in range(S)]
        # w_{p,s}: the total reads of the SNVs of segment s in sample p (1 if it has none)
        self.weight = np.array([[float(self.d[self.snvs_of[s], p].sum()) or 1.0 for s in range(S)]
                                for p in range(self.num_samples)])

    @property
    def num_samples(self):
        return self.v.shape[1]

    @property
    def num_segments(self):
        return len(self.segment_ids)

    def cn_events(self, s, min_mu=1e-9):
        """Copy-number states of segment s that need an event: not normal, called in some sample."""
        return sorted({b for p in range(self.num_samples) for b, m in self.mu[p][s].items() if b != NORMAL and m > min_mu})


def load_bulk(path):
    path = Path(path)
    snv, smp, sg, var, tot = [], [], [], [], []
    with open(path / "snvs.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            snv.append(int(r["snv"])); smp.append(int(r["sample"])); sg.append(int(r["segment"]))
            var.append(int(r["variant"])); tot.append(int(r["total"]))
    calls = {}
    with open(path / "cn.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (int(r["sample"]), int(r["segment"]))
            calls.setdefault(key, {})[(int(r["x"]), int(r["y"]))] = float(r["fraction"])
    samples = sorted({p for p, _ in calls} | set(smp))
    segment_ids = np.array(sorted({s for _, s in calls} | set(sg)))
    snv_ids = np.array(sorted(set(snv)))
    pi = {p: j for j, p in enumerate(samples)}
    si = {s: j for j, s in enumerate(segment_ids.tolist())}
    ni = {m: j for j, m in enumerate(snv_ids.tolist())}
    n, P = len(snv_ids), len(samples)
    v, d = np.full((n, P), -1), np.full((n, P), -1)
    seg = np.full(n, -1)
    for m, p, s, a, t in zip(snv, smp, sg, var, tot):
        j = ni[m]
        if seg[j] >= 0 and seg[j] != si[s]:
            raise ValueError(f"SNV {m} lies on two segments")
        seg[j], v[j, pi[p]], d[j, pi[p]] = si[s], a, t
    if (v < 0).any():
        raise ValueError("snvs.tsv must have every SNV in every sample")
    mu = [[calls.get((p, s), {NORMAL: 1.0}) for s in segment_ids.tolist()] for p in samples]
    cluster = None
    if (path / "clusters.tsv").exists():
        cl = {}
        with open(path / "clusters.tsv", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                cl[int(r["snv"])] = int(r["cluster"])
        if set(cl) != set(ni):
            raise ValueError("clusters.tsv must have every SNV of snvs.tsv")
        ids = sorted(set(cl.values()))
        ci = {c: j for j, c in enumerate(ids)}
        cluster = np.array([ci[cl[m]] for m in snv_ids.tolist()])
    return Bulk(v=v, d=d, seg=seg, mu=mu, snv_ids=snv_ids, segment_ids=segment_ids, cluster=cluster)


def _read_clonesim(path):
    """Topology (node 0 is the root, nodes in preorder), states, copies of mutated alleles,
    and proportions of a clonesim instance, with its node, segment and SNV ids."""
    path = Path(path)
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

    rows = [l.split() for l in open(path / "proportions.tsv").read().splitlines()[1:] if l.strip()]
    u = np.zeros((len(rows[0]) - 1, k))
    for r in rows:
        u[:, idx[int(r[0])]] = [float(x) for x in r[1:]]
    u /= u.sum(axis=1, keepdims=True)
    return dict(parent=parent, labels=order, states=states, mult=mult, snv_seg=snv_seg, u=u)


def load_truth(path, bulk):
    """The true tree of a clonesim instance, indexed like `bulk`: parent (node 0 is the
    normal clone), states[s][v], copies M (n x k) of every SNV, usages u (P x k), the node
    where each SNV arises, and the node where each called copy-number state first appears."""
    t = _read_clonesim(path)
    parent, k = t["parent"], len(t["parent"])
    states = [t["states"][s] for s in bulk.segment_ids.tolist()]
    M = np.array([t["mult"][m] for m in bulk.snv_ids.tolist()])
    arises = np.array([next(w for w in range(k) if M[i, w] > 0 and (parent[w] < 0 or M[i, parent[w]] == 0))
                       for i in range(len(M))])
    event_node = []
    for s in range(bulk.num_segments):
        first = {}
        for w in range(1, k):
            b = states[s][w]
            if b != states[s][parent[w]] and b not in first:
                first[b] = w
        event_node.append({b: first[b] for b in bulk.cn_events(s) if b in first})
    return dict(parent=parent, labels=t["labels"], states=states, M=M, u=t["u"], arises=arises,
                event_node=event_node)


def load_clonesim(path, sample=0):
    """One sample of an instance with snvs.tsv and cn.tsv (from generatebulk), as Data,
    with SNV clusters and the truth from the clonesim files."""
    bulk = load_bulk(path)
    t = _read_clonesim(path)
    parent, labels = t["parent"], t["labels"]
    k = len(parent)
    states = [t["states"][s] for s in bulk.segment_ids.tolist()]
    M = np.array([t["mult"][m] for m in bulk.snv_ids.tolist()])
    keep = np.flatnonzero(M.any(axis=1))
    M, seg = M[keep], bulk.seg[keep]
    v, d = bulk.v[keep, sample], bulk.d[keep, sample]
    mu = bulk.mu[sample]
    chat = bulk.chat[sample]
    u = t["u"][sample]

    # clusters: SNVs grouped by the node where they arise
    intro = np.array([next(w for w in range(k) if M[i, w] > 0 and (parent[w] < 0 or M[i, parent[w]] == 0))
                      for i in range(len(M))])
    cl_nodes = sorted(set(intro.tolist()))
    clusters = [np.flatnonzero(intro == w) for w in cl_nodes]
    cn = [{s: states[s][w] for s in range(len(states)) if w > 0 and states[s][w] != states[s][parent[w]]}
          for w in range(k)]
    truth = dict(parent=parent, states=states, mults=M, u=u, labels=labels,
                 cluster=[cl_nodes.index(w) if w in cl_nodes else -1 for w in range(k)], cn=cn)
    return Data(v=v, d=d, seg=seg, mu=mu, chat=chat, clusters=clusters, truth=truth)
