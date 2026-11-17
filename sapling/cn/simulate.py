"""Add an SNV cluster after every copy-number change of a clonesim instance.

    python -m sapling.cn.simulate <clonesim dir> <output dir> [--seed 0] [--size N]

Every edge that changes the copy number of some segment but has no SNV arising on it gets
a new cluster of SNVs (by default the median size of the existing clusters). Each new SNV
lies on a uniformly random segment, arises after the change on its edge on one copy of an
allele (chosen in proportion to its copies), and is propagated down the tree with
clonesim's rules: if every copy of the mutated allele carries the SNV, so do all copies in
the child; an amplification of the mutated allele adds every new copy to the mutated
ones; a deletion keeps a uniformly random feasible number of mutated copies. Topology,
copy-number states and proportions are copied unchanged; the files keep clonesim's format.
"""
import argparse
import csv
import shutil
from pathlib import Path

import numpy as np


def _child_count(xp, xc, c, rng):
    if c == 0 or xp == xc:
        return c
    if c == xp:
        return xc
    if xc > xp:
        return c + (xc - xp)
    lo, hi = max(0, xc - (xp - c)), min(c, xc)
    return int(rng.integers(lo, hi + 1))


def add_clusters(src, dst, seed=0, size=None):
    src, dst = Path(src), Path(dst)
    rng = np.random.default_rng(seed)
    edges = [tuple(map(int, l.split())) for l in open(src / "tree.tsv") if l.strip()]
    par = {c: p for p, c in edges}
    nodes = sorted({x for e in edges for x in e})
    root = next(x for x in nodes if x not in par)
    ch = {x: [] for x in nodes}
    for p, c in edges:
        ch[p].append(c)

    state, snv_seg, snv_cl, cnt = {}, {}, {}, {}
    with open(src / "node.tsv", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v, s, m = int(r["node"]), int(r["segment"]), int(r["m"])
            state[(v, s)] = (int(r["x"]), int(r["y"]))
            if m >= 0:
                snv_seg[m], snv_cl[m] = s, int(r["clusterIDofMutation"])
                cnt[(v, m)] = (int(r["xbar"]), int(r["ybar"]))
    segments = sorted({s for _, s in state})
    snvs = sorted(snv_seg)

    def arises_at(m):
        return next(v for v in nodes if sum(cnt[(v, m)]) > 0 and (v == root or sum(cnt[(par[v], m)]) == 0))

    origin = {arises_at(m) for m in snvs}
    sizes = [sum(1 for m in snvs if snv_cl[m] == c) for c in sorted(set(snv_cl.values()))]
    size = size or int(np.median(sizes))
    cn_edges = [v for v in nodes if v != root and v not in origin
                and any(state[(v, s)] != state[(par[v], s)] for s in segments)]

    next_m, next_c = max(snvs) + 1, max(snv_cl.values()) + 1
    added = []
    for w in cn_edges:
        for _ in range(size):
            m, s = next_m, int(rng.choice(segments))
            next_m += 1
            x, y = state[(w, s)]
            h = 0 if rng.random() < x / (x + y) else 1  # assumes x + y > 0 on the edge into w
            for v in nodes:
                cnt[(v, m)] = (0, 0)
            cnt[(w, m)] = (1, 0) if h == 0 else (0, 1)
            stack = [w]
            while stack:
                u = stack.pop()
                for v in ch[u]:
                    cu, cv = cnt[(u, m)][h], _child_count(state[(u, s)][h], state[(v, s)][h], cnt[(u, m)][h], rng)
                    cnt[(v, m)] = (cv, 0) if h == 0 else (0, cv)
                    stack.append(v)
            snv_seg[m], snv_cl[m] = s, next_c
            snvs.append(m)
        added.append((w, next_c, size))
        next_c += 1

    dst.mkdir(parents=True, exist_ok=True)
    shutil.copy(src / "tree.tsv", dst / "tree.tsv")
    shutil.copy(src / "proportions.tsv", dst / "proportions.tsv")
    by_seg = {s: [m for m in snvs if snv_seg[m] == s] for s in segments}
    with open(dst / "node.tsv", "w") as out:
        out.write("node\tsegment\tx\ty\tm\txbar\tybar\tclusterIDofMutation\n")
        for v in nodes:
            for s in segments:
                x, y = state[(v, s)]
                for m in by_seg[s] or [-1]:
                    xb, yb = cnt[(v, m)] if m >= 0 else (-1, -1)
                    out.write(f"{v}\t{s}\t{x}\t{y}\t{m}\t{xb}\t{yb}\t{snv_cl[m] if m >= 0 else -1}\n")
    return added


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--size", type=int, help="SNVs per new cluster (default: median of the existing clusters)")
    a = ap.parse_args()
    for w, c, n in add_clusters(a.src, a.dst, a.seed, a.size):
        print(f"clone {w}: new cluster {c} with {n} SNVs")


if __name__ == "__main__":
    main()
