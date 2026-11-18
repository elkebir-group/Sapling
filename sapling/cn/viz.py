"""Drawings of Sapling-CN trees with graphviz (needs `dot` on the path).

    python -m sapling.cn.viz truth <instance dir> <output prefix> [--loss L]
    python -m sapling.cn.viz run <instance dir> <run json> <output prefix> [--desc TEXT]

`truth` draws the true tree with the true multiplicities and usages. `run` draws a tree
inferred by sapling.cn.grow (its --json output): green events have the same ancestral events
as in the true tree, red ones do not. Both write <output prefix>.dot and .png.

Each copy-number line on an edge that carries clusters reads
    s9: (1,1) [C6 12] →[C6 5] (2,2) [C6 3]
for the SNVs of each cluster of the edge that lie on the segment: before the change (after the
parent state), undetermined (on the arrow) and after the change (after the child state), as
computed by grow.order_counts. Clusters with a count of 0 are left out of a bracket, and an
empty bracket is left off.
"""
import argparse
import json
import subprocess

import numpy as np

from .data import load_bulk, load_truth
from .grow import all_events, order_counts, true_tree
from .run import per_event_match
from .tree import Tree


def _bracket(counts, j):
    """[C0 8, C2 1], leaving out clusters with a count of 0; empty if every count is 0."""
    items = [f"C{c} {n[j]}" for c, n in sorted(counts.items()) if n[j] > 0]
    return "[" + ", ".join(items) + "]" if items else ""


def cn_line(bulk, a0, B, s, counts):
    """The label of copy-number event s: a0 -> B, with the before/after counts if any."""
    name = f"s{bulk.segment_ids[s]}: "
    counts = counts or {}
    parts = [f"({a0[0]},{a0[1]})", _bracket(counts, 0), "→" + _bracket(counts, 1), f"({B[0]},{B[1]})", _bracket(counts, 2)]
    return name + " ".join(x for x in parts if x)


def draw(bulk, tree, u, M, title, out, node_word="clone", match=None, missing=None):
    """Writes out.dot and out.png. u: usages (samples x clones); M: multiplicities (SNVs x
    clones); match: event -> True/False to color events green/red (None: black); missing:
    labels of events not in the tree."""
    counts = order_counts(bulk, tree, M)
    unmatched = sum(n[3] for per in counts.values() for n in per.values())
    if unmatched:
        print(f"  {unmatched} SNVs have multiplicities that are not among their options (left out of the counts)")
    st = tree.states(bulk.num_segments)
    L = ['digraph T {', '  node [shape=box, style=rounded, fontname="Helvetica", fontsize=10];', '  edge [fontname="Helvetica", fontsize=9];',
         f'  labelloc=t; label="{title}";']
    for v in range(tree.k):
        lab = ("normal" if v == 0 else f"{node_word} {v}") + "\\nu = " + " / ".join(f"{x:.3f}" for x in u[:, v])
        style = ', style="rounded,dashed", color=gray, fontcolor=gray' if u[:, v].max() < 1e-3 else ""
        L.append(f'  {v} [label="{lab}"{style}];')
    col = (lambda e: "black") if match is None else (lambda e: "darkgreen" if match.get(e) else "red")
    for v in range(1, tree.k):
        parts = [f'<font color="{col(("c", c))}"><b>C{c} ({int((bulk.cluster == c).sum())} SNVs)</b></font>' for c in sorted(tree.clusters[v])]
        for s, B in sorted(tree.cn[v].items()):
            line = cn_line(bulk, st[s][tree.parent[v]], B, s, counts.get((v, s)))
            parts.append(f'<font color="{col(("s", s, B))}">{line}</font>')
        L.append(f'  {tree.parent[v]} -> {v} [label=<{"<br/>".join(parts) or "—"}>];')
    if missing:
        L.append(f'  missing [shape=plaintext, label="not placed: {", ".join(missing)}"];')
    L.append("}")
    open(out + ".dot", "w").write("\n".join(L))
    subprocess.run(["dot", "-Tpng", "-Gdpi=120", out + ".dot", "-o", out + ".png"], check=True)


def draw_truth(instance, out, loss=None):
    """The true tree with the true multiplicities and usages."""
    b = load_bulk(instance)
    t = load_truth(instance, b)
    title = instance.rstrip("/").split("/")[-1] + ", ground truth (true placement, multiplicities and usages)"
    if loss is not None:
        title += f": loss {loss:.1f}"
    draw(b, true_tree(b, t), t["u"], t["M"], title, out, node_word="node")


def draw_run(instance, run_json, out, desc="inferred by growing"):
    """A tree from the --json output of sapling.cn.grow, colored by agreement with the truth."""
    b = load_bulk(instance)
    t = load_truth(instance, b)
    r = json.load(open(run_json))
    tr = r["tree"]
    if "M" not in tr:
        raise SystemExit(f"{run_json} has no multiplicities (M); rerun sapling.cn.grow to get them")
    tree = Tree(tr["parent"], [set(x) for x in tr["clusters"]], [{int(s): tuple(v) for s, v in x.items()} for x in tr["cn"]])
    events = all_events(b)
    match = per_event_match(tree, true_tree(b, t), events)
    missing = [f"s{b.segment_ids[e[1]]} {e[2]}" for e in events if e[0] == "s" and e not in tree.node_of()]
    title = (f"{instance.rstrip('/').split('/')[-1]}, {desc}: loss {r['loss']:.1f}, ground truth {r['gt_loss']:.1f} "
             "(green: same ancestral events as in the truth, red: different)")
    draw(b, tree, np.array(tr["u"]), np.array(tr["M"]), title, out, match=match, missing=missing)


def main():
    ap = argparse.ArgumentParser(description="Draw Sapling-CN trees with graphviz")
    sub = ap.add_subparsers(dest="what", required=True)
    p = sub.add_parser("truth", help="the true tree of an instance")
    p.add_argument("instance")
    p.add_argument("out", help="output prefix")
    p.add_argument("--loss", type=float, help="loss to show in the title")
    p = sub.add_parser("run", help="a tree inferred by sapling.cn.grow")
    p.add_argument("instance")
    p.add_argument("json", help="the --json output of sapling.cn.grow")
    p.add_argument("out", help="output prefix")
    p.add_argument("--desc", default="inferred by growing", help="how the tree was inferred, for the title")
    a = ap.parse_args()
    if a.what == "truth":
        draw_truth(a.instance, a.out, a.loss)
    else:
        draw_run(a.instance, a.json, a.out, a.desc)


if __name__ == "__main__":
    main()
