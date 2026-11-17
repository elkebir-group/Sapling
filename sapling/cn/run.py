"""Run Sapling-CN on a clonesim instance and compare with the true tree.

    python -m sapling.cn.run <instance dir> --clones k [--beam 25] [--workers 8]

`--clones true` uses the number of clones of the true tree (the normal clone included).
fastppm must be the build of the `copy-number` branch (it provides regress_copy_number).
"""
import argparse
import json
import time

import numpy as np

from .data import load_clonesim
from .dot import write_dot
from .score import score
from .search import Search, all_events, tied
from .tree import from_truth


def relations(t, events):
    """(e, f) -> whether e is at or above f, for the events of `events` that are in t."""
    node, par = t.node_of(), t.parent
    up = {}
    for e in events:
        if e in node:
            v, path = node[e], set()
            while v >= 0:
                path.add(v)
                v = par[v]
            up[e] = path
    return {(e, f): node[e] in up[f] for e in up for f in up if e != f}


def ancestor_agreement(t, truth, events):
    """Fraction of ordered pairs of events placed in both trees with the same relation, by pair type."""
    r1, r2 = relations(t, events), relations(truth, events)
    agree = {"cc": [0, 0], "ss": [0, 0], "cs": [0, 0]}
    for (e, f), a in r1.items():
        if (e, f) in r2:
            typ = "".join(sorted(e[0] + f[0]))
            agree[typ][0] += a == r2[(e, f)]
            agree[typ][1] += 1
    return {k: (a / n if n else float("nan")) for k, (a, n) in agree.items()}


def per_event_match(t, truth, events):
    """event -> whether its set of ancestral events is the same in both trees."""
    def ancestors(tr):
        node, par = tr.node_of(), tr.parent
        res = {}
        for e in events:
            if e in node:
                v, s = par[node[e]], set()
                while v > 0:
                    s |= {f for f in tr.events_at(v) if f in events}
                    v = par[v]
                res[e] = frozenset(s)
        return res
    a1, a2 = ancestors(t), ancestors(truth)
    return {e: a1[e] == a2.get(e) for e in a1}


def consensus(trees, truth, events):
    """Ordered event pairs on which all trees agree (among pairs placed in all of them), and how many are correct."""
    rels = [relations(t, events) for t, *_ in trees]
    true = relations(truth, events)
    pairs = [p for p in true if all(p in r for r in rels)]
    agreed = [p for p in pairs if all(r[p] == rels[0][p] for r in rels)]
    return dict(resolved=len(agreed) / max(len(pairs), 1),
                correct=float(np.mean([rels[0][p] == true[p] for p in agreed])) if agreed else float("nan"))


def fmt(ag):
    return " ".join(f"{k}={v:.3f}" for k, v in ag.items())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("instance")
    ap.add_argument("--clones", type=str, required=True,
                    help="number of clones, the normal clone included (or 'true' for the count of the true tree)")
    ap.add_argument("--depth", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--beam", type=int, default=25)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--kappa", type=float, default=-1.0, help="penalty scale of fastppm-CN (default: computed from each instance)")
    ap.add_argument("--eps", type=float, default=1.0)
    ap.add_argument("--cn-weight", type=str, default="reads",
                    help="w_s: 'reads' (total reads of the SNVs of segment s) or a number")
    ap.add_argument("--sweeps", type=int, default=3)
    ap.add_argument("--json", type=str)
    ap.add_argument("--dot", type=str, help="prefix for <prefix>.true.dot, <prefix>.found.dot and <prefix>.beam<r>.dot")
    a = ap.parse_args()

    data = load_clonesim(a.instance, depth=a.depth, seed=a.seed)
    if a.cn_weight == "reads":
        data.weight = np.array([float(data.d[data.seg == s].sum()) or 1.0 for s in range(data.num_segments)])
    else:
        data.weight = np.full(data.num_segments, float(a.cn_weight))
    events = all_events(data)
    truth = from_truth(data.truth)
    k = truth.k if a.clones == "true" else int(a.clones)
    kw = dict(kappa=a.kappa, eps=a.eps)
    true_fixed = score(data, truth, data.truth["u"], fixed=(np.arange(len(data.v)), data.truth["mults"]), **kw)
    true_free = score(data, truth, data.truth["u"], **kw)
    print(f"{a.instance}: {truth.k} clones, {len(data.clusters)} clusters, {len(events) - len(data.clusters)} "
          f"copy-number events; true tree loss {true_free['loss']:.1f} "
          f"({true_fixed['loss']:.1f} with the true copy counts); searching for {k} clones, beam {a.beam}")

    t0 = time.time()
    srch = Search(data, beam=a.beam, workers=a.workers, kappa=a.kappa, eps=a.eps)
    beam = srch.grow(k)
    t_grow = time.time() - t0
    tree, u, loss, _ = beam[0]
    grown = ancestor_agreement(tree, truth, events)
    ties = tied(beam)
    cons = consensus(ties, truth, events)
    tree, u, loss = srch.refine(tree, u, loss, max_sweeps=a.sweeps)
    t_all = time.time() - t0
    srch.close()
    final = ancestor_agreement(tree, truth, events)
    missing = [e for e in events if e not in tree.node_of()]

    print(f"  grown: loss {beam[0][2]:.1f}, agreement {fmt(grown)}")
    print(f"  {len(ties)} of {len(beam)} beam trees tied with the best (lower bound <= best loss); "
          f"pairs on which all tied trees agree: {cons['resolved']:.3f} of pairs, {cons['correct']:.3f} of them correct")
    print(f"  refined: loss {loss:.1f} (true {true_free['loss']:.1f}), {tree.k} clones (true {truth.k}), "
          f"agreement {fmt(final)}, {len(missing)} copy-number events not placed, {srch.solves} solves, {t_all:.0f}s")

    out = dict(instance=a.instance, clones=k, clones_true=truth.k, beam=a.beam, loss_true=true_free["loss"],
               loss_true_fixed=true_fixed["loss"], loss_grown=beam[0][2], loss=loss, solves=srch.solves,
               time_grow=t_grow, time=t_all, agree_grown=grown, agree=final, tied=len(ties), consensus=cons,
               missing=[list(map(str, e)) for e in missing],
               tree=dict(parent=tree.parent, clusters=[sorted(x) for x in tree.clusters],
                         cn=[{s: list(b) for s, b in x.items()} for x in tree.cn], u=list(map(float, u))))
    if a.dot:
        open(a.dot + ".true.dot", "w").write(write_dot(truth, true_free["u"], data, f"true tree: loss {true_free['loss']:.1f}"))
        for r, (t2, u2, l2, _) in enumerate(beam[1:], start=2):
            open(a.dot + f".beam{r}.dot", "w").write(write_dot(
                t2, u2, data, f"beam tree {r} before refinement: loss {l2:.1f} (best before refinement {beam[0][2]:.1f})",
                per_event_match(t2, truth, events)))
        open(a.dot + ".found.dot", "w").write(write_dot(
            tree, u, data, f"inferred tree: loss {loss:.1f} (green: same ancestors as in the true tree)",
            per_event_match(tree, truth, events)))
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
