"""Test of the placement step on the true topology.

    python -m sapling.cn.fixtopo <instance dir> ... [--restarts 5] [--seed 0]

Each instance is a clonesim directory with snvs.tsv and cn.tsv (from generatebulk). The
placement loop gets the true topology only. It starts once from the true usages and the
true event placement (the reference), and --restarts times from a fit to the copy-number
calls alone (--init calls) or from random usages and placement (--init random). A restart succeeds when its loss is at most the reference loss
plus the tolerance. Accuracy is the fraction of SNVs on an edge equivalent to their true edge,
and of called copy-number states on an edge equivalent to where they first appear in the
true tree (equivalent: the same clones with cells below it).
"""
import argparse
import json
import time

import numpy as np

from .data import load_bulk, load_truth
from .place import Placer


def accuracy(res, truth, placer):
    """Fraction of SNVs, and of called copy-number states, on an edge equivalent to their true
    edge. Two edges are equivalent when the same clones with cells (positive true usage in
    some sample) lie below them, so the data cannot tell them apart."""
    alive = truth["u"].sum(axis=0) > 0
    below = [frozenset(v for v in placer.subtree[w] if alive[v]) for w in range(placer.k)]
    snv = float(np.mean([below[a] == below[b] for a, b in zip(res["edge"], truth["arises"])]))
    hits = tot = 0
    for s, nodes in enumerate(res["nodes"]):
        for b, v in truth["event_node"][s].items():
            tot += 1
            hits += b in nodes and below[nodes[b]] == below[v]
    return snv, hits / max(tot, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("instances", nargs="+")
    ap.add_argument("--restarts", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--kappa", type=float, default=-1.0)
    ap.add_argument("--eps", type=float, default=1.0)
    ap.add_argument("--init", choices=["calls", "random"], default="calls",
                    help="start from a fit to the copy-number calls (default) or from random usages and placement")
    ap.add_argument("--json", type=str)
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    log = print if a.verbose else None

    out = []
    for path in a.instances:
        bulk = load_bulk(path)
        truth = load_truth(path, bulk)
        placer = Placer(bulk, truth["parent"], kappa=a.kappa, eps=a.eps)
        tol = 2 * a.eps * bulk.num_samples
        unplaced = sum(len(placer.events[s]) - len(truth["event_node"][s]) for s in range(bulk.num_segments))
        t0 = time.time()
        true_cl = None
        if bulk.cluster is not None:
            true_cl = np.array([truth["arises"][bulk.cluster == c][0] for c in range(bulk.cluster.max() + 1)])
        ref = placer.run(truth["u"].copy(), [dict(x) for x in truth["event_node"]], true_cl, log=log)
        ref_time = time.time() - t0
        sa, ea = accuracy(ref, truth, placer)
        name = path.rstrip("/").split("/")[-1]
        print(f"{name}: {placer.k} nodes, {len(bulk.seg)} SNVs, {bulk.num_segments} segments, "
              f"{sum(map(len, placer.events))} copy-number events ({unplaced} not in the true tree), {bulk.num_samples} samples")
        print(f"  from the truth:  loss {ref['loss']:.1f}, {ref['rounds']} rounds, SNVs on true edge {sa:.3f}, "
              f"events on true edge {ea:.3f}, {ref_time:.1f}s")
        rng = np.random.default_rng(a.seed)
        rows = []
        for r in range(a.restarts):
            t0 = time.time()
            if a.init == "calls":
                u0, nodes0, cl0 = placer.init_from_calls(rng, log=log)
            else:
                u0, nodes0, cl0 = placer.random_start(rng)
            res = placer.search(u0, nodes0, cl0, log=log) if bulk.cluster is not None else placer.run(u0, nodes0, cl0, log=log)
            dt = time.time() - t0
            sa, ea = accuracy(res, truth, placer)
            ok = res["loss"] <= ref["loss"] + tol
            rows.append(dict(loss=res["loss"], lower=res["lower"], rounds=res["rounds"], snv_acc=sa, event_acc=ea,
                             success=bool(ok), time=dt, swaps=res.get("swaps", 0)))
            print(f"  restart {r + 1}: loss {res['loss']:.1f} ({res['loss'] - ref['loss']:+.1f}), {res['rounds']} rounds, "
                  f"{res.get('swaps', 0)} swaps, SNVs on true edge {sa:.3f}, events on true edge {ea:.3f}, {'success' if ok else 'FAIL'}, {dt:.1f}s", flush=True)
        out.append(dict(instance=path, nodes=placer.k, ref_loss=ref["loss"], ref_rounds=ref["rounds"], restarts=rows))
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
