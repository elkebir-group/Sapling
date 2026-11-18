"""Score a clone tree with fastppm-CN.

The loss is the fastppm-CN objective over the SNVs of the placed clusters and the
copy-number calls of every segment. The copy counts of each SNV are chosen per SNV
(given the usages they are separable), alternating with fastppm-CN solves.
"""
import numpy as np

import fastppm

from .data import NORMAL


def _nll(v, d, p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -(v[:, None] * np.log(p) + (d - v)[:, None] * np.log1p(-p))


def assign(data, tree, u, options):
    """Best copy-count vector for every SNV of every placed cluster, given usages u."""
    rows, snvs = [], []
    for (c, s), O in options.items():
        idx = data.groups[(c, s)]
        p = (O @ u) / data.chat[s]
        best = np.argmin(_nll(data.v[idx], data.d[idx], p[None, :]), axis=1)
        rows.append(O[best])
        snvs.append(idx)
    if not rows:
        return np.zeros(0, dtype=int), np.zeros((0, tree.k), dtype=int)
    return np.concatenate(snvs), np.concatenate(rows)


def regress(data, tree, snvs, M, kappa=-1.0, eps=1.0, max_iterations=20000):
    S = data.num_segments
    st = tree.states(S)
    res = fastppm.regress_copy_number(
        parent=tree.parent,
        states=[list(st[s]) for s in range(S)],
        multiplicity=np.ascontiguousarray(M, dtype=np.int32),
        segment=np.ascontiguousarray(data.seg[snvs], dtype=np.int32),
        variant=np.ascontiguousarray(data.v[snvs], dtype=np.int32),
        total=np.ascontiguousarray(data.d[snvs], dtype=np.int32),
        state_fractions=[[(a, b, float(m)) for (a, b), m in data.mu[s].items()] for s in range(S)],
        mean_copy_number=data.chat.tolist(), segment_weight=data.weight.tolist(),
        kappa=kappa, eps=eps, max_iterations=max_iterations)
    # calls of states that the tree does not have yet: their squared loss is mu^2
    missing = sum(data.weight[s] * m * m for s in range(S) for b, m in data.mu[s].items() if b not in set(st[s]))
    return res["loss"] + missing, res["lower_bound"] + missing, np.array(res["usages"]), res["converged"]


def score(data, tree, u0, rounds=3, fixed=None, **kw):
    """Returns dict(loss, lower, u, snvs, M). With `fixed` = (snvs, M), copy counts are not chosen."""
    S = data.num_segments
    if fixed is not None:
        snvs, M = fixed
        L, lb, u, ok = regress(data, tree, snvs, M, **kw)
        return dict(loss=L, lower=lb, u=u, snvs=snvs, M=M, converged=ok)
    options = {}
    for w in range(1, tree.k):
        for c in tree.clusters[w]:
            for s in data.cluster_segments[c]:
                O = tree.snv_options(w, s, S)
                if len(O) == 0:
                    return dict(loss=np.inf, lower=np.inf, u=u0, snvs=None, M=None, converged=False)
                options[(c, s)] = O
    u = np.asarray(u0, dtype=float)
    prev = None
    for _ in range(rounds):
        snvs, M = assign(data, tree, u, options)
        if prev is not None and np.array_equal(M, prev):
            break
        L, lb, u, ok = regress(data, tree, snvs, M, **kw)
        prev = M
    return dict(loss=L, lower=lb, u=u, snvs=snvs, M=M, converged=ok)
