# Sapling-CN: notes for continuing the prototype

Sapling-CN infers a clone tree from bulk read counts together with copy-number calls.
It grows the tree one clone at a time, as Sapling does, and scores every candidate tree with fastppm-CN.
This note explains how to run it, what the code does, what we have seen so far, and what to do next.

## Setup

Sapling-CN needs the Python module of fastppm built from the `copy-number` branch of `elkebir-group/fastppm`.
That branch adds `fastppm.regress_copy_number`.

```bash
git clone -b copy-number git@github.com:elkebir-group/fastppm.git fastppm-cn
cd fastppm-cn && mkdir build-std && cd build-std
cmake .. -DCMAKE_BUILD_TYPE=Release -DSPDLOG_USE_STD_FORMAT=ON
make -j8 fastppm
```

The flag `-DSPDLOG_USE_STD_FORMAT=ON` is needed with recent Apple clang, which fails to compile the copy of fmt bundled with spdlog.
The module lands in `build-std/src/fastppm.cpython-<version>-<platform>.so`.
Put that directory on `PYTHONPATH` and use the same Python version it was built for.
Sapling-CN itself needs only numpy.

## Running

From the root of this repository:

```bash
PYTHONPATH=<fastppm-cn>/build-std/src:. python -m sapling.cn.run <clonesim instance> --clones true \
    --dot out/101 --json out/101.json
dot -Tpng out/101.found.dot -o out/101.found.png
```

An instance is a clonesim output directory with `tree.tsv`, `node.tsv` and `proportions.tsv`.
The 16 test trees are in `fastppm-dollo/lean/implementation/seeded-test-trees/`, and tree 101 is the one we studied.
The main options are the following.

| option | meaning |
|---|---|
| `--clones k` | number of clones, the normal clone included (required); `true` takes the count of the true tree |
| `--beam` | beam width (default 25) |
| `--workers` | processes that score candidate trees in parallel (default 8) |
| `--depth`, `--seed` | read depth for a diploid segment, and the seed of the read simulation |
| `--cn-weight` | weight w_s of the copy-number calls of segment s (default `reads`, see below) |
| `--kappa`, `--eps` | penalty scale and duality-gap tolerance of fastppm-CN (default: kappa computed from each instance, and a gap of 1) |
| `--sweeps` | passes of event moves after growing (default 3) |
| `--dot`, `--json` | output prefix for Graphviz files, and a JSON summary |

The run prints the loss of the true tree and of the inferred tree, and the agreement between them.
Agreement is the fraction of ordered pairs of events whose ancestor relation is the same in both trees, split into cluster pairs (`cc`), copy-number pairs (`ss`) and mixed pairs (`cs`).
It also reports which trees of the final beam are tied with the best.
A tree is tied when its lower bound is at most the best loss, so the solver cannot rank the two.
The DOT files show the true tree, the inferred tree, and the other trees of the final beam.
Events in green have the same ancestor events as in the true tree.

To simulate an instance with an SNV cluster after every copy-number change, run

```bash
python -m sapling.cn.simulate <clonesim instance> <output dir> [--size N] [--seed 0]
```

It keeps the topology, copy-number states and proportions, and adds a cluster on every edge that changes copy number but has no SNVs.
New clusters have the median size of the existing ones unless `--size` is given, so the total number of SNVs grows.
Several copy-number changes on one edge share one new cluster after them.

## Model

The input is the variant and total reads of every SNV, the segment of every SNV, and the copy-number calls of every segment.
The calls are the fraction mu_{s,B} of cells in each allele-specific state B, and the mean copy number of the segment.
SNV clusters are taken from the simulation for now, grouping SNVs by the clone where they arise.

The tree has the normal clone as its root, in state (1,1) on every segment.
Every other edge carries SNV clusters (any number) and copy-number events.
A copy-number event (s, B) sets segment s to state B, with at most one event per segment on an edge.
The state of a clone on s is set by the nearest event on s at or above it.
Each state appears on one edge only (no homoplasy), and an allele with zero copies cannot come back.
A copy-number event may be missing from a tree, and its call then counts as unexplained.

Each SNV is gained once, on the edge of its cluster, on one copy of one allele.
If that edge also changes the segment, the SNV may arise before or after the change.
Below it, an amplification of the mutated allele adds between 0 and all of the new copies, and a deletion keeps any feasible number of mutated copies.
So an SNV can be lost any number of times, but only where its allele is deleted, and it never comes back.
`Tree.snv_options` lists the resulting copy-count vectors.

## Scoring a tree

The loss is the binomial negative log-likelihood of all SNVs plus w_s times the squared error of the state fractions of every segment.
By default w_s is the total read depth of the SNVs on segment s.
Without this weight the calls count for almost nothing next to a likelihood over hundreds of thousands of reads, and copy-number events are placed badly.

Given the usages, the copy counts of each SNV are chosen separately, by trying every option and keeping the best.
Given the copy counts, fastppm-CN solves for the usages.
`score.score` alternates the two, at most 3 rounds, starting from usages inherited from the parent tree.
This is coordinate descent, so the result is a local optimum.

## Search

1. Start from the normal clone and one tumor clone that carries every SNV cluster.
2. To go from k to k+1 clones, add one clone as a new leaf or by splitting an edge, and give it one event.
   The event is moved from another edge, or it is a copy-number event not yet in the tree.
   All such trees are scored and the best `--beam` are kept.
3. On each kept tree, the absorb step moves events onto the new clone and adds missing copy-number events to any edge.
   It scores every single move, applies the improving ones greedily, and repeats until nothing improves.
4. Stop at the requested number of clones.
   Then move one event at a time to its best edge, never exceeding that number of clones.

Every tree within a step has the same number of clones and contains every SNV, so their losses compare directly and no penalty on the number of clones is needed.

## Results on tree 101

Both runs use `--clones true` (7 clones), a beam of 25, and depth 100 with one sample.

| instance | loss found | loss of true tree | cluster pairs | copy-number pairs | mixed pairs | time |
|---|---|---|---|---|---|---|
| tree 101 | 142935.6 | 142944.9 | 1.00 | 0.89 | 0.94 | 29 s |
| tree 101 with SNVs after every copy-number change | 192742.5 | 192740.2 | 0.93 | 0.96 | 0.93 | 144 s |

On tree 101 the inferred tree fits slightly better than the true tree.
Its trunk and the edge of C1 match the true tree.
The differences involve clones with tiny usage, or placements that the sum condition allows with one sample.
Segment 8 has the only nested copy-number change, (1,1) to (2,1) to (2,2), and the inferred tree does not recover it.

With SNVs after every copy-number change, almost every edge is recovered.
The true tree is not in the final beam, though.
Its loss is below the lower bound of the best beam tree, so the search missed a better tree here.
The remaining errors are C5 placed beside C1 instead of below C3, and the order of the two changes on segment 8 reversed.
Fixing them needs two events to move together, which single-event moves cannot do.

## Known problems

- With one sample, many trees fit almost equally well, for example sibling placements allowed by the sum condition.
  The tie report shows these.
- Absorb and refinement move one event at a time and need a full solve for every trial.
  They miss improvements that need several events to move together.
- Each SNV picks its own best copy counts.
  This freedom is not penalized, and trees with more copy-number events below a cluster gain options for its SNVs.
  Sharing one choice per cluster and segment would remove it.
- The copy-count options are capped at 256 per cluster and segment, which could matter in deep subtrees.
- SNV clusters come from the simulation.
- The beam keeps the 25 best trees at every step, and the trees on the way to the true tree can rank low before the tree is complete.
- fastppm-CN computes its penalty scale kappa from each instance, as 1 / sqrt(lambda_min lambda_max) of a k by k matrix, and polishes the usages when the gap stalls.
  Rules that adapt kappa while solving did not help (see `fastppm-dollo/paper/adaptive-penalty.md`).
- Every solve starts cold.

## Next steps

1. Optimize the placement of events on a fixed topology by alternating, as k-means does.
   With the usages fixed, each SNV cluster picks its best edge independently.
   Each segment then picks its best assignment of states to edges, which is a small enumeration.
   One fastppm-CN solve then updates the usages, and the loop repeats until nothing changes.
   This costs one solve per round instead of one per trial move, and it moves all events at once.
   A first test is to fix the true topology of tree 101, scramble the events, and check that the loop recovers the true placement.
   It would then replace absorb inside the search.
2. Run all 16 test trees, with and without SNVs after every copy-number change.
3. Infer the SNV clusters instead of taking them from the simulation.
4. Use several samples, which fastppm-CN does not support yet.
5. Warm-start fastppm-CN from the parent tree, to make wider beams affordable.
