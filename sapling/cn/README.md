# Sapling-CN: notes for continuing the prototype

Sapling-CN infers a clone tree from bulk read counts of several samples together with copy-number calls.
It grows the tree one clone at a time, as Sapling does, and scores every candidate tree with fastppm-CN.
This note explains how to build and run it, what the code does, what we have seen so far, and what to do next.

## Setup

Sapling-CN needs fastppm-CN, the `copy-number` branch of `elkebir-group/fastppm`, with the warm-start option `initial_usages` of `CopyNumberADMM::Options`.
Its scoring core is a C++ extension built against the fastppm-CN sources.

```bash
git clone -b copy-number git@github.com:elkebir-group/fastppm.git fastppm-cn
cmake -S cpp -B build -DFASTPPM_CN_DIR=<fastppm-cn> -DPython_EXECUTABLE=<python>
cmake --build build -j8
```

This writes `sapling/cn/_core.cpython-<version>-<platform>.so`.
The build sets `SPDLOG_USE_STD_FORMAT`, since recent Apple clang fails to compile the copy of fmt bundled with spdlog.
`grow.py` also imports the Python module of fastppm-CN, for the reference scoring in Python.
Build it in the fastppm-cn checkout as below and put `build-std/src` on `PYTHONPATH`, with the same Python version.

```bash
cd fastppm-cn && mkdir build-std && cd build-std
cmake .. -DCMAKE_BUILD_TYPE=Release -DSPDLOG_USE_STD_FORMAT=ON
make -j8 fastppm
```

Python needs numpy, and graphviz (`dot`) draws the trees.

## Input

An instance is a directory with three files, as written by `generatebulk` of clonesim (branch `bounded-sampling`).

| file | columns |
|---|---|
| `snvs.tsv` | snv, sample, segment, variant reads, total reads |
| `cn.tsv` | segment, sample, x, y, fraction of cells in state (x, y) |
| `clusters.tsv` | snv, cluster |

A cluster is the set of SNVs introduced on one edge of the true tree.
For evaluation, the directory also holds the clonesim output (`tree.tsv`, `node.tsv`, `proportions.tsv`).
The test instances were made with

```bash
simulate -r -S <cnatrees> -purity .99 -minProp 0.01/l -threshold 0.2/max(k,l) -kk 2 -f -s <seed> -l {10,20} -k {20,40} -n {2000,4000} -m 3
generatebulk -in_dir <dir> -depth 100 -s <seed>
```

## Running

From the root of this repository:

```bash
PYTHONPATH=<fastppm-cn>/build-std/src:. python -m sapling.cn.grow <instance> --clones true --beam 10 --json out.json
python -m sapling.cn.viz truth <instance> out/truth
python -m sapling.cn.viz run <instance> out.json out/found
```

| option of `grow` | meaning |
|---|---|
| `--clones k` | number of clones, the normal clone included; `true` takes the count of the true tree |
| `--beam` | beam width (default 10) |
| `--workers` | threads of the scoring core (default 8) |
| `--screen` | the best `screen` × beam candidates of a step are scored fully (default 3) |
| `--kappa`, `--eps` | penalty scale and duality-gap tolerance of fastppm-CN (default: kappa computed from each instance, and a gap of 1) |
| `--json` | output: the tree, its usages and multiplicities, and the best tree of every step (`history`) |

The run prints the loss of the true tree (with the true multiplicities) and of the inferred tree, and the agreement between them.
Agreement is the fraction of ordered pairs of events whose ancestor relation is the same in both trees, over all pairs and split into cluster pairs (`cc`), copy-number pairs (`ss`) and mixed pairs (`cs`).

In a drawing, green events have the same ancestor events as in the true tree.
Each copy-number line on an edge with clusters reads `s9: (1,1) [C6 12] →[C6 5] (2,2) [C6 3]`.
The brackets count the SNVs of each cluster of the edge on that segment that arose before the change, that cannot be classified, and that arose after it.
Counts of 0 are left out.

`run.py`, `place.py` and `fixtopo.py` are earlier experiments: one sample with reads drawn from the clonesim proportions, and placement of events on a fixed topology.

## Model

The tree has the normal clone as its root, in state (1,1) on every segment.
Every other edge carries SNV clusters (any number) and copy-number events, and at least one event.
A copy-number event (s, B) sets segment s to state B, with at most one event per segment on an edge.
The state of a clone on s is set by the nearest event on s at or above it.
Each state appears on one edge only (no homoplasy), and an allele with zero copies cannot come back.
A copy-number event may be missing from a tree, and its call then counts as unexplained.

Each SNV is gained once, on the edge of its cluster, on one copy of one allele.
If that edge also changes the segment, the SNV may arise before or after the change.
Below it, an amplification of the mutated allele adds between 0 and all of the new copies, and a deletion keeps any feasible number of mutated copies.
`tree.snv_option_kinds` lists the resulting multiplicity vectors, each marked as coming from an SNV before or after the change.

## Scoring a tree

The loss sums over the samples the binomial negative log-likelihood of all SNVs plus w_{p,s} times the squared error of the state fractions of every segment.
w_{p,s} is the total read depth of the SNVs on segment s in sample p.
The denominator of each VAF is the mean copy number from the calls, as in fastppm-CN.

Given the usages, every SNV takes the allowed multiplicities that fit its reads in all samples best.
Given the multiplicities, fastppm-CN solves for the usages of each sample, warm-started from the current usages.
Scoring alternates the two, at most 3 rounds.
This is coordinate descent, so the result is a local optimum.

`grow.score` does this in Python, and the C++ core (`cpp/core.cpp`) does the same for a batch of trees on several threads.
The two agree to the last digit, and the search uses the core.

## Search

1. Start from the normal clone and one tumor clone that carries every SNV cluster.
2. To go from k to k+1 clones, add one clone as a new leaf or by splitting an edge, and give it one event.
   The event is moved from another edge, or it is a copy-number event not yet in the tree.
   All candidates are screened with one round of alternation and at most 10 iterations of fastppm-CN per sample.
   The best `screen` × beam are scored fully, and the best `--beam` are kept.
3. Each kept tree is polished.
   A polishing step tries every move that keeps the number of clones, and applies the best one.
   The moves reattach a clone with its subtree, move everything on a clone's edge to a new clone, swap everything on two edges, and move a single event.
   A move is accepted only when its loss is below the lower bound of the current tree, so it lowers the minimum for the chosen multiplicities.
4. Stop at the requested number of clones, and polish the best tree once more.

## Results

All instances have 3 samples at depth 100, and every run uses `--clones true`.
The loss of the true tree is given with refitted multiplicities, scored exactly as the search scores trees.

| instance | clones | beam | loss found − loss of true tree | agreement (all pairs) | time |
|---|---|---|---|---|---|
| l10s41 | 14 | 5 | −0.6 | 0.995 | 92 s |
| l10s45 | 16 | 10 | −11.5 | 0.932 | 160 s |

On both instances the inferred tree fits at least as well as the true tree, so the remaining differences come from the model and the data, not from the search.
On l10s45 the differences are the order of two-step changes ((1,1) to (1,2) to (2,2) inferred as (1,1) to (2,2) to (1,2)), a leaf cluster placed at the bottom of another chain, and a subtree attached one level off.
The gain of the inferred tree comes entirely from SNVs with several allowed multiplicities.
Misplacing C1 below C4 and C5 gives their SNVs on s19 an extra option, which they use to fit noise.

An l20 tree (28 clones) takes about 8 minutes at beam 5 with the earlier absorb step instead of polishing.

## Known problems

- Each SNV picks its own best multiplicities, with no prior.
  Trees that give SNVs more options gain from it.
  A prior over multiplicities, or summing over them, is the principled fix, but the prior is not defined yet.
  It follows from the copy-number history (an allele is hit in proportion to its copies, and later gains and losses act on random copies), except for the order of an SNV and a copy-number change on the same edge.
- The number of clones is given.
  The loss keeps falling as clones are added, so it cannot choose it.
  The marginal likelihood with the usages integrated out (Laplace approximation with a correction for usages at 0) peaks at 14 clones on l10s45, against 16 in the true tree.
- The copy-number term is a weighted squared error, not a likelihood, and its weight sets how much an extra clone is penalized in the marginal likelihood.
- The multiplicity options are capped at 256 per cluster and segment, which could matter in deep subtrees.
- SNV clusters come from the simulation.
  Placing single SNVs on the true tree puts only 55–93% on their true edge, so inferring clusters needs to pool SNVs.

## Next steps

1. Read how clonesim generates multiplicities, and define the prior over multiplicities from it.
2. Score the trees of a run (the true tree and the best tree of every step) with the marginal likelihood on all l10 instances, before changing the search.
3. Reduce the candidates per growth step for the l20 trees, for example by adding new clones only as leaves and letting polishing move them.
