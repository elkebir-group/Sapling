# Sapling

Sapling is a method that infers a small set of backbone trees on a smaller subset of mutations that collectively summarize the entire set of possible phylogenies.
Sapling can also grow given backbone trees into full phylogeny. Finally, Sapling can directly output full phylogenies up to a specified fraction `1-rho` away from optimality.

![overview](docs/overview.png)

The data repository containing all simulation and real data and results used in the paper can be found [here](https://github.com/elkebir-group/Sapling-data).

## Installation

To install using `conda`, execute the following command:

```bash
conda install -c bioconda sapling
```

This will also install `fastppm`.

## Dependencies

Sapling requires the following python packages:
* numpy
* pandas
* tqdm
* one of [cvxopt](https://cvxopt.org)/[gurobipy](https://pypi.org/project/gurobipy/)/[fastppm](https://github.com/elkebir-group/fastppm) (recommended)

If you choose to use [fastppm](https://github.com/elkebir-group/fastppm), please make sure that the corresponding python library is installed in `$PYTHONPATH`. This includes the current directory. For example:

```bash
ls -alFh fastppm*
lrwxr-xr-x 1 melkebir 72 Jan 28 14:33 fastppm.cpython-311-darwin.so -> /Users/melkebir/Projects/fastppm/build/src/fastppm.cpython-311-darwin.so*
```

## I/O Format

Sapling takes a TSV (tab-separated values) file as input.
The first line includes the names of the columns. The following columns are required:

- `sample_index`: (0, 1, ..., m-1) for m samples.
- `mutation_index`: (0, 1, ..., n-1) for n mutations.
- `var`: The number of variant reads supporting the mutation.
- `depth`: The total read depth at the locus.
- `cluster_index` (optional): (0, 1, ..., k-1) for k mutation clusters. If not provided, cluster_index defaults to mutation_index.

[Here](example/example_input.tsv) is an example input file. 

The output is a TSV file.
The first line includes the names of the columns, including:

- `tree`: The index of the tree.
- `llh`: The log-likelihood value of the tree.
- `pi_i`: The parent index for mutation `i` (where `i` is the mutation index). A value of `-1` indicates that node `i` is a root node.
- `f_p_i`: The inferred frequency of mutation `i` in sample `p`.

[Here](example/example_output.tsv) is an example output file. 

## Usage instructions

```
      usage: sapling [-h] [-f F] [--init_trees INIT_TREES] [-o O] [--sep SEP] [-a RHO] [-t TAU] [-l ELL] [--big_expand] [-b BEAM_WIDTH] [-L {cvxopt,pLP,fastppm}] [--alt_roots] [-m] [--use_clusters]

      Sapling is an algorithm for summarizing and inferring tumor phylogenies from bulk DNA sequencing data

      options:
      -h, --help            show this help message and exit
      -f F                  Input filename with mutation read counts (default: STDIN)
      --init_trees INIT_TREES
                              Input filename with initial backbone trees to expand
      -o O                  Output filename to store trees and frequencies (default: STDOUT)
      --sep SEP             Input/output column separator (default: \t)
      -a RHO, --rho RHO     Rho parameter, minimum deviation allowed from max likelihood (default: 0.9, ignored when beam_width specified)
      -t TAU, --tau TAU     Tau parameter, maximum number of backbone trees (default: 5)
      -l ELL, --ell ELL     Ell parameter, minimum number of mutations (default: -1, unlimited)
      --big_expand          Use big expand (new mutations are anywhere, not just as leaves or splitting a single edge)
      -b BEAM_WIDTH, --beam_width BEAM_WIDTH
                              Maximum beam width (default: -1, limited only by --rho)
      -L {cvxopt,pLP,fastppm}
                              Regression method (default: fastppm)
      --alt_roots           Explore alternative root nodes
      -m, --poly_clonal_root
                              Allow poly clonal root node
      --use_clusters        Use provided clustering (taking median read depth and using average frequency for variant counts)
```

Example command:

```bash
sapling --tau 5 --rho 0.9 -f example/example_input.tsv > example/example_output.tsv
```

This will output up to `tau=5` backbone trees with a likelihood cut-off of `rho==0.9`. [The output](example/example_output.tsv) of the above command on the example input.

#### Infer full trees

To infer full trees use the following options:

```bash
sapling --tau -1 --ell -1 --beam_width 100 -f example/example_input.tsv > example/example_full_trees.tsv`
```

This will use a beam width of 100 to enumerate up to a 100 trees containg all mutations. [Here](example/example_full_trees) is the output.

#### Expand given backbone trees into full trees

Example command:

```bash
sapling --tau -1 --ell -1 --rho 0.9 --init_trees example/example_output.tsv -f example/example_input.tsv -o example/example_expand.tsv
```

This will expand the [given backbone trees](example/example_output.tsv) into full trees (no restrictions on `tau` and `ell`) that are a factor of `1-rho=1-0.9=0.1` away from optimality. [Here](example/example_expand.tsv) is the output of the above command.
