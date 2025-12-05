import sys
import pandas as pd
import numpy as np
from . import bbt_solver
import argparse
from math import log
from math import ceil

from . import __version__

def parse_input_trees(filename, m, n, sep):
    """
    Parses a CSV file containing tree structures and their associated log-likelihood values.

    The CSV file is expected to have columns representing tree edges in the format 'pi_X',
    where 'X' is the target node, and the value in the column represents the source node.
    Additionally, the file must contain a 'llh' column representing the log-likelihood
    of each tree.

    Args:
        filename (str): Path to the CSV file containing the tree data.
        sep (str): Delimiter used in the CSV file (e.g., ',' for comma-separated values).
        m (int): Number of samples
        n (int): Number of total mutations (in the actual input read count data)

    Returns:
        list: A list of tuples, where each tuple contains:
            - A list of edges, where each edge is represented as a tuple (source, target).
            - A tuple of the frequency matrix (F) and the log-likelihood (llh) value associated with the tree.

    Raises:
        SystemExit: If the input file does not contain a 'llh' column, the function will
                   terminate the program with an error message.

    Example:
        Given a CSV file with the following content:
        ```
        llh,pi_1,pi_2,pi_3
        -10.5,0,1,1
        -12.3,0,0,2
        ```

        The function will return:
        ```
        [
            ([(0, 1), (1, 2), (1, 3)], -10.5),
            ([(0, 1), (0, 2), (2, 3)], -12.3)
        ]
        ```
    """
    df = pd.read_csv(filename, sep=sep)
    
    if "llh" not in df.columns:
        sys.stderr.write("Error: '%s' does not contain 'llh' column" % filename)
        sys.exit(1)
        
    tree_columns = [col for col in df.columns if col.startswith('pi_')]
    freq_columns = [tuple(map(int, col.split('_')[1:])) for col in df.columns if col.startswith('f_')]
    trees = []
    for idx, row in df.iterrows():
        F = np.zeros((m,n))
        edges = []
        for col in tree_columns:
            target = int(col.split("_")[1])
            source = int(row[col])
            edges.append((source, target))
        for col in freq_columns:
            F[col] = row["f_%d_%d" % col]
        llh = float(row["llh"])
        trees.append((edges, (F, llh)))
    
    return trees

def parse_input(filename, sep, ignore_clusters):
    """
    Parses a CSV file containing mutation data and computes matrices for variant counts,
    reference counts, and observed variant frequencies.

    The CSV file is expected to contain the following columns:
    - `mutation_index` or `cluster_index`: Identifies the mutation or cluster.
    - `sample_index`: Identifies the sample.
    - `var`: The number of reads supporting the variant.
    - `depth`: The total read depth at the locus.

    If `cluster_index` is not present, it is assumed to be the same as `mutation_index`.

    Args:
        filename (str): Path to the CSV file containing the mutation data.
        sep (str): Delimiter used in the CSV file (e.g., ',' for comma-separated values).
        ignore_clusters (bool): Indicator to ignore provided clustering

    Returns:
        tuple: A tuple containing four numpy arrays:
            - V (numpy.ndarray): A matrix of variant counts, where `V[i][p]` represents the
              number of variant reads for sample `i` and mutation/cluster `p`.
            - R (numpy.ndarray): A matrix of reference counts, where `R[i][p]` represents the
              number of reference reads for sample `i` and mutation/cluster `p`.
            - hat_F (numpy.ndarray): A matrix of observed variant frequencies, where `hat_F[i][p]`
              represents the frequency of the variant for sample `i` and mutation/cluster `p`.

    Raises:
        SystemExit: If the input file is missing any of the required columns (`mutation_index`,
                   `cluster_index`, `sample_index`, `var`, or `depth`), the function will
                   terminate the program with an error message.

    Example:
        Given a CSV file with the following content:
        ```
        sample_index,mutation_index,var,depth
        0,0,10,100
        0,1,5,50
        1,0,20,200
        1,1,10,100
        ```

        The function will return:
        ```
        V = [[10, 5],
             [20, 10]]
        R = [[90, 45],
             [180, 90]]
        hat_F = [[0.1, 0.1],
                 [0.1, 0.1]]
        ```
    """
    if filename == "-":
        df = pd.read_csv(sys.stdin, sep=sep)
    else:
        df = pd.read_csv(filename, sep=sep)
    
    if "mutation_index" not in df.columns and "cluster_index" not in df.columns:
        sys.stderr.write("Error: '%s' does not contain 'mutation_index' or 'cluster_index' column" % filename)
        sys.exit(1)
    
    if "sample_index" not in df.columns:
        sys.stderr.write("Error: '%s' does not contain 'sample_index' column" % filename)
        sys.exit(1)
        
    if "sample_index" not in df.columns:
        sys.stderr.write("Error: '%s' does not contain 'sample_index' column" % filename)
        sys.exit(1)
    
    if "var" not in df.columns:
        sys.stderr.write("Error: '%s' does not contain 'var' column" % filename)
        sys.exit(1)
    
    if "depth" not in df.columns:
        sys.stderr.write("Error: '%s' does not contain 'depth' column" % filename)
        sys.exit(1)
        
    df["ff"] = df.apply(lambda x:float(x["var"])/x["depth"],axis = 1)
    
    
    if "cluster_index" not in df.columns or ignore_clusters:
        n = len(set(df["mutation_index"]))
        m = len(set(df["sample_index"]))
        
        D = np.zeros((m,n),dtype=int)
        V = np.zeros((m,n),dtype=int)
        hat_F = np.zeros((m,n))
        for idx, row in df.iterrows():
            p = int(row["sample_index"])
            i = int(row["mutation_index"])
            D[p,i] = row["depth"]
            V[p,i] = row["var"]
            
        hat_F = V / D
        R = D - V
        
        assert D.shape == V.shape == R.shape == hat_F.shape
        
        return V, R, hat_F
    else:
        n = len(set(df["cluster_index"]))
        m = len(set(df["sample_index"]))
        
        D = np.zeros((m,n),dtype=int)
        V = np.zeros((m,n),dtype=int)
        hat_F = np.zeros((m,n))
        df["ff"] = df.apply(lambda x:float(x["var"])/x["depth"],axis = 1)
        for p in range(m):
            for i in range(n):
                D[p,i] = df[(df["sample_index"]==p) & (df["cluster_index"]==i)]["depth"].median()
                hat_F[p,i] = df[(df["sample_index"]==p) & (df["cluster_index"]==i)]["ff"].mean()
                V[p,i] = int(ceil(D[p,i]*hat_F[p,i]))

        R = D-V
        
        assert D.shape == V.shape == R.shape == hat_F.shape
    
        return V, R, hat_F

def process_output(BBTs, V, R):
    """
    Processes a list of Backbone Trees (BBTs) and their associated log-likelihood values
    into a structured DataFrame. The DataFrame includes tree indices, log-likelihood values,
    parent indices (`pi`) for each mutation, and the inferred frequency matrix (`F`) for
    each sample and mutation.

    Args:
        BBTs (list): A list of tuples, where each tuple contains:
            - A tree object (with `edges` attribute representing the tree structure).
            - The log-likelihood (`llh`) value associated with the tree.

    Returns:
        pandas.DataFrame: A DataFrame containing the following columns:
            - `tree`: The index of the tree.
            - `llh`: The log-likelihood value of the tree.
            - `pi_X`: The parent index for mutation `X` (where `X` is the mutation index).
            - `f_p_X`: The inferred frequency of mutation `X` in sample `p`.

    Example:
        Given an input `BBTs` list:
        ```
        [
            (tree1, -10.5),  # tree1 has edges [(0, 1), (1, 2)]
            (tree2, -12.3)   # tree2 has edges [(0, 1), (0, 2)]
        ]
        ```

        The function will return a DataFrame like:
        ```
           tree   llh  pi_1  pi_2  f_0_1  f_0_2  f_1_1  f_1_2
        0     0 -10.5     0     1    0.1    0.2    0.3    0.4
        1     1 -12.3     0     0    0.2    0.3    0.4    0.5
        ```
    """
    # Extract mutations and number of samples
    mutations = BBTs[0][0].vertices[1:]
    nr_samples = V.shape[0]

    # Create a list to store rows of data
    rows = []

    # Iterate over BBTs and collect data
    for i, (t, (F, llh)) in enumerate(BBTs):
        # Create a dictionary for the current row
        row = {"tree": i, "llh": llh}
        
        # Add pi values (parent indices) for each mutation
        for e in sorted(t.edges, key=lambda x: int(x[1])):
            row[f"pi_{e[1]}"] = e[0]
        
        for p in range(nr_samples):
            for mut_idx, mut in enumerate(mutations):
                row[f"f_{p}_{mut}"] = F[p][mut_idx]

        rows.append(row)

    # Create a DataFrame from the list of rows
    df = pd.DataFrame(rows)

    # Reorder columns to match the original format
    columns = ["tree", "llh"] + [f"pi_{mut}" for mut in mutations] + [f"f_{p}_{mut}" for p in range(nr_samples) for mut in mutations]
    df = df[columns]
    
    return df

def main():
    def handle_args():
        parser = argparse.ArgumentParser(description="Sapling is an algorithm for summarizing and inferring tumor phylogenies from bulk DNA sequencing data")
        
        # Define arguments
        parser.add_argument("-f", type=str, required=True, help="Input filename with mutation read counts (use '-' for STDIN)")
        parser.add_argument("--init_trees", type=str, help="Input filename with initial backbone trees to expand")
        parser.add_argument("-o", type=str, help="Output filename to store trees and frequencies (default: STDOUT)")
        parser.add_argument("--sep", type=str, default="\t", help="Input/output column separator (default: \\t)")
        parser.add_argument("-a", "--rho", type=float, default=0.9, help="Rho parameter, minimum deviation allowed from max likelihood (default: %(default)s, ignored when beam_width specified)")
        parser.add_argument("-t", "--tau", type=int, default=5, help="Tau parameter, maximum number of backbone trees (default: %(default)s)")
        parser.add_argument("-l", "--ell", type=int, default=-1, help="Ell parameter, minimum number of mutations (default: %(default)s, unlimited)")
        parser.add_argument("--big_expand", action="store_true", help="Use big expand (new mutations are anywhere, not just as leaves or splitting a single edge)")
        parser.add_argument("-b", "--beam_width", type=int, default=-1, help="Maximum beam width (default: %(default)s, limited only by --rho)")
        parser.add_argument("-L", type=str, default="fastppm", help="Regression method (default: %(default)s)", choices=bbt_solver.choices)
        parser.add_argument("--alt_roots", action="store_true", help="Explore alternative root nodes")
        parser.add_argument("-m", "--poly_clonal_root", action="store_true", help="Allow poly clonal root node")
        parser.add_argument("--use_clusters", action="store_true", help="Use provided clustering (taking median read depth and using average frequency for variant counts)")
        parser.add_argument("--version", action="version", version=f"Sapling {__version__}", help="Show program's version number and exit")

        # TODO: add multi-threading
        
        # Parse arguments
        return parser.parse_args()

    # Handle arguments
    parameter = handle_args()

    # Use the parameters
    bbt_solver.LLH_method(parameter.L)
    use_big_expand = parameter.big_expand
    
    if parameter.ell > 0:
        #given ell overrides tau
        parameter.tau = -1

    EPS = 1e-3
    llh_EPS = 1e-4

    # Parse input
    V, R, hat_F = parse_input(parameter.f, parameter.sep,not(parameter.use_clusters))
    m, n = V.shape

    BBT_solver = bbt_solver.BackboneTreeSolver(V, R, log(parameter.rho),
                                               EPS=EPS, llh_EPS=llh_EPS, poly_clonal_root=parameter.poly_clonal_root, 
                                               use_big_expand=use_big_expand, beam_width=parameter.beam_width,
                                               alt_roots=parameter.alt_roots)

    if parameter.init_trees is None:
        BBT_solver.init()
    else:
        init_BBTs = parse_input_trees(parameter.init_trees, m, n, parameter.sep)
        BBT_solver.init(init_BBTs)
    
    BBTs = BBT_solver.main(parameter.ell, parameter.tau)
    df = process_output(BBTs, V, R)
    
    # Write the DataFrame to a CSV file
    if parameter.o is None:
        df.to_csv(sys.stdout, index=False, sep=parameter.sep)
    else:
        df.to_csv(parameter.o, index=False, sep=parameter.sep)
    
if __name__ == "__main__":
    main()
