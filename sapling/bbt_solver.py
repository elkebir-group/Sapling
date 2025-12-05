import sys
from .tree_utils import *
from tqdm import tqdm

choices = ["cvxopt", "pLP", "fastppm"]

def LLH_method(method):
    """
    Configures the likelihood optimization method to be used by the solver.

    Args:
        method (str): The optimization method to use. Must be one of the following:
                     - "cvxopt": Uses the CVXOPT-based likelihood optimizer.
                     - "pLP": Uses the pLP-based likelihood optimizer.
                     - "fastppm": Uses the fastppm-based likelihood optimizer.

    Raises:
        AssertionError: If the specified method is not in the list of supported choices.
    """
    assert method in choices
    global Likelihood_optimizer 
    if method == "cvxopt":
        from . import likelihood_optimizer_cvxopt
        Likelihood_optimizer = likelihood_optimizer_cvxopt.Likelihood_optimizer
    elif method == "pLP":
        from . import likelihood_optimizer_LP
        Likelihood_optimizer = likelihood_optimizer_LP.Likelihood_optimizer
    elif method == "fastppm":
        from . import likelihood_optimizer_fastppm
        Likelihood_optimizer = likelihood_optimizer_fastppm.Likelihood_optimizer


class BackboneTreeSolver:
    """
    A solver for inferring backbone trees from variant and reference read counts.

    The solver uses a beam search algorithm to iteratively expand and evaluate candidate trees
    based on their log-likelihood scores. It supports multiple likelihood optimization backends.

    Attributes:
        V (numpy.ndarray): A 2D array of variant read counts (samples x mutations).
        R (numpy.ndarray): A 2D array of reference read counts (samples x mutations).
        m (int): The number of samples.
        n (int): The number of mutations.
        D (numpy.ndarray): A 2D array of total read counts (V + R).
        log_rho (float): Approximation threshold in log space for pruning candidate trees.
        EPS (float): A small value for numerical stability.
        llh_EPS (float): A small value for log-likelihood comparison tolerance.
        hatF (numpy.ndarray): A 2D array of observed variant frequencies (V / D).
        sof (numpy.ndarray): A 1D array of the sum of frequencies across samples.
        rank (list): A list of mutation indices sorted by descending sum of frequencies.
        solver (Likelihood_optimizer): The likelihood optimizer instance.
        poly_clonal_root (bool): Whether to allow polyclonal roots in the tree.
        filters (set or None): A set of allowed edges for tree expansion.
        use_big_expand (bool): Whether to use the `big_expand` method for tree expansion.
        beam_width (int): The maximum number of candidate trees to retain during beam search.
        BBTs (dict): A dictionary storing candidate trees at each iteration step.

    Methods:
        __init__: Initializes the solver with input data and parameters.
        init: Initializes the candidate trees for the solver.
        iteration: Performs one iteration of the beam search algorithm.
        main: Runs the solver for a specified number of mutations or trees.
    """
    
    def __init__(self, V, R, log_rho, EPS, llh_EPS, poly_clonal_root, use_big_expand, beam_width, alt_roots):
        """
        Initializes the BackboneTreeSolver with input data and parameters.

        Args:
            V (numpy.ndarray): A 2D array of variant read counts (samples x mutations).
            R (numpy.ndarray): A 2D array of reference read counts (samples x mutations).
            log_rho (float): Approximation threshold in log space for pruning candidate trees.
            EPS (float): A small value for numerical stability.
            llh_EPS (float): A small value for log-likelihood comparison tolerance.
            poly_clonal_root (bool): Whether to allow polyclonal roots in the tree.
            use_big_expand (bool): Whether to use the `big_expand` method for tree expansion.
            beam_width (int): The maximum number of candidate trees to retain during beam search.
            alt_roots (bool): Whether to explore alternative roots.
        """
        # variable read counts
        self.V = V
        # reference read counts
        self.R = R
        # (number of samples, number of mutations)
        self.m, self.n = self.V.shape
        # total read counts
        self.D = self.V+self.R
        # approximation threshold (in log space)
        self.log_rho = log_rho
        
        self.EPS = EPS
        
        self.llh_EPS = llh_EPS
        
        self.hatF = self.V / self.D
        
        # sum of frequencies
        self.sof = self.hatF.sum(axis=0)
        
        self.rank = list(range(self.n))
        
        self.rank.sort(key = lambda x:-self.sof[x])
        
        self.solver = Likelihood_optimizer(V, R, EPS)
        
        self.poly_clonal_root = poly_clonal_root
        
        self.filters = None
        
        self.use_big_expand = use_big_expand
        
        self.beam_width = beam_width
        
        self.alt_roots = alt_roots

    def init(self, init_BBTs=[]):
        """
        Initializes the candidate trees for the solver.

        Args:
            init_BBTs (list): A list of initial backbone trees. Each tree is represented as a tuple
                            (edges, (F, log-likelihood)). If empty, the solver starts with n trees,
                            each containing a distinct mutation.
        """
        if len(init_BBTs) == 0:
            if self.alt_roots:
                self.BBTs = {1 : []}
                for i in range(self.n):
                    init_tree = Tree([(-1, self.rank[i])], self.use_big_expand)
                    F, llh = self.solver.regress(init_tree)
                    self.BBTs[1].append((init_tree, (F, llh)))
            else:
                init_tree = Tree([(-1, self.rank[0])], self.use_big_expand)
                F, llh = self.solver.regress(init_tree)
                self.BBTs = {1:[(init_tree, (F, llh))]}
        else:
            nr_muts = len(init_BBTs[0][0])
            self.BBTs = {nr_muts: []}
            for t, (F, llh) in init_BBTs:
                self.BBTs[nr_muts].append((Tree(t, self.use_big_expand), (F, llh)))

    def iteration(self):
        """
        Performs one iteration of the beam search algorithm.

        During each iteration, the solver expands all candidate trees, evaluates their log-likelihood
        scores, and retains the top candidates based on the beam width and approximation threshold.
        """
        iter_index = sorted(self.BBTs.keys())[-1]
        partial_trees = []
        
        for t,_ in self.BBTs[iter_index]:
            new_mut = None
            blacklist_mut = set(t.vertices)
            for mut in self.rank:
                if mut not in blacklist_mut:
                    new_mut = mut
                    break
            partial_trees.extend(t.expand(new_mut, self.poly_clonal_root))

        # sys.stderr.write("iteration %d: exploring %d trees\n" % (iter_index,len(partial_trees)))
        # ress = [self.solver.regress(t) for t in partial_trees]
        ress = self.solver.regress_multiple_trees(partial_trees)
        max_llh = max([llh for (_, llh) in ress])

        if self.beam_width != -1:
            # apply beam width
            self.BBTs[iter_index+1] = [(t, (F,ll)) for (t,(F,ll)) in zip(partial_trees,ress)]
            self.BBTs[iter_index+1] = sorted(self.BBTs[iter_index+1], key=lambda x: -x[1][1])
            self.BBTs[iter_index+1] = self.BBTs[iter_index+1][:self.beam_width]
        else:
            self.BBTs[iter_index+1] = [(t, (F,ll)) for (t,(F,ll)) in zip(partial_trees,ress) if ll >=max_llh + self.log_rho - self.llh_EPS]

    def main(self, ell, tau=-1):
        """
        Runs the solver for a specified number of mutations or trees.

        Args:
            ell (int): The maximum number of mutations to include in the backbone trees.
                    If <= 0, the solver runs until all mutations are included.
            tau (int): The maximum number of trees to return. If <= 0, the solver returns all trees.

        Returns:
            list: A list of candidate trees, each represented as a tuple (edges, log-likelihood).
        """
        assert len(self.BBTs) > 0
        start = sorted(self.BBTs.keys())[0]
        for i in tqdm(range(start + 1, self.n + 1)):
            if ell > 0:
                if i-1 >= ell:
                    return self.BBTs[i-1]
            self.iteration()
            if tau > 0:
                if len(self.BBTs[i]) > tau:
                    return self.BBTs[i-1]
            # clear previous results
            if i > start + 1:
                self.BBTs[i-1].clear()
        return self.BBTs[self.n]
