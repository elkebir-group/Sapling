import sys
from tree_utils import *

choices = ["cvxopt", "pLP", "fastppm"]

# import likelihood_optimizer_LP import *
def LLH_method(method):
    assert method in choices
    global Likelihood_optimizer 
    if method == "cvxopt":
        import likelihood_optimizer
        Likelihood_optimizer = likelihood_optimizer.Likelihood_optimizer
    elif method == "pLP":
        import likelihood_optimizer_LP
        Likelihood_optimizer = likelihood_optimizer_LP.Likelihood_optimizer
    elif method == "fastppm":
        import likelihood_optimizer_fastppm
        Likelihood_optimizer = likelihood_optimizer_fastppm.Likelihood_optimizer


class BackboneTreeSolver:
    def __init__(self, V, R, log_rho, EPS, llh_EPS, poly_clonal_root, use_big_expand, beam_width):
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

    def init(self, init_BBTs=[]):
        if len(init_BBTs) == 0:
            init_tree = Tree([(-1, self.rank[0])], self.use_big_expand)
            llh = self.solver.llh(init_tree)
            self.BBTs = {1:[(init_tree, llh)]}
        else:
            nr_muts = len(init_BBTs[0][0])
            self.BBTs = {nr_muts: []}
            for t, llh in init_BBTs:
                self.BBTs[nr_muts].append((Tree(t, self.use_big_expand), llh))
        
    def greedy_expand(self, tree):
        partial_trees = tree.expand(self.rank[tree.n], self.poly_clonal_root)
        ress = [self.solver.llh(t) for t in partial_trees]
        max_llh = max(ress)
        next_trees = [(_,t) for t,_ in zip(partial_trees, ress)]
        next_trees.sort()
        return next_trees[-1]

    def iteration(self):
        iter_index = sorted(self.BBTs.keys())[-1]
        partial_trees = []
        for t,_ in self.BBTs[iter_index]:
            partial_trees.extend(t.expand(self.rank[iter_index], self.poly_clonal_root))

        sys.stderr.write("iteration %d: exploring %d trees\n" % (iter_index,len(partial_trees)))
        ress = [self.solver.llh(t) for t in partial_trees]
        max_llh = max(ress)

        self.BBTs[iter_index+1] = [(t, ll) for t,ll in zip(partial_trees,ress) if ll >=max_llh + self.log_rho - self.llh_EPS]
        if self.beam_width != -1:
            # additionally apply beam width
            self.BBTs[iter_index+1] = sorted(self.BBTs[iter_index+1], key=lambda x: -x[1])
            self.BBTs[iter_index+1] = self.BBTs[iter_index+1][:self.beam_width]
            
        sys.stderr.write("iteration %d: %d trees selected\n" % (iter_index,len(self.BBTs[iter_index+1])))

    def main(self, ell, tau=-1):
        assert len(self.BBTs) > 0
        start = sorted(self.BBTs.keys())[0]
        for i in range(start + 1, self.n+1):
            if ell > 0:
                if i-1 >= ell:
                    return self.BBTs[i-1]
            self.iteration()
            if tau > 0:
                if len(self.BBTs[i]) > tau:
                    return self.BBTs[i-1]
        return self.BBTs[self.n]