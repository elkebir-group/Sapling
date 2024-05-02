from tree_utils import *
from likelihood_optimizer_LP import *

class BBT_solver:
    def __init__(self,V,R,lapprox,EPS,llh_EPS,neg):
        self.V = V
        self.R = R
        self.m,self.n = self.V.shape
        self.T = self.V+self.R
        self.lapprox = lapprox
        self.EPS = EPS
        self.llh_EPS = llh_EPS
        self.ave_f = self.V/self.T
        self.sof = self.ave_f.sum(axis=0)
        self.rank = list(range(self.n))
        self.rank.sort(key = lambda x:-self.sof[x])
        self.solver = Likelihood_optimizer(V,R,EPS)
        self.neg = neg
        self.filters = None

    def init(self):
        init_tree = Tree([(-1,self.rank[0])])
        llh = self.solver.llh(init_tree)
        self.BBTs = {1:[(init_tree,llh)]}

    def greedy_expand(self,tree):
        partial_trees = tree.expand(self.rank[tree.n],self.neg)
        ress = [self.solver.llh(t) for t in partial_trees]
        max_llh = max(ress)
        next_trees = [(_,t) for t,_ in zip(partial_trees,ress)]
        next_trees.sort()
        return next_trees[-1]

    def iteration(self):
        iter_index = len(self.BBTs)
        partial_trees = []
        for t,_ in self.BBTs[iter_index]:
            partial_trees.extend(t.expand(self.rank[iter_index],self.neg))
        print("iteration %d: exploring %d trees"%(iter_index,len(partial_trees)))
        ress = [self.solver.llh(t) for t in partial_trees]
        max_llh = max(ress)
        self.BBTs[iter_index+1] = [(t,_) for t,_ in zip(partial_trees,ress) if _ >= max_llh + self.lapprox - self.llh_EPS]
        print("iteration %d: %d trees selected"%(iter_index,len(self.BBTs[iter_index+1])))

    def main(self,ell,tau=-1):
        self.init()
        for i in range(2,self.n+1):
            if ell > 0:
                if i-1 >= ell:
                    return self.BBTs[i-1]
            self.iteration()
            if tau > 0:
                if len(self.BBTs[i]) > tau:
                    return self.BBTs[i-1]
        return self.BBTs[self.n]