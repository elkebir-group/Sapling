from likelihood_function import mll
from tree_utils import *

# parameters = {"zero_adjustment":1e-3}

class Likelihood_optimizer:
    def __init__(self,V,R,EPS):
        """
        V,R: A 2d np.array of size (m*n)
        T: A tree of size n in form of edges"
        """
        self.m,self.n = V.shape
        self.V = V+EPS
        self.R = R+EPS
        self.llh_ref = {}

    def llh(self,T):
        if T.edges in self.llh_ref:
            return self.llh_ref[T.edges]
        init_values = [0.1*p/self.n for p in T.subtree_size]*self.m
        llh_res = mll(self.V[:,T.vertices[1:]],self.R[:,T.vertices[1:]],T.edges_relabel_from_zero,init_values)
        self.llh_ref[T.edges] = -llh_res["primal objective"]
        return -llh_res["primal objective"]
        
        

