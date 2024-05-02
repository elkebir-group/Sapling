import gurobipy as gp
from gurobipy import GRB
from likelihood_function_LP import mll

class Likelihood_optimizer:
    
    def __init__(self,V,R,EPS):
        """
        V,R: A 2d np.array of size (m*n)
        T: A tree of size n in form of edges"
        """
        self.m,self.n = V.shape
        self.llh_ref = {}
        self.V = V
        self.R = R
        
        self.env = gp.Env(empty=True)
        self.env.setParam("OutputFlag",0)
        self.env.start()

    def llh(self,T):
        if T.edges in self.llh_ref:
            return self.llh_ref[T.edges]
        llh_res = mll(self.V[:,T.vertices[1:]],self.R[:,T.vertices[1:]],T.children,self.env)
        return llh_res
