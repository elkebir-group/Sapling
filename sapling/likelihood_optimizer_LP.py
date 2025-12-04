import gurobipy as gp
from gurobipy import GRB
from likelihood_function_LP import mll

class Likelihood_optimizer:
    """
    A class for optimizing the log-likelihood of a tree structure given variant and reference
    read counts. It uses the Gurobi solver to perform the optimization.

    Attributes:
        m (int): The number of samples.
        n (int): The number of mutations.
        llh_ref (dict): A dictionary to store precomputed log-likelihood values for reference.
        V (numpy.ndarray): A 2D array of size (m x n) representing the variant read counts.
        R (numpy.ndarray): A 2D array of size (m x n) representing the reference read counts.
        env (gurobipy.Env): A Gurobi environment configured for optimization.
    """
    
    def __init__(self,V,R,EPS):
        """
        Initializes the Likelihood_optimizer with variant and reference read counts.

        Args:
            V (numpy.ndarray): A 2D array of size (m x n) representing the variant read counts
                               for each sample and mutation.
            R (numpy.ndarray): A 2D array of size (m x n) representing the reference read counts
                               for each sample and mutation.
            EPS (float): A small value used for numerical stability (not directly used in this class).
        """
        self.m,self.n = V.shape
        self.llh_ref = {}
        self.V = V
        self.R = R
        
        self.env = gp.Env(empty=True)
        self.env.setParam("OutputFlag",0)
        self.env.start()

    def llh(self,T):
        """
        Computes the log-likelihood of a given tree structure. If the log-likelihood for the tree
        has already been computed and cached, it returns the cached value; otherwise, it computes
        the log-likelihood using the Gurobi solver.

        Args:
            T (Tree): A Tree object representing the tree structure.

        Returns:
            float: The log-likelihood of the tree.
        """
        if T.edges in self.llh_ref:
            return self.llh_ref[T.edges]
        llh_res = mll(self.V[:,T.vertices[1:]],self.R[:,T.vertices[1:]],T.children,self.env)
        return llh_res
