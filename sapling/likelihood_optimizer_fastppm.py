import fastppm
import numpy as np

class Likelihood_optimizer:
    """
    A class for optimizing the log-likelihood of a tree structure given variant and reference
    read counts. It uses the `fastppm` library to perform regression and compute the likelihood.

    Attributes:
        m (int): The number of samples.
        n (int): The number of mutations.
        llh_ref (dict): A dictionary to store precomputed log-likelihood values for reference.
        V (numpy.ndarray): A 2D array of size (m x n) representing the variant read counts.
        R (numpy.ndarray): A 2D array of size (m x n) representing the reference read counts.
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
        self.m, self.n = V.shape
        self.llh_ref = {}
        self.V = V
        self.R = R

    def llh(self, T):
        """
        Performs regression to compute the log-likelihood of a given tree structure.

        Args:
            T (Tree): A Tree object representing the tree structure.

        Returns:
            float: The log-likelihood of the tree (to be maximized).
        """
        adj_list = [[] for i in range(len(T.vertices[1:]))]     # create adjacency list for fastppm
        for i, j in T.edges_relabel_from_zero:
            if i == -1: continue
            adj_list[i].append(j)

        # subset read count matrices
        V_prime = self.V[:,T.vertices[1:]] 
        R_prime = self.R[:,T.vertices[1:]]
        T_prime = V_prime + R_prime
        res = fastppm.regress_counts(adj_list, V_prime, T_prime, loss_function="binomial") 
        return -res['objective']
    
    def regress(self, T):
        """
        Performs regression to estimate the frequency matrix for a given tree structure.

        Args:
            T (Tree): A Tree object representing the tree structure.

        Returns:
            numpy.ndarray, float: A frequency matrix of size (m x n), where m is the number of samples
                                  and n is the number of mutations in the tree. Each entry represents
                                  the estimated frequency of a mutation in a sample.
                                  The log-likelihood of the tree (to be maximized).
        """
        adj_list = [[] for i in range(len(T.vertices[1:]))]     # create adjacency list for fastppm
        for i, j in T.edges_relabel_from_zero:
            if i == -1: continue
            adj_list[i].append(j)
            
        # subset read count matrices
        V_prime = self.V[:,T.vertices[1:]] 
        R_prime = self.R[:,T.vertices[1:]]
        T_prime = V_prime + R_prime
        res = fastppm.regress_counts(adj_list, V_prime, T_prime, loss_function="binomial") 
        return np.clip(0, 1, res['frequency_matrix']), -res['objective']
    
    def regress_multiple_trees(self, trees):
        """
        Performs regression to estimate the frequency matrix for given tree structures.

        Args:
            trees (list): A list of Tree objects representing the tree structures.

        Returns:
            numpy.ndarray, float: A frequency matrix of size (m x n), where m is the number of samples
                                  and n is the number of mutations in the tree. Each entry represents
                                  the estimated frequency of a mutation in a sample.
                                  The log-likelihood of the tree (to be maximized).
        """
        adj_lists = []
        for T in trees:
            adj_list = [[] for i in range(len(T.vertices[1:]))]     # create adjacency list for fastppm
            for i, j in T.edges_relabel_from_zero:
                if i == -1: continue
                adj_list[i].append(j)
            adj_lists.append(adj_list)
            
        # subset read count matrices
        V_prime = self.V[:,T.vertices[1:]] 
        R_prime = self.R[:,T.vertices[1:]]
        T_prime = V_prime + R_prime
        
        res = fastppm.regress_counts_multiple_trees(adj_lists, V_prime, T_prime, loss_function="binomial") 
        
        results = []
        for idx in range(len(trees)):
            results.append((np.clip(0, 1, res[idx]['frequency_matrix']), -res[idx]['objective']))
            
        return results