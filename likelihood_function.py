from cvxopt import solvers, matrix, spdiag, log, spmatrix, mul, div

solvers.options["show_progress"]=False

def mll(V,R,T,init_values):
    """
    Solves the maximum log-likelihood (MLL) optimization problem for a given tree structure
    using the CVXOPT solver. The function computes the optimal variant frequencies that
    maximize the log-likelihood of the observed data under the given tree constraints.

    Args:
        V (numpy.ndarray): A 2D array of size (m x n) representing the variant read counts
                           for each sample and mutation.
        R (numpy.ndarray): A 2D array of size (m x n) representing the reference read counts
                           for each sample and mutation.
        T (list): A list of tuples representing the edges of the tree. Each tuple contains
                  two integers (parent, child) defining the tree structure.
        init_values (numpy.ndarray or None): Initial values for the optimization variables.
                                             If None, default values of 0.25 are used.

    Returns:
        dict: A dictionary containing the results of the optimization, including:
              - `x`: The optimal solution (variant frequencies).
              - `status`: The status of the optimization (e.g., 'optimal').
              - Other solver-specific information (e.g., iterations, objective value).

    Notes:
        - The function uses the CVXOPT solver to perform constrained optimization.
        - The tree structure is used to define constraints on the variant frequencies,
          ensuring that child nodes do not exceed their parent nodes in frequency.
        - If `init_values` is not provided, the solver starts with default values of 0.25
          for all variables.

    Example:
        Given:
        ```
        V = np.array([[10, 5], [20, 10]])
        R = np.array([[90, 95], [80, 90]])
        T = [(-1, 0), (0, 1)]
        init_values = None
        ```

        The function will return the optimal variant frequencies and solver status.
    """
    m,n = V.shape
    if init_values is None:
        init_values = matrix(0.25,(m*n,1))
    else :
        init_values = matrix(init_values)
        
    def F(x=None, z=None):
        """
        Local function required by CVXOPT for evaluating the objective, gradient, and Hessian.

        Args:
            x (cvxopt.matrix): The current point in the optimization.
            z (cvxopt.matrix): The scaling factor for the Hessian (used in second-order methods).

        Returns:
            tuple: Depending on the input, returns:
                   - If `x` is None: The number of constraints and initial values.
                   - If `z` is None: The objective value and gradient.
                   - Otherwise: The objective value, gradient, and Hessian.
        """
        if x is None: 
            return 0, init_values
        if min(x) <= 0 or max(x) >= 1: return None
        val = 0
        d1 = [None for i in range(n*m)]
        d2 = [None for i in range(n*m)]
        for i in range(m):
            for mut in range(n):
                idx_x = i*n+mut
                val+=-V[i][mut]*log(x[idx_x])-R[i][mut]*log(1.-x[idx_x])
                d1[idx_x]=[-V[i][mut]/x[idx_x]+R[i][mut]/(1.-x[idx_x])]
                if z is not None:
                    d2[idx_x] = (V[i][mut]/(x[idx_x]*x[idx_x]) + R[i][mut]/((1.-x[idx_x])*(1.-x[idx_x])))
                    d2[idx_x] *= z[0] 
        if z is None: 
            return val,matrix(d1)
        return val,matrix(d1),spdiag(d2)

    children = [[] for _ in range(n+1)]
    for e in T:
        children[e[0]].append(e[1])

    leq_idx = 0
    leq_idx_list = []
    idx_x_list = []
    coef_list = []
    for mut in range(n):
        if len(children[mut])>0:
            for i in range(m):
                idx_x = i*n+mut
                coef_list.append(-1)
                idx_x_list.append(idx_x)
                leq_idx_list.append(leq_idx)
                for ch in children[mut]:
                    idx_ch = i*n+ch
                    coef_list.append(1)
                    idx_x_list.append(idx_ch)
                    leq_idx_list.append(leq_idx)
                leq_idx+=1
                
    #GL mutiple mutations: support Pairtree type of output
    for i in range(len(V)):
        for ch in children[-1]:
            idx_ch = i*n+ch
            coef_list.append(1)
            idx_x_list.append(idx_ch)
            leq_idx_list.append(leq_idx)
        leq_idx+=1
    G=spmatrix(coef_list,leq_idx_list,idx_x_list)
    h=matrix(0.,(leq_idx,1) )
    for i in range(len(V)):
        h[-1-i]=0.5
    res_solver = solvers.cp(F,G=G,h=h)
    return res_solver