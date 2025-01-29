import gurobipy as gp
from gurobipy import GRB
import numpy as np
from math import log

def log_eps(x,eps=1e-6,s_n=3):
    """
    Computes the natural logarithm of `x` with a smooth approximation for values below `eps`.

    Args:
        x (float): The input value.
        eps (float): The threshold below which the approximation is used. Default is 1e-6.
        s_n (int): The number of terms in the Taylor series approximation. Default is 3.

    Returns:
        float: The natural logarithm of `x`, or its approximation if `x < eps`.
    """
    if x < eps:
        return_val = log(eps)
        iter_val = (eps-x)/eps
        for i in range(1,s_n):
            return_val -= (iter_val**i)/i
        return return_val
    return log(x)

def mll(V,R,T,env):
    """
    Solves the maximum log-likelihood (MLL) optimization problem for a given tree structure
    using the Gurobi solver. The function computes the optimal variant frequencies that
    maximize the log-likelihood of the observed data under the given tree constraints.

    Args:
        V (numpy.ndarray): A 2D array of size (m_ x n) representing the variant read counts
                           for each sample and mutation.
        R (numpy.ndarray): A 2D array of size (m_ x n) representing the reference read counts
                           for each sample and mutation.
        T (list): A list of lists representing the tree structure. Each sublist contains
                  the children of a node, where `T[-1]` represents the root node's children.
        env (gurobipy.Env): A Gurobi environment configured for optimization.

    Returns:
        float: The total log-likelihood value for the given data and tree structure.

    Notes:
        - The function uses piecewise linear approximation to handle the logarithmic terms
          in the objective function.
        - The optimization is performed iteratively, refining the search range for the
          variant frequencies in each iteration.
        - The Gurobi solver is used to handle the constrained optimization problem.

    Example:
        Given:
        ```
        V = np.array([[10, 5], [20, 10]])
        R = np.array([[90, 95], [80, 90]])
        T = [[1], []]  # Tree structure where node 0 is the parent of node 1
        env = gp.Env(empty=True)
        env.setParam("OutputFlag", 0)
        env.start()
        ```

        The function will return the total log-likelihood value for the given data and tree.
    """
    tot_answer = 0.
    m_,n = V.shape
    n_intervals = 6
    m = 1
    M = np.zeros((1,n),dtype=float)
    
    model = gp.Model(env=env)

    f_vars = [[model.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1) 
               for i in range(n)] 
              for p in range(m)]
    
    f_i_vars = [[[model.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1) 
                  for break_point in range(n_intervals+1)] 
                 for i in range(n)]
                for p in range(m)]
    
    lin_expr = gp.LinExpr()
    lin_expr2 = gp.LinExpr()
    
    for p in range(m):
        for i in T[-1]:
            lin_expr += f_vars[p][i]
        model.addConstr(lin_expr<=0.5)
        lin_expr.clear()
    for p in range(m):
        for i in range(n):
            for j in T[i]:
                lin_expr += f_vars[p][j]
            model.addConstr(lin_expr<=f_vars[p][i])
            lin_expr.clear()
    
    update_constrs = [[None for i in range(n)] for p in range(m)]
    for p in range(m):
        for i in range(n):
            for k in range(n_intervals+1):
                lin_expr += f_i_vars[p][i][k];
                lin_expr2 += f_i_vars[p][i][k] # * break_points[p][i][k] placeholder only
            model.addConstr(lin_expr == 1)
            update_constrs[p][i] = model.addConstr(lin_expr2 == f_vars[p][i])
            lin_expr.clear()
            lin_expr2.clear()

    for _ in range(m_):
        M.fill(0.25)
        Range = 0.25
        for __ in range(15):
            L = M - Range
            U = M + Range
            break_points = [[[(1-break_point/n_intervals)*l+break_point/n_intervals*u 
                              for break_point in range(n_intervals+1) ]
                             for l,u in zip(ll,uu)]
                            for ll,uu in zip(L,U)]

            for p in range(m):
                for i in range(n):
                    f_vars[p][i].LB=L[p][i]
                    f_vars[p][i].UB=U[p][i]

            for p in range(m):
                for i in range(n):
                    for k in range(n_intervals+1):
                        model.chgCoeff(update_constrs[p][i], f_i_vars[p][i][k], newvalue=break_points[p][i][k])

            for p in range(m):
                for i in range(n):
                    for k in range(n_intervals+1):
                        lin_expr += f_i_vars[p][i][k] * (V[_][i]*log_eps(break_points[p][i][k])+
                                                  R[_][i]*log_eps(1-break_points[p][i][k]))
            
            model.setObjective(lin_expr, GRB.MAXIMIZE)
            
            lin_expr.clear()
            model.optimize()
            Range*=4/n_intervals
            for p in range(m):
                for i in range(n):
                    M[p][i] = f_vars[p][i].X
                    
        tot_answer += model.ObjVal

    return tot_answer
        
