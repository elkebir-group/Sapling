import gurobipy as gp
from gurobipy import GRB
import numpy as np
from math import log
def log_eps(x,eps=1e-6,s_n=3):
    if x < eps:
        return_val = log(eps)
        iter_val = (eps-x)/eps
        for i in range(1,s_n):
            return_val -= (iter_val**i)/i
        return return_val
    return log(x)

def mll(V,R,T,env):
    M = np.full_like(V,0.25,dtype=float)
    Range = 0.25
    n_intervals = 10

    #build the model for once
    m,n = V.shape
    L = M - Range
    U = M + Range
    break_points = [[[(1-break_point/n_intervals)*l+break_point/n_intervals*u 
                      for break_point in range(n_intervals+1) ]
                     for l,u in zip(ll,uu)]
                    for ll,uu in zip(L,U)]

    model = gp.Model(env=env)

    f_vars = [[model.addVar(vtype=GRB.CONTINUOUS,lb=L[p][i],ub=U[p][i]) 
               for i in range(n)] 
              for p in range(m)]
    f_i_vars = [[[model.addVar(vtype=GRB.CONTINUOUS,lb=0,ub=1) 
                  for break_point in range(n_intervals+1)] 
                 for i in range(n)]
                for p in range(m)]

    lin_expr = gp.LinExpr()
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

    lin_expr2 = gp.LinExpr()
    update_constrs = [[None for i in range(n)] for p in range(m)]
    for p in range(m):
        for i in range(n):
            for k in range(n_intervals+1):
                lin_expr += f_i_vars[p][i][k];
                lin_expr2 += f_i_vars[p][i][k] * break_points[p][i][k]
            model.addConstr(lin_expr == 1)
            update_constrs[p][i] = model.addConstr(lin_expr2 == f_vars[p][i])
            lin_expr.clear()
            lin_expr2.clear()

    for p in range(m):
        for i in range(n):
            for k in range(n_intervals+1):
                lin_expr += f_i_vars[p][i][k] * (V[p][i]*log_eps(break_points[p][i][k])+
                                          R[p][i]*log_eps(1-break_points[p][i][k]))
    model.setObjective(lin_expr, GRB.MAXIMIZE)
    lin_expr.clear()

    model.optimize()

    #update model #10 iterations are enough
    for __ in range(10):
        Range*=6/n_intervals
        for p in range(m):
            for i in range(n):
                M[p][i] = f_vars[p][i].X
        L = M - Range
        U = M + Range
        for p in range(m):
            for i in range(n):
                f_vars[p][i].LB=L[p][i]
                f_vars[p][i].UB=U[p][i]


        break_points = [[[(1-break_point/n_intervals)*l+break_point/n_intervals*u 
                          for break_point in range(n_intervals+1) ]
                         for l,u in zip(ll,uu)]
                        for ll,uu in zip(L,U)]

        for p in range(m):
            for i in range(n):
                for k in range(n_intervals+1):
                    model.chgCoeff(update_constrs[p][i], f_i_vars[p][i][k], newvalue=break_points[p][i][k])
        for p in range(m):
            for i in range(n):
                for k in range(n_intervals+1):
                    lin_expr += f_i_vars[p][i][k] * (V[p][i]*log_eps(break_points[p][i][k])+
                                              R[p][i]*log_eps(1-break_points[p][i][k]))
        model.setObjective(lin_expr, GRB.MAXIMIZE)
        lin_expr.clear()
        model.optimize()
    return model.ObjVal
