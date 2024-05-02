
from cvxopt import solvers, matrix, spdiag, log, spmatrix, mul, div

solvers.options["show_progress"]=False

def mll(V,R,T,init_values):
    """
    V,R: A 2d np.array of size (m*n)
    T: A tree of size n in form of edges"
    """

    m,n = V.shape
    if init_values is None:
        init_values = matrix(0.25,(m*n,1))
    else :
        init_values = matrix(init_values)
        
    def F(x=None, z=None):
        """
        Local function required by cvxopt
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