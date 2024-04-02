from cvxopt import solvers, matrix, spdiag, log, spmatrix, mul, div
#import pyBBT
import numpy as np
import random

solvers.options["maxiters"]=200
solvers.options["show_progress"]=False

dnt_cvx_reapply = True

# parameters = {"prune type":"leaf","prune para":None,"log":True}
###############################################50
parameters = {"prune type":"pair","prune para":50,"log":True,"LPN":100,"LPR":10}

# parameters = {"prune type":"mix","prune para":(1.1,0.25),"log":True}
#type:"none" no prune
#type:"pair" unrestricted freq, llh diff
#type:"top" unrestricted freq, 2 children, top-llh-change
#type

def get_edges(seq):
    n=len(seq)+2
    al = list(range(n))
    edges = []
    for i in range(n-2):
        piv = min( set(al).difference(set(seq[i:])) )
        edges.append( (piv,seq[i]) )
        al.remove (piv)
    edges.append( (al[0],al[1]) )
    return edges


def pair_single(v1,r1,v2,r2):
    f1 = v1/(v1+r1)
    f2 = v2/(v2+r2)
    fm = (v1+v2)/(v1+v2+r1+r2)
    if f1 >= f2:
        pass
    else :
        f1,f2 = fm,fm
    return_val = 0.0
    for mul_v,log_v in zip([v1,r1,v2,r2],[f1,1-f1,f2,1-f2]):
        if mul_v > 0:
            return_val += mul_v*log(log_v)
    return return_val

# def triple_llh(vp,rp,v1,r1,v2,r2):
#     def g(vp,rp,v1,r1,v2,r2):
#         fp = vp/(vp+rp)
#         f1 = v1/(v1+r1)
#         f2 = v2/(v2+r2)
#         # if (fp>=f1+f2):
#         #     pass
#         if (v1 == 0):
#             return pair_single(vp,rp,v2,r2)
#         if (v2 == 0):
#             return pair_single(vp,rp,v1,r1)
#         if (vp == 0):
#             return -1e300
#         vs = matrix(np.array([vp,v1,v2]))
#         rs = matrix(np.array([rp,r1,r2]))
#         def F(x=None, z=None):
#             if x is None: return 0,matrix(0.5,(3,1))
#             if min(x) <= 0 or max(x) >= 1.0: return None
#             f = -sum(mul(vs,log(x)))-sum(mul(rs,log(1-x)))
#             df = (-div(vs,x)+div(rs,1-x)).T
#             if z is None: return f,df
#             ddf = mul(vs,(x**-2))+mul(rs,((1-x)**-2))
#             return f,df,spdiag(z[0]*ddf)
#         G = matrix([-1.,1.,1.]).T
#         # print(G.size)
#         h = matrix([0.])
#         a = solvers.cp(F,G=G,h=h)
#         return a["primal objective"]
#     sum_ = 0.
#     for v,r,vv,rr,vvv,rrr in zip(vp,rp,v1,r1,v2,r2):
#         sum_ += g(v,r,vv,rr,vvv,rrr)
#     return -sum_
        

def pair_llh(v1,r1,v2=None,r2=None):
    if v2 is None:
        v2=np.zeros(v1.shape,int)
        r2=np.ones(v1.shape,int)
    # s1,s2=0.0,0.0
    s1=0.0
    for i in range(len(v1)):
        s1 += pair_single(v1[i],r1[i],v2[i],r2[i])
        # s2 += g(v2[i],r2[i],v1[i],r1[i])
    return s1#-s2

class solver:
    def __init__(self,V,R,lapprox,EPS=1e-4,llh_EPS=0.01,neg=False):
        self.V = np.array(V,dtype = int)
        self.R = np.array(R,dtype = int)
        self.m,self.n = self.V.shape
        self.T = self.V+self.R
        self.lapprox = lapprox
        self.EPS = EPS
        self.llh_EPS = llh_EPS
        self.ave_f = self.V/self.T
        self.sof = self.ave_f.sum(axis=0)
        self.rank = list(range(self.n))
        self.rank.sort(key = lambda x:-self.sof[x])
        if neg:self.rank = [-1]+self.rank
        self.BBTs = {}
        self.BBTs_set = {}
        self.GF = np.zeros((self.n,self.n),int)
        # self.GF2 = np.zeros((self.n,self.n,self.n),int)
#         i,j=13,30
#         print(self.ave_f[:,[i,j]])
        
#         s1 = pair_llh(V[:,i],R[:,i],V[:,j],R[:,j])
#         s2 = pair_llh(V[:,j],R[:,j],V[:,i],R[:,i])
#         print(s1,s2)
        if parameters["prune type"] == "none" or parameters["prune type"] == "single": # or parameters["prune type"]=="leaf":
            self.GF.fill(1)
            # self.GF2.fill(1)
        # if parameters["prune type"] == "mix":
        #     diff_factor = parameters["prune para"][0]
        # else :
        diff_factor = parameters["prune para"]
        if parameters["prune type"] == "pair": # or parameters["prune type"] == "mix":
            # if parameters["prune type"] == "pair":
            #     self.GF2.fill(1)
            for i in range(self.n):
                for j in range(i,self.n):
                    if i==j: continue
                    s1 = pair_llh(V[:,i],R[:,i],V[:,j],R[:,j])
                    s2 = pair_llh(V[:,j],R[:,j],V[:,i],R[:,i])
                    if s1-s2 >= diff_factor*lapprox:
                        self.GF[i][j] = 1
                    if s2-s1 >= diff_factor*lapprox:
                        self.GF[j][i] = 1
        # if parameters["prune type"] == "top" or parameters["prune type"] == "mix":
        #     if(parameters["prune type"] == "mix"):
        #         soc = self.GF.sum(axis=0)
        #         mix_flag = True
        #         threshold = max(soc)*parameters["prune para"][1]
        #     single_llhs = [pair_llh(self.V[:,i],self.R[:,i]) for i in range(self.n)]
        #     double_llhs = [[pair_llh(self.V[:,i],self.R[:,i],self.V[:,j],self.R[:,j]) if i!=j else -1e300 for j in range(self.n)] for i in range(self.n)]
        #     opt_llh_delta = [-1e300 for _ in range(self.n)]
        #     for i in range(self.n):
        #         if mix_flag and soc[i] <= threshold: 
        #             self.GF2[:,:,i]=1
        #             continue
        #         # print("recalculating %d"%i)
        #         deltas = [double_llhs[j][i]-single_llhs[j] for j in range(self.n)]
        #         max_delta = max(deltas)
        #         for j in range(self.n):
        #             if deltas[j] >= max_delta + diff_factor*lapprox:
        #                 self.GF[j][i] = 1
        #         for p in range(self.n):
        #             if not self.GF[p][i]: continue
        #             for c in range(self.n):
        #                 if not self.GF[p][c]: continue
        #                 # print(p,c)
        #                 triple_llh_pci = triple_llh(self.V[:,p],self.R[:,p],self.V[:,c],self.R[:,c],self.V[:,i],self.R[:,i])
        #                 triple_delta = triple_llh_pci - double_llhs[p][c]
        #                 if triple_delta >= max_delta + diff_factor*lapprox:
        #                     self.GF2[p][c][i] = 1
                        

    def init_tree(self,init_n):
        def dfs(r,linked_list,res_parent):
            for ch in linked_list[r]:
                if res_parent[r]==ch: continue
                res_parent[ch] = r
                dfs(ch,linked_list,res_parent)

        trees = [None for i in range(init_n**(init_n-1))]
        if self.rank[0]<0:
            rlist = [0]
        else :
            rlist = range(init_n)
        for prufer_id in range(init_n**(init_n-2)):
            seq = [None for i in range(init_n-2)]
            prufer_number = prufer_id
            for i in range(init_n-2):
                seq[i] = prufer_id%init_n
                prufer_id//=init_n
            edges = get_edges(seq)
            link_list = [[] for i in range(init_n)]
            for e in edges:
                link_list[e[0]].append(e[1])
                link_list[e[1]].append(e[0])
            for r in rlist:
                parent = [-1 for i in range(init_n)]
                dfs(r, link_list, parent)
                arcs = [(self.rank[parent[i]],self.rank[i]) for i in range(init_n) if parent[i]>=0]
                arcs.sort()
                trees[prufer_number*init_n+r] = tuple(arcs)
        # print(trees)
        if self.rank[0]<0:
            trees = [_ for _ in trees if _ is not None]
        # print(trees)
        return trees

    def expanded_trees(self,partial_T,expanded_vertex=None):
        if expanded_vertex is None:
            expanded_vertex = self.rank[len(partial_T)+1]
            free_append = False
        else :
            free_append = True
        trees = []
        if parameters["prune type"] == "none" or free_append:
            parent = {}
            children = {}
            for e in partial_T:
                parent[e[1]]=e[0]
                if e[0] in children:
                    children[e[0]].append(e[1])
                else :
                    children[e[0]]=[e[1]]
            setall = set(parent.keys()).union(children.keys())
            for v in setall:
                if v not in children:
                    nt = list(partial_T)
                    nt.append((v,expanded_vertex))
                    nt.sort()
                    trees.append(tuple(nt))
                else :
                    for i in range(1<<len(children[v])):
                        set_split = []
                        for idx in range(len(children[v])):
                            if ((1<<idx)&i):
                                set_split.append(children[v][idx])
                        set_split = set(set_split)
                        nt = [e for e in partial_T if ((e[0]!=v) or (e[1] not in set_split))]\
                                    +[(v,expanded_vertex)]+[(expanded_vertex,nv) for nv in set_split]
                        nt.sort()
                        trees.append(tuple(nt))
                            
        if parameters["prune type"] != "none" and not free_append:
            i = len(partial_T)+1
            for k in self.rank[:i]:
                if (k<0 or self.GF[k][self.rank[i]]):
                    # sibs = [e[1] for e in partial_T if e[0]==k]
                    # flag = True
                    # for sib in sibs:
                    #     if not self.GF2[k][sib][self.rank[i]]:
                    #         flag = False
                    #         break
                    # if (not flag): continue
                    nt = list(partial_T)
                    nt.append((k,self.rank[i]))
                    nt.sort()
                    trees.append(tuple(nt))
            # if parameters["prune type"] == "leaf": return trees
            for e in partial_T:
                if ((e[0]<0 or self.GF[e[0]][self.rank[i]]) and self.GF[self.rank[i]][e[1]]):
                    # sibs = [e_[1] for e_ in partial_T if e_[0]==e[0] and e_[1]!=e[1]]
                    # flag = True
                    # for sib in sibs:
                    #     if not self.GF2[e[0]][sib][self.rank[i]]:
                    #         flag = False
                    #         break
                    # if (not flag): continue
                    nt = list(partial_T)
                    nt.append((e[0],self.rank[i]))
                    nt.append((self.rank[i],e[1]))
                    nt.remove(e)
                    nt.sort()
                    trees.append(tuple(nt))
        # print(trees)
        return trees

def mll(V,R,muts,partial_T):
    tmp_idx = [-1 for _ in V[0]]
    muts = set(np.ravel(np.array(list(partial_T))).tolist())
    if (-1 in muts): muts.remove(-1)
    # muts = list(muts)
    # muts.sort()
    for i,mut in enumerate(muts):
        tmp_idx[mut] = i
    def F(x=None, z=None):
        if x is None: 
            return 0,matrix(0.25,(len(muts)*len(V),1))
        if min(x) <= 0 or max(x) >= 1: return None
        val = 0
        d1 = [None for i in range(len(muts)*len(V))]
        d2 = [None for i in range(len(muts)*len(V))]
        for i in range(len(V)):
            for mut in muts:
                idx_x = i*len(muts)+tmp_idx[mut]
                val+=-V[i][mut]*log(x[idx_x])-R[i][mut]*log(1.-x[idx_x])
                d1[idx_x]=[-V[i][mut]/x[idx_x]+R[i][mut]/(1.-x[idx_x])]
                if z is not None:
                    d2[idx_x] = (V[i][mut]/(x[idx_x]*x[idx_x]) + R[i][mut]/((1.-x[idx_x])*(1.-x[idx_x])))
                    d2[idx_x] *= z[0] 
        if z is None: 
            return val,matrix(d1)
        return val,matrix(d1),spdiag(d2)
    
    children = [[] for _ in range(len(V[0])+1)]
    for e in partial_T:
        children[e[0]].append(e[1])
    # print(partial_T)
    # print(children)
    leq_idx = 0
    leq_idx_list = []
    idx_x_list = []
    coef_list = []
    for mut in muts:
        if len(children[mut])>0:
            for i in range(len(V)):
                # leq_idx =
                idx_x = i*len(muts)+tmp_idx[mut]
                coef_list.append(-1)
                idx_x_list.append(idx_x)
                leq_idx_list.append(leq_idx)
                for ch in children[mut]:
                    idx_ch = i*len(muts)+tmp_idx[ch]
                    coef_list.append(1)
                    idx_x_list.append(idx_ch)
                    leq_idx_list.append(leq_idx)
                leq_idx+=1
    if (len(children[-1]))>0:
        for i in range(len(V)):
            # idx_x = i*len(muts)+tmp_idx[mut]
            # coef_list.append(-1)
            # idx_x_list.append(idx_x)
            # leq_idx_list.append(leq_idx)
            for ch in children[-1]:
                idx_ch = i*len(muts)+tmp_idx[ch]
                coef_list.append(1)
                idx_x_list.append(idx_ch)
                leq_idx_list.append(leq_idx)
            leq_idx+=1
    G=spmatrix(coef_list,leq_idx_list,idx_x_list)
    h=matrix(0.,(leq_idx,1) )
    if (len(children[-1]))>0:
        for i in range(len(V)):
            h[-1-i]=0.5
    res_solver = solvers.cp(F,G=G,h=h)
    # print(res_solver["x"])
    return res_solver

class solver_mix(solver):
    def __init__(self,V,R,lapprox,EPS=1e-4,llh_EPS=0.01,brange=10000,neg=False):
        solver.__init__(self,V,R,lapprox,EPS,llh_EPS,neg)
        self.mll_V=self.V+EPS
        self.mll_R=self.R+EPS
        self.brange = 10000
        # self.LPA_solver = pyBBT.maxllh(self.V.tolist(),(self.V+self.R).tolist(),self.rank,brange,1)
        # self.LPA_pre_solver = pyBBT.maxllh(self.V.tolist(),(self.V+self.R).tolist(),self.rank,parameters["LPN"],1)
        # print(self.GF.tolist())

    def solve_llh_cvxopt(self,partial_T,mut_set = None):
        if mut_set is None:
            solver_a = mll(self.mll_V,self.mll_R,self.rank[:len(partial_T)+1],partial_T)
        else :
            solver_a = mll(self.mll_V,self.mll_R,mut_set,partial_T)
        if solver_a["dual objective"]-solver_a["primal objective"] <= self.llh_EPS:
            return True,solver_a["primal objective"]
        else :
            return False,solver_a["primal objective"]

    def solve_llh_LPA(self,partial_T):
        return self.LPA_solver.llh_snv(list(partial_T))

    def solve_LPA_pre(self,partial_T):
        return self.LPA_pre_solver.llh_snv(list(partial_T))

    # def check_save_set(self,partial_trees):
    #     set_mut = set()
    #     for e in partial_trees[0]:
    #         set_mut.add(e[0])
    #         set_mut.add(e[1])
    #     set_mut = frozenset(set_mut)
    #     if parameters["log"]:
    #         print("%s: exploring %d trees."%(str(set_mut),len(partial_trees)))
    #     solved_to_opt = True
    #     # print(partial_trees,set_mut)
    #     ress = [self.solve_llh_cvxopt(t,set_mut) for t in partial_trees]
    #     llhs_solved_opt = [_[1] for _ in ress if _[0]]
    #     max_llh = min(llhs_solved_opt)
    #     tmp_status = {}
    #     status = True
    #     for t,_ in zip(partial_trees,ress):
    #         if _[1] < max_llh - self.lapprox + self.llh_EPS:
    #             tmp_status[t] = _[1]
    #             if (not _[0]):
    #                 status = False
    #     if status:
    #         self.BBTs_set[set_mut] = tmp_status
    #     else:
    #         if parameters["log"]:
    #             print("%s: cvxopt not solving to opt, fall back to LPA. exploring %d trees."%(str(set_mut),len(tmp_status)))
    #         tmp_status = [(self.solve_llh_LPA(t),t) for t in tmp_status]
    #         max_llh = max(tmp_status)[0]
    #         self.BBTs_set[set_mut] = {t:-_ for _,t in tmp_status if _ > max_llh + self.lapprox - self.llh_EPS}
    #     print("%s: %d trees selected"%(str(set_mut),len(self.BBTs_set[set_mut])))
    #     return self.BBTs_set[set_mut]

    def check_and_save(self,partial_trees):
        # if parameters["log"]:
        #     print("size %d: exploring %d trees."%(len(partial_trees[0])+1,len(partial_trees)))

        # preres = [self.solve_LPA_pre(t) for t in partial_trees]
        # llhs_solved_opt = max(preres)
        # partial_trees = [t for t,l in zip(partial_trees,preres) if l >= llhs_solved_opt + parameters["LPR"]*self.lapprox]
        if parameters["log"]:
            print("size %d: exploring %d trees."%(len(partial_trees[0])+1,len(partial_trees)))
        solved_to_opt = True
        ress = [self.solve_llh_cvxopt(t) for t in partial_trees]
        llhs_solved_opt = [_[1] for _ in ress if _[0]]
        max_llh = min(llhs_solved_opt)
        tmp_status = {}
        status = True
        for t,_ in zip(partial_trees,ress):
            if _[1] < max_llh - self.lapprox + self.llh_EPS:
                tmp_status[t] = _[1]
                if (not _[0]):
                    status = False
        if status or dnt_cvx_reapply:
            self.BBTs[len(partial_trees[0])+1] = tmp_status
        else:
            if parameters["log"]:
                print("%d: cvxopt not solving to opt, fall back to LPA. exploring %d trees."%(len(partial_trees[0])+1,len(tmp_status)))
            tmp_status = [(self.solve_llh_LPA(t),t) for t in tmp_status]
            max_llh = max(tmp_status)[0]
            self.BBTs[len(partial_trees[0])+1] = {t:-_ for _,t in tmp_status if _ > max_llh + self.lapprox - self.llh_EPS}
        print("size %d: %d trees selected"%(len(partial_trees[0])+1,len(self.BBTs[len(partial_trees[0])+1])))

    def greedy_expansion(self,partial_Ts,top=1):
        trees = [ t for partial_T in partial_Ts for t in self.expanded_trees(partial_T) ]
        llhs = [self.solve_llh_cvxopt(t)[1] for t in trees]
        zip_lt = [(llh,t) for t,llh in zip(trees,llhs)]
        zip_lt.sort()
        llhs = [_[0] for _ in zip_lt[:top]]
        ts = [_[1] for _ in zip_lt[:top]]
        return ts,llhs

    # def random_expansion(self,partial_T):
    #     trees = self.expanded_trees(partial_T)
    #     return trees#[idx]#,llhs[idx]

    def tl_expansion(self,partial_Ts,llh_opt=None):
        trees = [ t for partial_T in partial_Ts for t in self.expanded_trees(partial_T) ]
        llhs = [self.solve_llh_cvxopt(t)[1] for t in trees]
        if llh_opt is None:
            max_llh = min(llhs)
        else :
            max_llh = llh_opt
        zip_lt = [(llh,t) for t,llh in zip(trees,llhs) if llh <= max_llh - self.lapprox + self.llh_EPS]
        zip_lt.sort()
        llhs = [_[0] for _ in zip_lt]
        ts = [_[1] for _ in zip_lt]
        return ts,llhs

# def init_tree(init_n, mapping = (lambda x:x)):
#     def dfs(r,linked_list,res_parent):
#         for ch in linked_list[r]:
#             if res_parent[r]==ch: continue
#             res_parent[ch] = r
#             dfs(ch,linked_list,res_parent)

#     trees = [None for i in range(init_n**(init_n-1))]
#     for prufer_id in range(init_n**(init_n-2)):
#         seq = [None for i in range(init_n-2)]
#         prufer_number = prufer_id
#         for i in range(init_n-2):
#             seq[i] = prufer_id%init_n
#             prufer_id//=init_n
#         edges = get_edges(seq)
#         link_list = [[] for i in range(init_n)]
#         for e in edges:
#             link_list[e[0]].append(e[1])
#             link_list[e[1]].append(e[0])
#         for r in range(init_n):
#             parent = [-1 for i in range(init_n)]
#             dfs(r, link_list, parent)
#             arcs = [(mapping(parent[i]),mapping(i)) for i in range(init_n) if parent[i]>=0]
#             trees[prufer_number*init_n+r] = arcs
#     return trees
