from BBT_func import *
import pandas as pd
import sys
import numpy as np
from cvxopt import solvers, matrix, spdiag, log, spmatrix

solvers.options["maxiters"]=200
solvers.options["show_progress"]=False

df = pd.read_csv(sys.argv[1],sep="\t")
approx = float(sys.argv[3])

max_BBT = int(sys.argv[2])


EPS = 1e-3
llh_EPS = .001

# print(set(df["cluster"]))
clusters = list(set(df["cluster"]))
cl_idx = {cl:i for i,cl in enumerate(clusters)}

n = len(clusters)
m = len(set(df["sample_index"]))

# root = -1
# if (len(sys.argv)>4):
#     root = int(sys.argv[4])
#     df = df[:n*m]

if (len(sys.argv)>4):
    m = int(sys.argv[4])
    df = df[:n*m]
# print(df)

depth = np.zeros((m,n),dtype=int)
var = np.zeros((m,n),dtype=int)
ave_f = np.zeros((m,n))
df["ff"] = df.apply(lambda x:float(x["var"])/x["depth"],axis = 1)
sof = np.zeros((n))
for i in range(m):
    for p in range(n):
        # depth[i][p] = df[(df["sample_index"]==i) & (df["cluster"]==clusters[p])]["depth"].sum()
        # var[i][p] = df[(df["sample_index"]==i) & (df["cluster"]==clusters[p])]["var"].sum()
        # ave_f[i][p]=var[i][p]/depth[i][p]
        depth[i][p] = df[(df["sample_index"]==i) & (df["cluster"]==clusters[p])]["depth"].median()
        ave_f[i][p] = df[(df["sample_index"]==i) & (df["cluster"]==clusters[p])]["ff"].mean()
        var[i][p] = int(depth[i][p]*ave_f[i][p])
        sof[p]+=ave_f[i][p]


ref = depth-var

BBT_solver = solver_mix(var,ref,log(approx),EPS=EPS,llh_EPS=llh_EPS,neg=1)#,root=root)
init_size = 2 
init_trees = BBT_solver.init_tree(init_size)

BBT_solver.check_and_save(init_trees)
# print("init tree size %d: %d"%(init_size,len(BBT_solver.BBTs[init_size])))

for i in range(init_size,len(BBT_solver.rank)):
    tmp_trees = []
    for t in BBT_solver.BBTs[i]:
        tmp_trees.extend(BBT_solver.expanded_trees(t))
    BBT_solver.check_and_save(tmp_trees)
    print(BBT_solver.BBTs[i+1])
    if len(BBT_solver.BBTs[i+1]) > max_BBT:# and i >= 0.5*n:
        break

print(BBT_solver.BBTs)

print(clusters)
