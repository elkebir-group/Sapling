import BBT_func
from math import log
import pandas as pd
import sys
import numpy as np

df = pd.read_csv(sys.argv[1],sep="\t")

max_BBT = int(sys.argv[2])



EPS = 1e-3
llh_EPS = 0.005

# print(set(df["cluster"]))
clusters = list(set(df["cluster"]))
cl_idx = {cl:i for i,cl in enumerate(clusters)}

n = len(clusters)
m = len(set(df["sample_index"]))

# root = -1
# if (len(sys.argv)>4):
#     root = int(sys.argv[4])
#     df = df[:n*m]

with open(sys.argv[3],"r") as fin:
    lines = [line for line in fin]
BBTs = eval(lines[-2])

if (len(sys.argv)>4):
    m = int(sys.argv[4])
    df = df[:n*m]


depth = np.zeros((m,n),dtype=int)
var = np.zeros((m,n),dtype=int)
ave_f = np.zeros((m,n))
df["ff"] = df.apply(lambda x:float(x["var"])/x["depth"],axis = 1)
sof = np.zeros((n))
for i in range(m):
    for p in range(n):
        # depth[i][p] = df[(df["sample_index"]==i) & (df["cluster"]==clusters[p])]["depth"].sum()
        # var[i][p] = df[(df["sample_index"]==i) & (df["cluster"]==clusters[p])]["var"].sum()
        depth[i][p] = df[(df["sample_index"]==i) & (df["cluster"]==clusters[p])]["depth"].median()
        ave_f[i][p] = df[(df["sample_index"]==i) & (df["cluster"]==clusters[p])]["ff"].mean()
        var[i][p] = int(depth[i][p]*ave_f[i][p])
        sof[p]+=ave_f[i][p]

ref = depth-var

approx = 0.9

BBT_solver = BBT_func.solver_mix(var,ref,log(approx),EPS=EPS,llh_EPS=llh_EPS,neg=0)#,root=root)

sorted_keys = list(BBTs.keys())
sorted_keys.sort()
for i in sorted_keys:
    if len(BBTs[i])> max_BBT:
        i=i-1
        break
return_res = BBTs[i]

for t in BBTs[i]:
    llh = BBTs[i][t]
    while(len(t)<len(BBT_solver.rank)-1):
        (t,),(llh,) = BBT_solver.greedy_expansion([t])
    print(t,":",llh)
