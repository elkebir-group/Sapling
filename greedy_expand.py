import pandas as pd
import numpy as np
import bbt_solver
from args_handler import *
from math import log

parameter = handle_args()
bbt_solver.LLH_method(parameter["L"])
df = pd.read_csv(parameter["f"],sep="\t")
approx = float(parameter["a"])
out = parameter["o"]
EPS = 1e-3
llh_EPS = 0.005

if "cluster_index" not in df.columns:
    df["cluster_index"] = df["mutation_index"]

n = len(set(df["cluster_index"]))
m = len(set(df["sample_index"]))

depth = np.zeros((m,n),dtype=int)
var = np.zeros((m,n),dtype=int)
ave_f = np.zeros((m,n))
df["ff"] = df.apply(lambda x:float(x["var"])/x["depth"],axis = 1)
for i in range(m):
    for p in range(n):
        depth[i][p] = df[(df["sample_index"]==i) & (df["cluster_index"]==p)]["depth"].median()
        ave_f[i][p] = df[(df["sample_index"]==i) & (df["cluster_index"]==p)]["ff"].mean()
        var[i][p] = int(depth[i][p]*ave_f[i][p])

ref = depth-var

BBT_solver = bbt_solver.BBT_solver(var,ref,log(approx),EPS=EPS,llh_EPS=llh_EPS,neg=parameter["m"])
ref = depth-var

with open(parameter["r"], "r") as fin:
    lines = [line for line in fin]

nbbts = int(lines[0].split()[1])
ell = int(lines[0].split()[4])
n_lines = ell+1
trees = []
llhs = []
for i in range(nbbts):
    tree = []
    llhs.append(float(lines[1+n_lines*i].split()[-1]))
    for j in range(ell):
        v1,v2=lines[2+n_lines*i+j].split()
        tree.append((int(v1),int(v2)))
    trees.append(bbt_solver.Tree(tree))

full_trees = []
full_llhs = []
for i,(t,llh) in enumerate(zip(trees,llhs)):
    while (t.n < n):
        llh,t = BBT_solver.greedy_expand(t)
    full_trees.append(t)
    full_llhs.append(llh)

with open(parameter["o"],"w") as fout:
    fout.write("# %d full trees, %d mutations\n"%(len(full_trees), n))
    for i,(t,_) in enumerate(zip(full_trees,full_llhs)):
        fout.write("backbone tree %d, llh: %f\n"%(i,_))
        for e in t.edges:
            fout.write("%d %d\n"%e)

    