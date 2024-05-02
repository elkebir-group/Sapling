import pandas as pd
import numpy as np
import bbt_solver
from args_handler import *
from math import log

parameter = handle_args()
df = pd.read_csv(parameter["f"],sep="\t")
approx = float(parameter["a"])
tau = int(parameter["t"])
ell = int(parameter["l"])
out = parameter["o"]

EPS = 1e-3
llh_EPS = 1e-4

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

if ell > 0:
    #give ell overrides tau
    tau = -1
BBTs = BBT_solver.main(ell,tau)

with open(out,"w") as fout:
    fout.write("# %d backbone trees, %d mutations\n"%(len(BBTs), len(BBTs[0][0].edges)))
    for i,(t,_) in enumerate(BBTs):
        fout.write("backbone tree %d, llh: %f\n"%(i,_))
        for e in t.edges:
            fout.write("%d %d\n"%e)
