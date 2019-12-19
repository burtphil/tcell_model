# -*- coding: utf-8 -*-
from tcell_parameters import d_null, d_prec
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
import matplotlib.pyplot as plt
from test_module import run_exp, vary_param
import module_models as models
import numpy as np
import seaborn as sns
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

# =============================================================================
# define exp conditions
# =============================================================================
d_null = dict(d_null)
d_prec = dict(d_prec)

cond = [d_prec, d_prec]

cond_names = ["t1", "t2"]
time = np.arange(0,30,0.01)

model = models.branch_precursor
# =============================================================================
# run experiment
# =============================================================================
exp = run_exp(time, cond, cond_names, model = model)

g = sns.FacetGrid(exp, col = "cell", hue = "cond", height = 5)
g = g.map(plt.plot, "time", "value")
#df2 = generate_readouts(df)
param_arr = np.arange(10,50,1)
norm = 10
pname = "beta_p"
#df2 = vary_param(param_arr, pname, time, cond, cond_names, norm)

# keep only teff cells
#df_eff = df2[df2["cell"] == "teff"]

#fig = sns.FacetGrid(df_eff, col = "cond", hue = "readout", 
#                    height = 5., aspect = 1, palette = "deep")
#fig = fig.map(plt.plot, "x", "ylog")
#fig.set_ylabels("log2FC")
#fig.set_xlabels("beta p")
#fig.set(ylim = (0,10))
#fig.add_legend()
#fig.savefig("readouts.pdf")