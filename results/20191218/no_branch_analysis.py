# -*- coding: utf-8 -*-
from tcell_parameters import d_null, d_il7, d_il2, d_timer
from tcell_parameters import d_il7_death, d_il2_death, d_timer_death
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
import matplotlib.pyplot as plt
from test_module import run_exp, vary_param, multi_exp, update_dict, update_dicts
import module_models as models
import numpy as np
import seaborn as sns
import pandas as pd
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})
import matplotlib
# =============================================================================
# define exp conditions
# =============================================================================
cond1 = [d_null, d_il7, d_il2, d_timer]
cond2 = [d_null, d_il7_death, d_il2_death, d_timer_death]

cond_list = [cond1, cond2]
cond_names2 = ["prolif", "death"]
cond_names = ["Null", "K", "IL2", "Timer"]
time = np.arange(0,6,0.01)

model = models.th_cell_diff

arr = np.arange(15,30,2)
cond_list = [update_dicts(cond1, val, "beta_p") for val in arr]
cond_names2 = ["betap"+str(val) for val in arr]
# =============================================================================
# run experiment
# =============================================================================

# run experiment for impact on prolif
#exp1 = run_exp(time, cond1, cond_names, model = model)
#exp1 = exp1[exp1["cell"] == "teff"]
#exp1["cond2"] = "prolif"

#exp2 = run_exp(time, cond2, cond_names, model = model)
#exp2 = exp2[exp2["cell"] == "teff"]
#exp2["cond2"] = "death"

exp = multi_exp(time, cond_list, cond_names, cond_names2)

norm = matplotlib.colors.Normalize(
    vmin=np.min(arr),
    vmax=np.max(arr))

# choose a colormap
cm = matplotlib.cm.Blues

# create a ScalarMappable and initialize a data structure
sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])

g = sns.FacetGrid(exp, hue = "cond2", col = "cond", palette = "Blues", height = 4)
g = g.map(plt.plot, "time", "value")
g.set(ylim = (0,5))
g.set_titles("{col_name}")
#g.get_legend().remove()
g.fig.colorbar(sm)
#plt.tight_layout()

#df2 = generate_readouts(df)
#param_arr = np.arange(10,50,1)
#norm = 10
#pname = "beta_p"
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