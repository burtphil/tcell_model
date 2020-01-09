#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:36:25 2020

@author: burt
this should be run to analyze dynamic range of carrying capacity model

"""

import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
from tcell_parameters import d_null, d_il7, d_il2, d_timer
import matplotlib.pyplot as plt
from test_module import multi_param, array_from_dict
import module_models as models
import numpy as np
import seaborn as sns
import matplotlib.ticker as ticker
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

#colors = ["#3498db", "#95a5a6", "#e74c3c", "#34495e"]     
#sns.set_palette(colors)     

# =============================================================================
# define exp conditions
# =============================================================================

cond = [d_il7, d_il2, d_timer]

cond_names = ["Carr. Cap.", "IL2", "Timer"]

time = np.arange(0,30,0.01)

model = models.th_cell_diff

param_names = ["beta_p"]

param_arrays = [np.logspace(1,3,50)]

norm_list = [10]

# =============================================================================
# param scan
# =============================================================================
df_new = multi_param(param_arrays, param_names, time, cond,
                cond_names, norm_list, model = model)


figname = "prolif_only"
df_prolif_only = df_new[(df_new["pname"] == "beta_p") & (df_new["readout"] == "area")]

g = sns.relplot(x = "xnorm", y = "ylog", kind = "line", data = df_prolif_only, 
                hue = "cond", 
                facet_kws = {"margin_titles" : True, "sharex" : False, "sharey" : False},
                legend = "full",
                aspect = 1.2)

xlabel = r"$\beta_p$"
ylabel = "log2FC"
g.set(xlabel = xlabel, ylabel = ylabel, title = "resp. size")
g.set(xscale = "log", xlim = (1,100), ylim = (0, None))
ax = g.axes.flat[0]
loc_major = ticker.LogLocator(base = 10.0, numticks = 100)
loc_minor = ticker.LogLocator(base = 10.0, 
                              subs = np.arange(0.1,1,0.1),
                              numticks = 12)

ax.xaxis.set_major_locator(loc_major)
ax.xaxis.set_minor_locator(loc_minor)

g.savefig("beta_prolif.pdf")
#ax.xaxis.set_minor_formatter(ticker.NullFormatter())