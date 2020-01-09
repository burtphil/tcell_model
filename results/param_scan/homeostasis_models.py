#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:36:25 2020

@author: burt
vary important params for homeostasis models
for this, betap should be in unstable range (set this in tcell params)
"""

import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
from tcell_parameters import d_null, d_il7, d_il2, d_timer
import matplotlib.pyplot as plt
from test_module import multi_param, array_from_dict
import module_models as models
import numpy as np
import matplotlib.ticker as ticker
import seaborn as sns
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

#colors = ["#3498db", "#95a5a6", "#e74c3c", "#34495e"]     
#sns.set_palette(colors)     

# =============================================================================
# define exp conditions
# =============================================================================

cond = [d_il7, d_il2, d_timer]

cond_names = ["Carr. Cap.", "IL2", "Timer"]

time = np.arange(0,250,0.01)

model = models.th_cell_diff

param_names = ["beta", 
               "d_eff", 
               "beta_p"]

norm_list = [d_il7[name] for name in param_names]
param_arrays = [array_from_dict(d_il7, pname) for pname in param_names]

df_new = multi_param(param_arrays, param_names, time, cond,
                cond_names, norm_list, model = model)

figname = "pscan_readouts"
loc_major = ticker.LogLocator(base = 10.0, numticks = 100)
loc_minor = ticker.LogLocator(base = 10.0, 
                              subs = np.arange(0.2,1,0.2),
                              numticks = 12)

g = sns.relplot(x = "xnorm", y = "ylog", kind = "line", data = df_new, hue = "cond", 
                col = "readout", row = "pname", height = 4,
                facet_kws = {"margin_titles" : True, "sharex" : True, "sharey" : True},
                legend = "full")

for ax in g.axes.flat:
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(loc_major)
    ax.xaxis.set_minor_locator(loc_minor)    
    
[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')

# =============================================================================
# same plot but now facet for readouts
# =============================================================================
figname = "pscan_models"
g = sns.relplot(x = "xnorm", y = "ylog", kind = "line", data = df_new, hue = "readout", 
                col = "cond", row = "pname", height = 4,
                facet_kws = {"margin_titles" : True, "sharex" : True, "sharey" : True},
                legend = "full")

ylim = (None, None)
g.set(ylim = ylim, ylabel = "log2FC")
for ax in g.axes.flat:
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(loc_major)
    ax.xaxis.set_minor_locator(loc_minor)    
    
[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')


