#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 12:45:16 2020

@author: burt
this is to anaylze parameters of Null model, for this betap should be in stable range
"""

from tcell_parameters import d_null

import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
#sys.path.append("C:/Users/Philipp/Documents/projects/tcell_model/code")
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

cond = [d_null]

cond_names = ["Null Model"]

time = np.arange(0,20,0.01)

model = models.th_cell_diff

param_names = ["beta", 
               "d_eff", 
               "beta_p"]

norm_list = [d_null[name] for name in param_names]
param_arrays = [array_from_dict(d_null, pname) for pname in param_names]

df_new = multi_param(param_arrays, param_names, time, cond,
                cond_names, norm_list, model = model, adjust_time = False)

# =============================================================================
# same plot but now facet for readouts
# =============================================================================
figname = "pscan_models"
g = sns.relplot(x = "xnorm", y = "ylog", kind = "line", data = df_new, hue = "readout", 
                col = "pname", height = 4,
                facet_kws = {"margin_titles" : True, "sharex" : True, "sharey" : True},
                legend = "full")

ylim = (None, None)
g.set(ylim = ylim, ylabel = "log2FC")
for ax in g.axes.flat:
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(ticker.LogLocator(base = 10.0, numticks = 100))
    
[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
g.savefig("../figures/pscan_null.pdf")
