#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:36:25 2020

@author: burt
"""

from tcell_parameters import d_null, d_il7, d_il2, d_timer
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
import matplotlib.pyplot as plt
from test_module import multi_param
import module_models as models
import numpy as np
import seaborn as sns
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

colors = ["#3498db", "#95a5a6", "#e74c3c", "#34495e"]     
sns.set_palette(colors)     
# =============================================================================
# define exp conditions
# =============================================================================

cond = [d_null, d_il7, d_il2, d_timer]

cond_names = ["Null","Carr. Cap.", "IL2", "Timer"]

time = np.arange(0,6,0.01)

model = models.th_cell_diff

param_names = ["beta", 
               "d_eff", 
               "beta_p"]

param_arrays = [np.arange(5,15,1),
                np.arange(1,2,0.1),
                np.logspace(1,2,50)]

norm_list = [10,
             1.5,
             10]

# =============================================================================
# look at prolif rate
# =============================================================================
df_new = multi_param(param_arrays, param_names, time, cond,
                cond_names, norm_list, model = model)


df1 = df_new.loc[df_new["pname"] == r"$\beta_p$",:]

ylim = (0,5)
title = "Model"
xlabel = r"$\beta_p$"
g = sns.relplot(x = "x", y = "ylog", kind = "line", data = df1, hue = "cond", 
                col = "readout", height = 4,
                facet_kws = {"margin_titles" : True, "sharex" : True, "sharey" : True},
                legend = "full")


for ax in g.axes.flat:
    ax.set_xscale("log")
    
[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
g.set(ylim = ylim)
g.set(ylim = ylim, xlabel = xlabel)
#g.legend.set_title(title)
g.savefig("beta_p_models.svg")

g = sns.relplot(x = "x", y = "ylog", kind = "line", data = df1, hue = "readout", 
                col = "cond", height = 4, col_order = cond_names,
                facet_kws = {"margin_titles" : True, "sharex" : True, "sharey" : True},
                legend = "full")

for ax in g.axes.flat:
    ax.set_xscale("log")
    
[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
g.set(ylim = ylim, xlabel = xlabel)
g.savefig("beta_p_readouts.svg")