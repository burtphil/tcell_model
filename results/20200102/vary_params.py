#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:36:25 2020

@author: burt
"""

from tcell_parameters import d_null, d_il7, d_il2, d_timer, d_prec
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
import matplotlib.pyplot as plt
from test_module import vary_param, get_relative_readouts, multi_param
import module_models as models
import numpy as np
import seaborn as sns
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

# =============================================================================
# define exp conditions
# =============================================================================
#cond = [d_null, d_null]
cond = [d_prec, d_prec]

cond_names = ["a","b"]

time = np.arange(0,6,0.01)

#param_name = "beta"
param_name = "crit_timer"
param_arr = np.arange(1,10,1)
norm = 2

# =============================================================================
# run model
# =============================================================================
#model = models.th_cell_diff
model = models.branch_precursor
df = vary_param(param_arr, param_name, time, cond, cond_names, norm, model = model)

df_rel = get_relative_readouts(df)

param_arrays = [param_arr, param_arr]
param_names = ["crit_timer", "rate_il2"]
norm_list = [norm, norm]

df_new = multi_param(param_arrays, param_names, time, cond,
                cond_names, norm_list, model = model)


g = sns.relplot(x = "x", y = "y", kind = "line", data = df_new, hue = "readout(rel)", 
                row = "pname", col = "cond", height = 4, facet_kws = {"margin_titles" : True},
                legend = "full")

[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')

