#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:36:25 2020

@author: burt
"""

from tcell_parameters import d_il2, d_timer, d_il2_timer
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

cond = [d_il2, d_timer, d_il2_timer]

cond_names = ["IL2","Timer"]

time = np.arange(0,10,0.001)

model = models.th_cell_diff

param_names = ["beta",
               "beta_p"]

res = 50
param_arrays = [np.linspace(0,10, res),
                np.linspace(0.1,100, res)]

norm_list = [param_arrays[0][-1]/2,
             0.2]

# =============================================================================
# look at prolif rate
# =============================================================================
df_new = multi_param(param_arrays, param_names, time, cond,
                cond_names, norm_list, model = model, convert = True)

g = sns.relplot(x = "x", y = "ylog", kind = "line", col = "cond", hue = "readout",
                row = "pname",  data = df_new, 
                facet_kws = {"margin_titles" : True, "sharex":False, "sharey":False})

[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')