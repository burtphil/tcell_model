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
from test_module import vary_param, run_exp
import module_models as models
import numpy as np
import seaborn as sns
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

colors = ["#3498db", "#95a5a6", "#e74c3c", "#34495e"]     
sns.set_palette(colors)     
# =============================================================================
# define exp conditions
# =============================================================================

cond = [d_il7, d_il2, d_timer]

cond_names = ["Carr. Cap.", "IL2", "Timer"]

time = np.arange(0,6,0.01)

model = models.th_cell_diff

param_name = "beta"

param_arr = np.arange(5,20,0.5)

norm = 10

df = run_exp(time, cond, cond_names)

# =============================================================================
# look at prolif rate
# =============================================================================
df_new = vary_param(param_arr, param_name, time, cond, cond_names, norm, model = model)

xlabel = r"$\beta$"
g = sns.relplot(x = "x", y = "ylog", kind = "line", col = "cond", hue = "readout", data = df_new)
g.set(xlabel = xlabel)
g.savefig("vary_beta.svg")