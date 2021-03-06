#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 17:08:43 2020

@author: burt
vary n_div and analyse curve form
"""

from tcell_parameters import d_null, d_il7, d_il2, d_timer
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
import matplotlib.pyplot as plt
from test_module import run_exp, multi_exp, update_dict, update_dicts
import module_models as models
import numpy as np
import seaborn as sns
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})
import matplotlib
# =============================================================================
# define exp conditions
# =============================================================================
cond1 = [d_il7, d_il2, d_timer]
cond_names = ["K", "IL2", "Timer"]
time = np.arange(0,20,0.01)
model = models.th_cell_diff

arr = np.arange(1,2.05,0.05)
cond_list = [update_dicts(cond1, val, "n_div") for val in arr]
cond_names2 = ["ndiv"+str(val) for val in arr]

# =============================================================================
# run experiment
# =============================================================================
exp = multi_exp(time, cond_list, cond_names, cond_names2)

norm = matplotlib.colors.Normalize(
    vmin=np.min(arr),
    vmax=np.max(arr))

# choose a colormap
cm = matplotlib.cm.Blues

# create a ScalarMappable and initialize a data structure
sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])

g = sns.relplot(x = "time", y = "value", kind = "line", data = exp, hue = "cond2", 
                col = "cond", palette = "Blues", height = 5,legend = False)

   
g.set(ylim = (0, 25), xlim = (0,time[-1]))
ax = g.axes[0][0]
ax.set_ylabel("cell dens. norm.")
g.set_titles("{col_name}")
cbar = g.fig.colorbar(sm, ax = g.axes)
cbar.set_label("n_div")

g.savefig("../figures/prolif_tc_ndiv.pdf")
