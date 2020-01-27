#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 17:08:43 2020

@author: burt
vary beta and look at curve form
"""

import numpy as np
import seaborn as sns

from parameters import d
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/ode_models/")

from analysis_module import run_exp, multi_exp, generate_readouts, update_dict
from ode_models import il2_model, timer_model, null_model, C_model
import pandas as pd
import matplotlib
sns.set(context = "talk", style = "ticks")
# =============================================================================
# define exp conditions
# =============================================================================
model_list = [il2_model, timer_model, C_model]
model_names = ["IL2", "Timer", "K"]
time = np.arange(0,12.5,0.01)

arr = np.logspace(-1,1,51)

d_list = [update_dict(d, val, "r_diff") for val in arr]
sim_names = ["r_diff"+str(val) for val in arr]

# =============================================================================
# run experiment
# =============================================================================
exp = multi_exp(time, d_list, model_list, sim_names, model_names)

norm = matplotlib.colors.LogNorm(
    vmin=np.min(arr),
    vmax=np.max(arr))

# choose a colormap
cm = matplotlib.cm.Blues

# create a ScalarMappable and initialize a data structure
sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])

g = sns.relplot(x = "time", y = "cells", kind = "line", data = exp, hue = "sim", 
                col = "model", palette = "Blues", height = 5,legend = False)

   
#g.set(ylim = (0, 30), xlim = (0, time[-1]))
ax = g.axes[0][0]
ax.set_ylabel("cell dens. norm.")
g.set_titles("{col_name}")
cbar = g.fig.colorbar(sm, ax = g.axes, ticks = [arr[0], 1, arr[-1]])
cbar.set_label("r diff")

g.set(yscale = "log", ylim = (0.1, 100))
g.savefig("../figures/prolif_tc_r_diff_log.pdf")
