#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 17:08:43 2020

@author: burt
vary beta and look at curve form
"""

import numpy as np
import seaborn as sns

from parameters_menten import d
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/ode_models/")

from analysis_module import multi_exp, update_dict
from ode_models import C_model_menten, il2_model_menten, timer_model_menten, core_model_menten
import matplotlib
sns.set(context = "talk", style = "ticks")
# =============================================================================
# define exp conditions
# =============================================================================
model_list = [timer_model_menten, il2_model_menten, C_model_menten]
model_names = ["Timer", "IL2", "K"]
time = np.arange(0,20,0.01)

arr1 = np.logspace(-1,1,51)
arr2 = np.linspace(0,1,50)

pname = "r_p"
arr = arr1 if pname != "r_p" else arr2

d_list = [update_dict(d, val, pname) for val in arr]
sim_names = [pname+str(val) for val in arr]

# =============================================================================
# run experiment
# =============================================================================
exp = multi_exp(time, d_list, model_list, sim_names, model_names, core_model = core_model_menten)

if pname == "r_p":
    norm = matplotlib.colors.Normalize(
        vmin=np.min(arr),
        vmax=np.max(arr))
    
else:
    norm = matplotlib.colors.LogNorm(
        vmin=np.min(arr),
        vmax=np.max(arr))

# choose a colormap
cm = matplotlib.cm.Blues

# create a ScalarMappable and initialize a data structure
sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])

g = sns.relplot(x = "time", y = "cells", kind = "line", data = exp, hue = "sim", 
                col = "model", palette = "Blues", height = 5, legend = False)

   
#g.set(ylim = (0, 30), xlim = (0, time[-1]))
ax = g.axes[0][0]
ax.set_ylabel("cell dens. norm.")
g.set_titles("{col_name}")
cbar = g.fig.colorbar(sm, ax = g.axes, ticks = [arr[0], 1, arr[-1]])
cbar.set_label(pname)

g.set(yscale = "log", ylim = (0.1,None))

g.savefig("../figures/tc_"+pname+".pdf")
