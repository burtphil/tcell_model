#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:36:25 2020

@author: burt
this should be run to analyze dynamic range of carrying capacity model

"""
import numpy as np
import seaborn as sns

from parameters_menten import d
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/ode_models/")

from ode_models import C_model_menten, il2_model_menten, timer_model_menten, core_model_menten
from readout_module import get_area, get_peaktime
from analysis_module import get_heatmap
import matplotlib.pyplot as plt
import itertools
sns.set(context = "talk", style = "ticks", rc = {"lines.linewidth": 4})


# =============================================================================
# run experiment to get normalization conditions
# =============================================================================
time = np.arange(0, 200, 0.01)

# =============================================================================
# heatmap conditions
# =============================================================================

res = 30
arr1 = np.logspace(-1, 1, res)
arr2 = arr1

name1 = "n_div"
name2 = "r_diff"
label1 = name1
label2 = name2

pnames = [name1, name2]
model_list = [il2_model_menten, timer_model_menten]

# =============================================================================
# resp size heatmaps
# =============================================================================
titles = ["resp. size IL2", "resp. size Timer"]
readout_fun = get_area
norm_area = 18.5

vmin = -5
vmax = 5
color = "bwr"

figsize = (10,4)

fig, axes = plt.subplots(1,2, figsize = figsize)
     
for model, title, ax in zip(model_list, titles, axes):

    arr1, arr2, val = get_heatmap(arr1, arr2, name1, name2, time, d, model, readout_fun, norm = norm_area,
                                  core_model = core_model_menten)
    
    cmap = ax.pcolormesh(arr1, arr2, val, cmap = color, vmin = vmin, vmax = vmax,
                         rasterized = True)        
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(label1)
    ax.set_ylabel(label2)
    ax.set_title(title)
          
    cbar = plt.colorbar(cmap, ax = ax)
    cbar.set_label("log2FC")    
        
plt.tight_layout()

fig.savefig("../figures/heatmaps_area.pdf")

