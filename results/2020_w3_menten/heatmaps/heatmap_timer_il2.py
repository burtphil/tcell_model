#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:36:25 2020

@author: burt
do heatmaps for model specific parameters timer rmyc and ril2 and look at readouts

"""
import numpy as np
import seaborn as sns

from parameters import d
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/ode_models/")

from readout_module import get_area, get_peak, get_peaktime
from analysis_module import multi_param, arr_from_d, get_heatmap, plot_heatmap, run_exp
from ode_models import il2_model, timer_model, null_model, C_model, timer_il2_model
import pandas as pd
sns.set(context = "talk", style = "ticks", rc = {"lines.linewidth": 4})

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt

# =============================================================================
# run experiment to get normalization conditions
# =============================================================================
time = np.arange(0,100,0.01)

# =============================================================================
# heatmap conditions
# =============================================================================
res = 30
arr_myc = arr_from_d(d, "r_myc", res = res)
arr_il2 = arr_from_d(d, "r_il2", res = res)

name1 = "r_myc"
name2 = "r_il2"

model = timer_il2_model

# =============================================================================
# resp size heatmaps
# =============================================================================
readout_fun_list = [get_area, get_peak, get_peaktime]
norm_list = [40., 15, 3.4]

heatmap_list = []
for readout_fun, norm in zip(readout_fun_list, norm_list):
    
    heatmap_arr = get_heatmap(arr_myc, arr_il2, name1, name2, time, d, model, readout_fun, norm = norm)
    heatmap_list.append(heatmap_arr)


color = "bwr"
vlist = [(-3.5,3.5), (-3.5,3.5), (-3.5,3.5)]
titles = ["resp. size", "peak", "tau"]

fig, axes = plt.subplots(1,3, figsize = (14,4))


for heatmap_arr, vlim, title, ax in zip(heatmap_list, vlist, titles, axes):
    
    vmin, vmax = vlim
    arr1, arr2, val = heatmap_arr
    
    cmap = ax.pcolormesh(arr1, arr2, val, cmap = color, vmin = vmin, vmax = vmax,
                         rasterized = True)        
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(name1)
    ax.set_ylabel(name2)
    ax.set_title(title)
          
    cbar = plt.colorbar(cmap, ax = ax)
    cbar.set_label("log2FC")    
        
plt.tight_layout()

fig.savefig("../figures/heatmaps_timer_il2.pdf")
