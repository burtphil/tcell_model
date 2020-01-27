#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:36:25 2020

@author: burt
this should be run to analyze dynamic range of carrying capacity model

"""
import numpy as np
import seaborn as sns

from parameters import d
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/ode_models/")

from readout_module import get_area, get_peak, get_peaktime
from analysis_module import multi_param, arr_from_d, get_heatmap, plot_heatmap, run_exp
from ode_models import il2_model, timer_model, null_model, C_model
import pandas as pd
sns.set(context = "talk", style = "ticks", rc = {"lines.linewidth": 4})

import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import itertools


# =============================================================================
# run experiment to get normalization conditions
# =============================================================================
time = np.arange(0,200,0.01)

# =============================================================================
# heatmap conditions
# =============================================================================

pnames = ["n_div", "r_diff", "gamma"]
pnames2 = ["r_il2", "r_myc", "r_C"]
arr_spec = [arr_from_d(d, name, res = 20) for name in pnames2]
res = len(arr_spec[0])
arr_com = np.logspace(-1, 1, res)

model_list = [il2_model, timer_model, C_model]

# =============================================================================
# resp size heatmaps
# =============================================================================
titles = ["resp. size IL2", "resp. size Timer", "resp. size K"]
readout_fun = get_area
norm_area = 40

vmin = -3
vmax = 3
color = "bwr"

heatmap_list = []
fig, axes2 = plt.subplots(3,3, figsize = (14,10))

for name_com, axes in zip(pnames, axes2):
    label_com = name_com      
    
    for arr, name_spec, model, title, ax in zip(arr_spec, pnames2, model_list, titles, axes):
        label_spec = name_spec
        
        heatmap_arr = get_heatmap(arr, arr_com, name_spec, name_com, time, d, model, readout_fun, norm = norm_area)
        arr, arr_com, val = heatmap_arr
        cmap = ax.pcolormesh(arr, arr_com, val, cmap = color, vmin = vmin, vmax = vmax,
                             rasterized = True)        
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(label_spec)
        ax.set_ylabel(label_com)
        ax.set_title(title)
              
        cbar = plt.colorbar(cmap, ax = ax)
        cbar.set_label("log2FC")    
        
        heatmap_list.append(heatmap_arr)
        
plt.tight_layout()

fig.savefig("../figures/model_spec_heatmaps.pdf")


val_list = [heatmap_arr[2].T for heatmap_arr in heatmap_list]

y_arr = [[arr_com[np.argmin(np.abs(column))] for column in mesh] for mesh in val_list]

fig, axes = plt.subplots(1, 3, figsize = (14,4))

y1 = [y_arr[0], y_arr[1], y_arr[2]]
y2 = [y_arr[3], y_arr[4], y_arr[5]]
y3 = [y_arr[6], y_arr[7], y_arr[8]]
y_arr = [y1, y2, y3]

for ax, y, ylabel in zip(axes, y_arr, pnames):
    ax.scatter(arr_com[:-1], y[0], label = "IL2")
    ax.scatter(arr_com[:-1], y[1], label = "Timer")
    ax.scatter(arr_com[:-1], y[2], label = "K")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.1,10)
    ax.set_xlabel("model param norm.")
    ax.set_ylabel(ylabel)

axes[2].legend()
plt.tight_layout()

fig.savefig("../figures/normalized_params.pdf")
fig.savefig("../figures/normalized_params.svg")