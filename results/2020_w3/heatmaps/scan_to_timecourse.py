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
time = np.arange(0,20,0.01)
sim_name = "dummy"
model_name = "dummy"

# readouts should be normalized for default param conditions between models
model = C_model
df = run_exp(time, d, model, sim_name, model_name)

area_norm = get_area(df.time, df.cells)
# =============================================================================
# heatmap conditions
# =============================================================================

arr1 = np.logspace(-1, 1, 30)
arr2 = arr1

name1 = "n_div"
name2 = "r_diff"
name3 = "gamma"

pnames = [name1, name2, name3]
model_list = [il2_model, timer_model, C_model]

# =============================================================================
# resp size heatmaps
# =============================================================================
titles = ["resp. size IL2", "resp. size Timer", "resp. size K"]
readout_fun = get_area
norm_area = 40

vmin = -5
vmax = 5
color = "bwr"

fig, axes2 = plt.subplots(3,3, figsize = (14,10))

for names, axes in zip(list(itertools.combinations(pnames,2)), axes2):
    name1, name2 = names
    label1 = name1
    label2 = name2
      

    for model, title, ax in zip(model_list, titles, axes):

        arr1, arr2, val = get_heatmap(arr1, arr2, name1, name2, time, d, model, readout_fun, norm = norm_area)
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

# =============================================================================
# tau heatmaps
# =============================================================================
titles = ["tau IL2", "tau Timer", "tau K"]
readout_fun = get_peaktime
norm_area = 3.4

vmin = -2.
vmax = 2.
color = "bwr"

fig, axes2 = plt.subplots(3,3, figsize = (14,10))

for names, axes in zip(list(itertools.combinations(pnames,2)), axes2):
    name1, name2 = names
    label1 = name1
    label2 = name2
      

    for model, title, ax in zip(model_list, titles, axes):

        arr1, arr2, val = get_heatmap(arr1, arr2, name1, name2, time, d, model, readout_fun, norm = norm_area)
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

fig.savefig("../figures/heatmaps_tau.pdf")

# =============================================================================
# peak heatmaps
# =============================================================================
titles = ["peak IL2", "peak Timer", "peak K"]
readout_fun = get_peak
norm_area = 15.0

vmin = -5
vmax = 5
color = "bwr"

fig, axes2 = plt.subplots(3,3, figsize = (14,10))

for names, axes in zip(list(itertools.combinations(pnames,2)), axes2):
    name1, name2 = names
    label1 = name1
    label2 = name2
      

    for model, title, ax in zip(model_list, titles, axes):

        arr1, arr2, val = get_heatmap(arr1, arr2, name1, name2, time, d, model, readout_fun, norm = norm_area)
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

fig.savefig("../figures/heatmaps_peak.pdf")