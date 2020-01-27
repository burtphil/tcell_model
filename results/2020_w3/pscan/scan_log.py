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

from analysis_module import multi_param
from ode_models import il2_model, timer_model, null_model, C_model
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
sns.set(context = "talk", style = "ticks", rc = {"lines.linewidth": 4})
# =============================================================================
# define exp conditions
# =============================================================================

model_list = [il2_model, timer_model, C_model]
model_names = ["IL2", "Timer", "K"]

time = np.arange(0,200,0.001)
pnames = ["r_diff", "n_div", "gamma"]
res = 51
pnames_arr = [np.logspace(-1,1,res), np.logspace(-1,1,res), np.logspace(-1,1,res)]

# =============================================================================
# param scan
# =============================================================================
df_new = multi_param(pnames, pnames_arr, d, time, model_list, model_names)
 
g = sns.relplot(x = "sim", y = "ylog", kind = "line", data = df_new, 
                hue = "readout", col = "model", row = "pname",
                facet_kws = {"margin_titles" : False, "sharex" : False, "sharey" : True},
                legend = "full",
                aspect = 1.0,
                col_order = ["IL2", "Timer", "K"])
    
loc_major = ticker.LogLocator(base = 10.0, numticks = 100)
loc_minor = ticker.LogLocator(base = 10.0, 
                              subs = np.arange(0.1,1,0.1),
                              numticks = 12)

for ax, pname in zip(g.axes, pnames):
    for a in ax:
        a.set_xlabel(pname)
        a.set_ylabel("log2FC")
        a.xaxis.set_major_locator(loc_major)
        a.xaxis.set_minor_locator(loc_minor)
        a.set_title("")

for ax, title in zip(g.axes[0], model_names):
    ax.set_title(title)
   
g.set(xscale = "log")
g.set(ylim = (-5,5), ylabel = "log2FC")
#g.set_titles(row_template= "{row_var : ""}", col_template="{col_name}")
plt.tight_layout()
g.savefig("../figures/scan_log.svg")
g.savefig("../figures/scan_log.pdf")

# =============================================================================
# same plot with increased axis limits
# =============================================================================

g = sns.relplot(x = "sim", y = "ylog", kind = "line", data = df_new, 
                hue = "readout", col = "model", row = "pname",
                facet_kws = {"margin_titles" : False, "sharex" : False, "sharey" : True},
                legend = "full",
                aspect = 1.0,
                col_order = ["IL2", "Timer", "K"])
    
loc_major = ticker.LogLocator(base = 10.0, numticks = 100)
loc_minor = ticker.LogLocator(base = 10.0, 
                              subs = np.arange(0.1,1,0.1),
                              numticks = 12)

for ax, pname in zip(g.axes, pnames):
    for a in ax:
        a.set_xlabel(pname)
        a.set_ylabel("log2FC")
        a.xaxis.set_major_locator(loc_major)
        a.xaxis.set_minor_locator(loc_minor)
        a.set_title("")

for ax, title in zip(g.axes[0], model_names):
    ax.set_title(title)
   
g.set(xscale = "log")
g.set(ylabel = "log2FC")
#g.set_titles(row_template= "{row_var : ""}", col_template="{col_name}")
plt.tight_layout()
g.savefig("../figures/scan_log_high.svg")
g.savefig("../figures/scan_log_high.pdf")