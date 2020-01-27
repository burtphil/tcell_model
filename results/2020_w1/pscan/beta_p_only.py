#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:36:25 2020

@author: burt
this should be run to analyze dynamic range of carrying capacity model

"""
from tcell_parameters import d_null, d_il7, d_il2, d_timer

import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
#sys.path.append("C:/Users/Philipp/Documents/projects/tcell_model/code")

import matplotlib.pyplot as plt
from test_module import multi_param, array_from_dict
import module_models as models
import numpy as np
import seaborn as sns
import matplotlib.ticker as ticker
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

# =============================================================================
# define exp conditions
# =============================================================================

cond = [d_il7, d_il2, d_timer]
cond_names = ["K", "IL2", "Timer"]

time = np.arange(0,30,0.001)
model = models.th_cell_diff

param_names = ["beta_p"]
param_arrays = [np.logspace(1,3,100)]

norm_list = [10]

# =============================================================================
# param scan
# =============================================================================
df_new = multi_param(param_arrays, param_names, time, cond,
                cond_names, norm_list, model = model)

xlabel = r"$\beta_p$"
ylabel = "log2FC"

loc_major = ticker.LogLocator(base = 10.0, numticks = 100)
loc_minor = ticker.LogLocator(base = 10.0, 
                              subs = np.arange(0.1,1,0.1),
                              numticks = 12)

g = sns.relplot(x = "x", y = "ylog", kind = "line", data = df_new, 
                hue = "cond", col = "readout",
                facet_kws = {"margin_titles" : True, "sharex" : False, "sharey" : False},
                legend = "full",
                aspect = 1.2)

for ax in g.axes.flat:
    ax.xaxis.set_major_locator(loc_major)
    ax.xaxis.set_minor_locator(loc_minor)

g.set(xscale = "log", xlabel = xlabel, ylabel = ylabel)

    
g = sns.relplot(x = "x", y = "ylog", kind = "line", data = df_new, 
                hue = "readout", col = "cond",
                facet_kws = {"margin_titles" : True, "sharex" : False, "sharey" : False},
                legend = "full",
                aspect = 1.2)

for ax in g.axes.flat:
    ax.xaxis.set_major_locator(loc_major)
    ax.xaxis.set_minor_locator(loc_minor)
    
g.set(xscale = "log", xlabel = xlabel, ylabel = ylabel)

# =============================================================================
# beta prolif only look at peak
# =============================================================================
df_prolif_only = df_new[(df_new["pname"] == "beta_p") & (df_new["readout"] == "peak")]

g = sns.relplot(x = "x", y = "ylog", kind = "line", data = df_prolif_only, 
                hue = "cond", 
                facet_kws = {"margin_titles" : True, "sharex" : False, "sharey" : False},
                legend = "full",
                aspect = 1.2)


g.set(xlabel = xlabel, ylabel = ylabel, title = "resp. size")

ylim = (0, 10)
xlim = (None, 100)
g.set(xscale = "log", ylim = ylim, xlim = xlim)

ax = g.axes.flat[0]
ax.xaxis.set_major_locator(loc_major)
ax.xaxis.set_minor_locator(loc_minor)

g.savefig("../figures/beta_prolif_peak.pdf")
#ax.xaxis.set_minor_formatter(ticker.NullFormatter())