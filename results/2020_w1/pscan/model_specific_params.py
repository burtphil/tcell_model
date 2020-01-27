#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 12:45:16 2020

@author: burt
this is to anaylze model specific parameters
"""

from tcell_parameters import d_il2, d_il7, d_timer, d_null

import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
#sys.path.append("C:/Users/Philipp/Documents/projects/tcell_model/code")
import matplotlib.pyplot as plt
from test_module import multi_param, array_from_dict
import module_models as models
import numpy as np
import matplotlib.ticker as ticker
import seaborn as sns
import pandas as pd
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

# =============================================================================
# define exp conditions
# =============================================================================

cond2 = [[d_il2], [d_il7], [d_timer]]

cond_names2 = [["il2"], ["carr. cap."], ["timer"]]

time = np.arange(0, 50, 0.01)
model = models.th_cell_diff

param_names2 = [["rate_il2"], ["rate_il7"], ["crit_timer"]]

df_list = []

for cond, cond_names, param_names in zip(cond2, cond_names2, param_names2):

    norm_list = [d[name] for d, name in zip(cond, param_names)]
    param_arrays = [array_from_dict(d, pname) for d, pname in zip(cond, param_names)]
    
    df_new = multi_param(param_arrays, param_names, time, cond,
                    cond_names, norm_list, model = model, adjust_time = False)
    df_list.append(df_new)
    
df = pd.concat(df_list)

# =============================================================================
# plotting
# =============================================================================
loc_major = ticker.LogLocator(base = 10.0, numticks = 100)
loc_minor = ticker.LogLocator(base = 10.0, 
                              subs = np.arange(0.1,1,0.1),
                              numticks = 12)


g = sns.relplot(data = df, x = "xnorm", y = "ylog", kind = "line", hue = "cond",
                col = "readout")

for ax in g.axes.flat:
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(loc_major)
    ax.xaxis.set_minor_locator(loc_minor)  
    plt.setp(ax.texts, text="")

ylabel = "log2FC"
g.set(ylabel = ylabel, xscale = "log", ylim = (0, None))
g.savefig("pscan_model_specific2.pdf")
# =============================================================================
# new plot
# =============================================================================
g = sns.relplot(data = df, x = "xnorm", y = "ylog", kind = "line", hue = "readout",
                col = "cond")

ylabel = "log2FC"
g.set(ylabel = ylabel, xscale = "log", ylim = (0, None))

for ax in g.axes.flat:
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(loc_major)
    ax.xaxis.set_minor_locator(loc_minor)  
    plt.setp(ax.texts, text="")

    
g.savefig("../figures/pscan_model_specific.pdf")
#ax.xaxis.set_minor_formatter(ticker.NullFormatter())
