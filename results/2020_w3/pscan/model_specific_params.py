#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 12:45:16 2020

@author: burt
this is to anaylze model specific parameters
"""

import numpy as np
import seaborn as sns

from parameters import d
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/ode_models/")

from analysis_module import multi_param, arr_from_d, vary_param
from ode_models import il2_model, timer_model, null_model, C_model

sns.set(context = "talk", style = "ticks", rc = {"lines.linewidth": 4})
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import pandas as pd

# =============================================================================
# define exp conditions
# =============================================================================

model_list2 = [[il2_model], [C_model], [timer_model]]
model_names2 = [["il2"], ["K"], ["timer"]]
pnames = ["r_il2", "r_C", "r_myc"]

time = np.arange(0, 1000, 0.001)

df_list = []

for model_list, model_names, pname in zip(model_list2, model_names2, pnames):

    arr = arr_from_d(d, pname)    
    df_new = vary_param(arr, pname, d, time, model_list, model_names)
    
    df_list.append(df_new)
    
df = pd.concat(df_list)

# =============================================================================
# plotting
# =============================================================================
loc_major = ticker.LogLocator(base = 10.0, numticks = 100)
loc_minor = ticker.LogLocator(base = 10.0, 
                              subs = np.arange(0.1,1,0.1),
                              numticks = 12)


g = sns.relplot(data = df, x = "sim", y = "ylog", kind = "line", hue = "readout",
                col = "pname")

for ax in g.axes.flat:
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(loc_major)
    ax.xaxis.set_minor_locator(loc_minor)  
    plt.setp(ax.texts, text="")

ylabel = "log2FC"
g.set(ylabel = ylabel, xscale = "log", ylim = (0, None))
