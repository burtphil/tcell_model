#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:59:03 2020

@author: burt
vary model specific parameters and check how much I have to tune other parameters
to keep resp. size the same
"""

import numpy as np
import seaborn as sns

from parameters import d
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/ode_models/")

from analysis_module import multi_param, param_norm, arr_from_d, update_dict, vary_param
from ode_models import il2_model, timer_model, null_model, C_model

sns.set(context = "talk", style = "ticks", rc = {"lines.linewidth": 4})
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import pandas as pd
# =============================================================================
# normalize for gamma
# =============================================================================

time = np.arange(0,200,0.001)

arr_model_name = "r_myc"
arr_timer = arr_from_d(d, arr_model_name)

pname = "gamma"
guess_range = (0.1,4)

norm = 40.
model = timer_model

arr_gamma_timer = param_norm(arr_timer, arr_model_name, pname, 
                             guess_range, time, d, model, norm)


# =============================================================================
# normalize for C model
# =============================================================================

time = np.arange(0,500,0.001)

arr_model_name = "r_C"
arr_C = arr_from_d(d, arr_model_name)
pname = "gamma"
guess_range = (0.8,4)

norm = 40.
model = C_model

arr_gamma_C = param_norm(arr_C, arr_model_name, pname, 
                             guess_range, time, d, model, norm)

# =============================================================================
# normalize for IL2 model
# =============================================================================

time = np.arange(0,200,0.001)

arr_model_name = "r_il2"
arr_il2 = arr_from_d(d, arr_model_name)

pname = "gamma"
guess_range = (0.1,4)

norm = 40.
model = il2_model

arr_gamma_IL2 = param_norm(arr_il2, arr_model_name, pname, 
                             guess_range, time, d, model, norm)

df_il2 = pd.DataFrame({"x" : arr_il2, "y" : arr_gamma_IL2, 
                       "model" : "IL2"})

df_timer = pd.DataFrame({"x" : arr_il2, "y" : arr_gamma_timer, 
                        "model" : "Timer"})

df_C = pd.DataFrame({"x" : arr_il2, "y" : arr_gamma_C, "model" : "K"})

df = pd.concat([df_il2, df_timer, df_C])

g = sns.relplot(data = df, x = "x", y = "y", hue = "model")
g.set(xscale = "log", xlabel = "model param (norm.)", ylabel = r"$\gamma$")
g.savefig("../figures/norm_gamma.pdf")

# =============================================================================
# check how readouts change for the models with new param dicts
# =============================================================================

# generate dic_lists
pname = "gamma"
dic_list_il2 = [update_dict(d, a, pname) for a in arr_gamma_IL2]
dic_list_timer = [update_dict(d, a, pname) for a in arr_gamma_timer]
dic_list_C = [update_dict(d, a, pname) for a in arr_gamma_C]

# vary r myc and get readouts
pname = "r_myc"
arr = arr_timer
model_list = [timer_model]
model_names = ["Timer"]
dic_list = dic_list_timer

df_rmyc_reads = vary_param(arr, pname, d, time, model_list, model_names, dic_list)

# vary r il2 and get readouts
pname = "r_il2"
arr = arr_il2
model_list = [il2_model]
model_names = ["IL2"]
dic_list = dic_list_il2

df_il2_reads = vary_param(arr, pname, d, time, model_list, model_names, dic_list)


# vary r C and get readouts
pname = "r_C"
arr = arr_C
model_list = [C_model]
model_names = ["K"]
dic_list = dic_list_C

df_C_reads = vary_param(arr, pname, d, time, model_list, model_names, dic_list)
#df_C_reads = df_C_reads[df_C_reads.y.notnull()]
# plot this
df_reads = pd.concat([df_rmyc_reads, df_il2_reads, df_C_reads])

g = sns.relplot(data = df_reads, x = "sim", y = "ylog", hue = "readout", col = "model",
                facet_kws=dict(sharex=False))

g.set(xscale = "log", ylabel = "log2FC")

xlabels = ["decay myc", "rate il2", "rate K"]
for ax, xlabel in zip(g.axes.flat, xlabels):
    ax.set_xlabel(xlabel)
    if xlabel == "rate K":
        ax.set_xlim([0.01,1.0])

g.savefig("../figures/norm_gamma_reads.svg")












