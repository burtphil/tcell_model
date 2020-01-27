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

from analysis_module import multi_param, arr_from_d
from ode_models import C_model_menten, il2_model_menten, timer_model_menten, core_model_menten
import pandas as pd
sns.set(context = "talk", style = "ticks", rc = {"lines.linewidth": 4})

# =============================================================================
# define exp conditions
# =============================================================================

model_list = [timer_model_menten, il2_model_menten, C_model_menten]
model_names = ["Timer", "IL2", "K"]

time = np.arange(0,100,0.001)
pnames = ["r_il2", "r_myc", "r_C"]
pnames_arr = [arr_from_d(d, pname) for pname in pnames]

# =============================================================================
# param scan
# =============================================================================
df_new = multi_param(pnames, pnames_arr, d, time, model_list, model_names, core_model = core_model_menten)

df1 = df_new[(df_new["pname"] == "r_il2") & (df_new["model"] == "IL2")]
df2 = df_new[(df_new["pname"] == "r_C") & (df_new["model"] == "K")]
df3 = df_new[(df_new["pname"] == "r_myc") & (df_new["model"] == "Timer")]

df = pd.concat([df1, df2, df3])

g = sns.relplot(x = "sim", y = "ylog", kind = "line", data = df, 
                hue = "readout", col = "pname",
                facet_kws = {"margin_titles" : True, "sharex" : False, "sharey" : True},
                legend = "full", col_order = pnames,
                aspect = 1.0)
    
g.set(ylim = (-5,5))
g.set(xscale = "log")

for ax, xlabel in zip(g.axes.flat, pnames):
    ax.set_xlabel(xlabel)

g.savefig("../figures/pscan_specific.svg")
g.savefig("../figures/pscan_specific.pdf")

g = sns.relplot(x = "sim", y = "ylog", kind = "line", data = df, 
                hue = "readout", col = "pname",
                facet_kws = {"margin_titles" : True, "sharex" : False, "sharey" : True},
                legend = "full", col_order = pnames,
                aspect = 1.0)
    
g.set(xscale = "log")

for ax, xlabel in zip(g.axes.flat, pnames):
    ax.set_xlabel(xlabel)

g.savefig("../figures/pscan_specific_high.svg")
g.savefig("../figures/pscan_specific_high.pdf")
# make a new data frame only with area readout and make a fake column
# because model specific params are not on the same axis
df_area = df[df.readout == "area"]
df_area = df_area[["sim", "model", "ylog"]]

# IL2 x axis
xnew = pnames_arr[0]

# extract each data frame and append the IL2 x axis
df_list = []
for name in model_names:
    df = df_area[df_area.model == name].reset_index(drop = True)
    df = df.assign(sim2 = pd.Series(xnew))
    df_list.append(df)

df = pd.concat(df_list)

# plot all on same axis
g = sns.relplot(x = "sim2", y = "ylog", data = df,
                hue = "model")

g.set(xscale = "log", xlabel = "model param norm.")
g.savefig("../figures/area_model_specific.pdf")