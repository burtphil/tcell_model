#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:36:25 2020

@author: burt
vary rate parameter for branching models
show that in homeostasis models rates affect each other
"""

from tcell_parameters import d_prec, d_prec_il7, d_prec_il2, d_prec_timer
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
import matplotlib.pyplot as plt
from test_module import multi_param
import module_models as models
import numpy as np
import seaborn as sns
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

colors = ["#3498db", "#95a5a6", "#e74c3c", "#34495e"]     
sns.set_palette(colors)     

# =============================================================================
# define exp conditions
# =============================================================================
cond = [d_prec_il7, d_prec_il2, d_prec_timer]
cond_names = ["Carr. Cap.", "IL2", "Timer"]
time = np.arange(0, 6, 0.01)
model = models.branch_precursor

param_names = ["beta1", 
               "p1_def"]

param_arrays = [np.arange(5,15,0.2),
                np.arange(0, 1, 0.05)]

norm_list = [10,
             0.5]

# =============================================================================
# look at prolif rate
# =============================================================================
df_new = multi_param(param_arrays, param_names, time, cond,
                cond_names, norm_list, model = model)

df1 = df_new
#df1 = df_new.loc[df_new["pname"] == r"$\beta_1$", :] 
df1 = df1.loc[df1["readout"] == "area", :]

g = sns.relplot(x = "x", y = "ylog", kind = "line", data = df1, hue = "cell",
                col = "cond", row = "pname", height = 4, legend = "full",
                facet_kws = {"margin_titles" : True, "sharex" : False, "sharey" : True})

ylim = (-2,2)

xlabel1 = r"$\beta_1$"
xlabel2 = "$p_1$"
ylabel = "log2FC"

axes1 = g.axes[0]
axes2 = g.axes[1]
[ax.set_xlabel(xlabel1) for ax in axes1]
[ax.set_xlabel(xlabel2) for ax in axes2]

# =============================================================================
# 
# df1 = df_new.loc[df_new["pname"] == r"$\beta_1$", :]
# 
# #ylim = (0,5)
# 
# xlabel = r"$\beta_1$"
# ylabel = "log2FC"
# 
# g = sns.relplot(x = "x", y = "ylog", kind = "line", data = df1, hue = "cond", 
#                 col = "readout(rel)", height = 4,
#                 facet_kws = {"margin_titles" : True, "sharex" : True, "sharey" : True},
#                 legend = "full")
# 
# 
# #for ax in g.axes.flat:
# #    ax.set_xscale("log")
#     
# [plt.setp(ax.texts, text="") for ax in g.axes.flat]
# g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
# #g.set(ylim = ylim)
# g.set(xlabel = xlabel, ylabel = ylabel)
# #g.legend.set_title(title)
# g.savefig("beta_branching.svg")
# 
# g = sns.relplot(x = "x", y = "ylog", kind = "line", data = df1, hue = "readout(rel)", 
#                 col = "cond", height = 4, col_order = cond_names,
#                 facet_kws = {"margin_titles" : True, "sharex" : True, "sharey" : True},
#                 legend = "full")
# 
# #for ax in g.axes.flat:
# #    ax.set_xscale("log")
#     
# [plt.setp(ax.texts, text="") for ax in g.axes.flat]
# g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
# #g.set(ylim = ylim)
# g.set(xlabel = xlabel, ylabel = ylabel)
# g.savefig("beta_branching_readouts.svg")
# 
# =============================================================================
