#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:36:25 2020

@author: burt
"""

import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
from tcell_parameters import d_null, d_il7, d_il2, d_timer
import matplotlib.pyplot as plt
from test_module import multi_param
import module_models as models
import numpy as np
import seaborn as sns
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

colors = ["#3498db", "#95a5a6", "#e74c3c", "#34495e"]     
sns.set_palette(colors)     
savedir = "/home/burt/Documents/projects/2019/tcell_model/results/"
date = "20200108/"
savepath = savedir+date
# =============================================================================
# define exp conditions
# =============================================================================

cond = [d_null, d_il7, d_il2, d_timer]

cond_names = ["Null", "Carr. Cap.", "IL2", "Timer"]

time = np.arange(0,30,0.01)

model = models.th_cell_diff

param_names = ["beta", 
               "d_eff", 
               "beta_p"]

param_arrays = [np.arange(5,15,1),
                np.arange(1,2,0.1),
                np.logspace(1,2,50)]

norm_list = [10,
             1.5,
             10]

df_new = multi_param(param_arrays, param_names, time, cond,
                cond_names, norm_list, model = model)


figname = "pscan_readouts"
g = sns.relplot(x = "x", y = "ylog", kind = "line", data = df_new, hue = "cond", 
                col = "readout", row = "pname", height = 4,
                facet_kws = {"margin_titles" : True, "sharex" : False, "sharey" : False},
                legend = "full")


#for ax in g.axes.flat:
#    ax.set_xscale("log")
    
[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
g.savefig(savepath+figname+".pdf")
g.savefig(savepath+figname+".svg")

# =============================================================================
# same plot but now facet for readouts
# =============================================================================
figname = "pscan_models"
g = sns.relplot(x = "x", y = "ylog", kind = "line", data = df_new, hue = "readout", 
                col = "cond", row = "pname", height = 4,
                facet_kws = {"margin_titles" : True, "sharex" : False, "sharey" : False},
                legend = "full")


#for ax in g.axes.flat:
#    ax.set_xscale("log")
    
[plt.setp(ax.texts, text="") for ax in g.axes.flat]
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
g.savefig(savepath+figname+".pdf")
g.savefig(savepath+figname+".svg")