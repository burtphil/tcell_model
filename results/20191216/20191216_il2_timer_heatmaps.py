#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 12:52:14 2019

@author: burt

"""
import os
os.chdir("/home/burt/Documents/projects/2019/tcell_model/code")
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})
import module_no_branch as model 
import itertools

savepath = "/home/burt/Documents/projects/2019/tcell_model/results/20191216/"
# =============================================================================
# parameters
# =============================================================================
d_null = {
        "b" : 0,
        "alpha" : 10,
        "beta" : 10.,
        "alpha_p" : 10,
        "beta_p" : 40.,
        "d_eff" : 1.5,
        "d_prec" : 0,
        "rate_ifn" : 1.0,
        "K_ifn" : 1.0,
        "fb_ifn" : 0,
        "crit" : False,
        "t0" : None,
        "mode" : "Null",
        "crit_timer" : 7,
        "crit_il7" : 2.5,
        "crit_il2" : 0.0001,
        "rate_il2" : 10.0,
        "beta_sad" : 0,
        "fb_il2" : 1.0,
        "alpha_IL2" : 5,
        "decay_p" : 1.,
        }

d_timer = dict(d_null)
d_timer["mode"] = "timer"

d_il2 = dict(d_null)
d_il2["mode"] = "il2"

d_il2_timer = dict(d_null)
d_il2_timer["mode"] = "il2_timer"

# =============================================================================
# test il2 + timer model
# =============================================================================
labels = ["il2+timer"]
time = np.arange(0,30,0.01)
model.plot_timecourse(time, [d_il2_timer], labels = labels)

res = 30

# =============================================================================
# plot heatmaps
# =============================================================================
# =============================================================================
# params = dict(d_timer)
# arr1 = np.linspace(0,20,res)
# arr2 = np.linspace(15,50,res)
# name1 = "crit_timer"
# name2 = "beta_p"
# label1 = "timer t0"
# label2 = r"$\beta_p$"
#  
# readout_type = [0]
# beta_p_arr = [30]
# 
# figname = "timer_prolif.pdf"
# for readout, beta_p in itertools.product(readout_type, beta_p_arr):
#     
#     #params["beta_p"] = beta_p
#     heatmap_arr = model.get_heatmap(arr1, arr2, name1, name2, time, params, readout)
#     fig = model.plot_heatmap(heatmap_arr, vmin = 0, vmax = 10,
#                              label1 = label1, label2 = label2)
#     fig.savefig(savepath+figname)
#     
# params = dict(d_il2)
# arr1 = np.linspace(0,10000000000000,res)
# arr2 = np.linspace(15,50,res)
# name1 = "rate_il2"
# name2 = "beta_p"
# label1 = "IL2 turnover"
# label2 = r"$\beta_p$"
# readout_type = [0]
# beta_p_arr = [30]
# 
# figname = "il2_prolif.pdf" 
# for readout, beta_p in itertools.product(readout_type, beta_p_arr):
#     
#     #params["beta_p"] = beta_p
#     heatmap_arr = model.get_heatmap(arr1, arr2, name1, name2, time, params, readout)
#     fig = model.plot_heatmap(heatmap_arr, vmin = 0, vmax = 10,
#                              label1 = label1, label2 = label2)
#     fig.savefig(savepath+figname)
# 
# =============================================================================
params = dict(d_il2_timer)
arr1 = np.linspace(0,5,res)
arr2 = np.linspace(0,5,res)
name1 = "rate_il2"
name2 = "crit_timer"
label1 = "IL2 turnover"
label2 = "timer t0"
params["beta_p"] = 80

readout_type = [0]
beta_p_arr = [30]
 
figname = "timer_il2_tau.pdf"
for readout, beta_p in itertools.product(readout_type, beta_p_arr):
    
    #params["beta_p"] = beta_p
    heatmap_arr = model.get_heatmap(arr1, arr2, name1, name2, time, params, readout,
                                    normalize = False)
    fig = model.plot_heatmap(heatmap_arr, label1 = label1, label2 = label2)
    fig.savefig(savepath+figname)

