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

savepath = "/home/burt/Documents/projects/2019/tcell_model/results/20191216/"
# =============================================================================
# parameters
# =============================================================================
d_null = {
        "b" : 0,
        "alpha" : 10,
        "beta" : 10.,
        "alpha_p" : 10,
        "beta_p" : 10.,
        "d_eff" : 1.5,
        "d_prec" : 0,
        "rate_ifn" : 1.0,
        "K_ifn" : 1.0,
        "fb_ifn" : 0,
        "crit" : False,
        "t0" : None,
        "mode" : "Null",
        "crit_timer" : 1.5,
        "crit_il7" : 2.5,
        "crit_il2" : 0.2,
        "rate_il2" : 100.0,
        "beta_sad" : 0,
        "fb_il2" : 1.0,
        "alpha_IL2" : 5,
        "decay_p" : 1.,
        }

d_il2_timer = dict(d_null)
d_il2_timer["mode"] = "il2_timer"

# =============================================================================
# 
# =============================================================================

def get_heatmap_list(arr1, arr2, name1, name2, time, params, beta_p_arr, readout_types):
    
    heatmap_list = []
    for readout in readout_types:
        for beta_p in beta_p_arr:
            params["beta_p"] = beta_p
            heatmap_arr = model.get_heatmap(arr1, arr2, name1, name2, time, params, 
                                            readout)
            
            heatmap_list.append(heatmap_arr)
            
    return heatmap_list

res = 4
time = np.arange(0,30,0.01)

name1 = "rate_il2"
name2 = "crit_timer"
readout_types = [0,1]

arr1 = np.linspace(0,30,res)
arr2 = np.linspace(0,30,res)

# get heatmaps for different values of beta_p and for all readouts

# I need to normalize to same Null model
beta_p_arr = [25., 80.]
