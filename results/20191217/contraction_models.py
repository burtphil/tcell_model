#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 12:52:14 2019

@author: burt
"""
import os
savepath = os.getcwd()
os.chdir("/home/burt/Documents/projects/2019/tcell_model/code")
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})
import module_no_branch as model 

# =============================================================================
# parameters
# =============================================================================
d_null = {
        "b" : 0,
        "alpha" : 10,
        "beta" : 10.,
        "alpha_p" : 10,
        "beta_p" : 10.,
        "beta_sad" : 0,
        "d_eff" : 1.5,
        "d_prec" : 0,
        "mode" : "Null",
        "rate_il2" : 1.0,
        "fb_il2" : 1.0,
        "rate_ifn" : 1.0,
        "K_ifn" : 1.0,
        "fb_ifn" : 0,
        "crit" : False,
        "crit_timer" : 2.,
        "crit_il7" : 10.,
        "crit_il2" : 0.5,
        "t0" : None,
        "decay_p" : 1.,
        "death_mode" : False,
        }

d_il2 = dict(d_null)
d_timer = dict(d_null)
d_il7 = dict(d_null)

d_il2["mode"] = "il2"
d_il2["beta_sad"]= 1.

d_timer["mode"] = "timer"
d_il7["mode"] = "il7"

dicts = [d_null, d_il7, d_il2, d_timer]
labels = ["Null", "IL7", "IL2", "Timer"]
# =============================================================================
# time course
# =============================================================================

dicts = [d_il7]
labels = ["IL7"]
time = np.arange(0,10,0.01)
model.plot_timecourse(time, dicts, labels = labels)

# =============================================================================
# show time courses for different proliferation rates in the models
# =============================================================================
time = np.arange(0,30,0.01)

# =============================================================================
# show feedback effect for response time and single step model in one plot
# =============================================================================
dicts = [d_null, d_il7, d_il2, d_timer]
dicts_death = [dict(dic) for dic in dicts]
for dic in dicts_death:
    dic["death_mode"] = True
    
labels = ["Null", "IL7", "IL2", "Timer"]
time = np.arange(0,30,0.01)

ids = [0,1,4,5]
