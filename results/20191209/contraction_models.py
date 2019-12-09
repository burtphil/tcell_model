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

savepath = "/home/burt/Documents/projects/2019/tcell_model/results/20191209/"
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
        "d_eff" : 1.,
        "d_prec" : 0,
        "mode" : "Null",
        "rate_il2" : 1.0,
        "fb_il2" : 1.0,
        "rate_ifn" : 1.0,
        "K_ifn" : 1.0,
        "fb_ifn" : 0,
        "crit" : False,
        "crit_timer" : 2.,
        "crit_il7" : 10,
        "crit_il2" : 0.5,
        "t0" : None,
        "decay_p" : 1.,
        }

d_il2 = dict(d_null)
d_timer = dict(d_null)
d_il7 = dict(d_null)

d_il2["mode"] = "il2"
d_il2["beta_sad"]= 2.

d_timer["mode"] = "timer"
d_il7["mode"] = "il7"

# =============================================================================
# time course
# =============================================================================
dicts = [d_null, d_il7, d_il2, d_timer]
labels = ["Null", "IL7", "IL2", "Timer"]
#dicts = [d_il2]
time = np.arange(0,20,0.01)


#model.plot_timecourse(time, dicts, labels = labels)

# =============================================================================
# compare feedback effect in Null model for rtm and ssm
# =============================================================================
time = np.arange(0,10,0.01)


d_rtm = dict(d_null)
d_ssm = dict(d_rtm)
d_ssm["alpha"] = 1
d_ssm["beta"] = 1.

dicts1 = [d_rtm, d_ssm]
labels1 = ["SD $\Psi_0 = 0.1$", "SD $\Psi_0 = 1.0$"]
model.plot_timecourse(time, dicts1, labels = labels1, filename = "timecourse1")

d_rtm_p = dict(d_rtm)
d_rtm_p["alpha_p"] = 1
d_rtm_p["beta_p"] = 1.

dicts2 = [d_rtm, d_rtm_p]
labels2 = ["SD $\Psi_p = 0.1$", "SD $\Psi_p = 1.0$"]
model.plot_timecourse(time, dicts2, labels = labels2, filename = "timecourse2")


time = np.arange(0,30,0.01)

ids = [0,1,4,5]
param_names = ["fb_ifn", "beta"]
titles = labels1
p_labels = [r"fb strength $\beta_0$(IFN)", "beta"]
filename = savepath+"feedback1"
#ylim = [-0.5,0.5]
model.param_scan(param_names, dicts1, titles, ids, time, filename, p_labels, 
                 ylim = None)
# =============================================================================
# parameter scans, adjusted for ylim 
# =============================================================================
ids = [0,1,4,5]
titles = ["Null", "IL7", "IL2", "Timer"]
param_names = ["SD", "SD_p"]
p_labels = ["SD $\Psi$", "SD $\Psi_p$"]
filename = savepath+"distributions"
#ylim = [-3,3]
model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
                 ylim = None)

param_names = ["beta", "beta_p"]
p_labels = [r"$\beta$", r"$\beta_p$"]
filename = savepath+"proliferation"
ylim = [0,10]
model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
                 ylim = ylim)

param_names = ["beta", "d_eff"]
p_labels = [r"$\beta$", "death rate"]
filename = savepath+"rate_parameters"
#ylim = [-2,2]
model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
                 ylim = None)

param_names = ["fb_ifn", "beta"]
p_labels = [r"fb strength $\beta_0$(IFN)", "beta"]
filename = savepath+"feedback"
#ylim = [-0.5,0.5]
#model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
#                 ylim)

# =============================================================================
# effect of model specific parameters
# for this, put model in unstable state by increasing beta p
# =============================================================================
d_null = {
        "b" : 0,
        "alpha" : 10,
        "beta" : 10,
        "alpha_p" : 10,
        "beta_p" : 20,
        "beta_sad" : 0,
        "d_eff" : 1,
        "d_prec" : 0,
        "mode" : "Null",
        "rate_il2" : 1.0,
        "fb_il2" : 1.0,
        "rate_ifn" : 1.0,
        "K_ifn" : 1.0,
        "fb_ifn" : 0,
        "crit" : False,
        "crit_timer" : 2,
        "crit_il7" : 10,
        "crit_il2" : 0.5,
        "t0" : None,
        "decay_p" : 1.,
        }

d_il2 = dict(d_null)
d_timer = dict(d_null)
d_il7 = dict(d_null)

d_il2["mode"] = "il2"
d_il2["beta_sad"]= 1

d_timer["mode"] = "timer"
d_il7["mode"] = "il7"

dicts = [d_null, d_il7, d_il2, d_timer]

param_names = ["crit_timer"]
p_labels = ["timer t0"]
filename = savepath+"crit_timer"
titles = ["Timer"]
dicts = [d_timer]
ylim = [-3,3]


model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
                 ylim)

param_names = ["crit_il7"]
p_labels = ["carrying capacity"]
filename = savepath+"crit_il7"
titles = ["IL7 restriction"]
dicts = [d_il7]
ylim = [-3,3]
model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
                 ylim)

param_names = ["rate_il2"]
p_labels = ["IL2 secretion rate"]
filename = savepath+"crit_il2"
titles = ["IL2 restriction"]
dicts = [d_il2]
ylim = [-3,3]
model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
                 ylim)