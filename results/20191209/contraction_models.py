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
# compare feedback effect in Null model for rtm and ssm
# =============================================================================
time = np.arange(0,10,0.01)

d_rtm = dict(d_null)
d_ssm = dict(d_rtm)
d_ssm["alpha"] = 1
d_ssm["beta"] = 1.

dicts1 = [d_rtm, d_ssm]
labels1 = ["SD $\Psi_0 = 0.1$", "SD $\Psi_0 = 1.0$"]
#model.plot_timecourse(time, dicts1, labels = labels1, filename = "timecourse1")

d_rtm_p = dict(d_rtm)
d_rtm_p["alpha_p"] = 1
d_rtm_p["beta_p"] = 1.

dicts2 = [d_rtm, d_rtm_p]
labels2 = ["SD $\Psi_p = 0.1$", "SD $\Psi_p = 1.0$"]
#model.plot_timecourse(time, dicts2, labels = labels2, filename = "timecourse2")

# =============================================================================
# show time courses for different proliferation rates in the models
# =============================================================================
time = np.arange(0,10,0.01)

def make_dicts(beta_p, d):
    
    d_list = [dict(d) for i in beta_p]
    for dic, rate in zip(d_list, beta_p):
        dic["beta_p"] = rate
    return d_list

# =============================================================================
# plot prolif for low rates
# =============================================================================
#prepare dicts for different prolif rates
rates_beta_p = [10,15,20]

d_null_prolif = make_dicts(rates_beta_p, d_null)
d_il2_prolif = make_dicts(rates_beta_p, d_il2)
d_il7_prolif = make_dicts(rates_beta_p, d_il7)
d_timer_prolif = make_dicts(rates_beta_p, d_timer)

d_prolif = [d_null_prolif, d_il7_prolif, d_il2_prolif, d_timer_prolif]
labels = [r"$\beta_p=10$",r"$\beta_p=15$",r"$\beta_p=20$"]
titles = ["Null", "IL7", "IL2", "Timer"]
colors = ["k", "tab:grey", "tab:cyan", "tab:purple"]

fig, ax = plt.subplots(1,4, figsize = (12,4))

for a, dic, title in zip(ax, d_prolif, titles):
    
    all_cells = [model.get_cells2(time, d) for d in dic]
    teffs = [cells[-1] for cells in all_cells]
    
    for teff, label, c in zip(teffs, labels, colors):
        a.plot(time, teff, c = c, label = label)
    
    a.set_ylim([0,0.8])
  
    a.set_title(title)
    a.set_xlabel("time")
    a.set_ylabel("cell dens. norm.")
ax[0].legend()
    #a.set_ylim([0,1])

plt.tight_layout()
fig.savefig(savepath+"timecourse_prolif_low.pdf")

# =============================================================================
# plot prolif for higher rates
# =============================================================================
rates_beta_p = [20,30,40]

d_null_prolif = make_dicts(rates_beta_p, d_null)
d_il2_prolif = make_dicts(rates_beta_p, d_il2)
d_il7_prolif = make_dicts(rates_beta_p, d_il7)
d_timer_prolif = make_dicts(rates_beta_p, d_timer)

d_prolif = [d_null_prolif, d_il7_prolif, d_il2_prolif, d_timer_prolif]
labels = [r"$\beta_p=20$",r"$\beta_p=30$",r"$\beta_p=40$"]
titles = ["Null", "IL7", "IL2", "Timer"]
colors = ["k", "tab:grey", "tab:cyan", "tab:purple"]


fig, ax = plt.subplots(1,4, figsize = (12,4))

for a, dic, title in zip(ax, d_prolif, titles):
    
    all_cells = [model.get_cells2(time, d) for d in dic]
    teffs = [cells[-1] for cells in all_cells]
    
    for teff, label, c in zip(teffs, labels, colors):
        a.plot(time, teff, c = c, label = label)
    
    a.set_ylim([0,18])
  
    a.set_title(title)
    a.set_xlabel("time")
    a.set_ylabel("cell dens. norm.")
ax[0].legend()
    #a.set_ylim([0,1])

plt.tight_layout()
fig.savefig(savepath+"timecourse_prolif_high.pdf")

# =============================================================================
# show feedback effect for response time and single step model in one plot
# =============================================================================
dicts = [d_null, d_il7, d_il2, d_timer]
labels = ["Null", "IL7", "IL2", "Timer"]
time = np.arange(0,30,0.01)

ids = [0,1,4,5]
titles = ["SD $\Psi_0 = 0.1$", "SD $\Psi_0 = 1.0$"]
p_labels = [r"fb strength $\beta_0$(IFN)"]
filename = savepath+"feedback1"
ylim = [-0.5,0.5]
arr = np.linspace(0,10,30)
param_names = ["fb_ifn"]
# get normal cell numbers for normalization
#model.param_scan(param_names, dicts1, titles, ids, time, filename, p_labels, 
#                 ylim = ylim)
# =============================================================================
# parameter scans, adjusted for ylim 
# =============================================================================
ids = [0,1,4,5]
titles = ["Null", "IL7", "IL2", "Timer"]
param_names = ["SD_p"]
p_labels = ["SD $\Psi_p$"]
filename = savepath+"SD_proliferation"
ylim = [-0.3,1.4]
#model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
#                 ylim = ylim)

param_names = ["SD"]
p_labels = ["SD $\Psi$"]
filename = savepath+"SD_differentiation"
ylim = [-2.2,0.5]
#model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
#                 ylim = ylim)

param_names = ["beta_p"]
p_labels = [r"$\beta_p$"]
filename = savepath+"rate_proliferation"
ylim = [0,10]
#model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
#                 ylim = ylim)

param_names = ["beta", "d_eff"]
p_labels = [r"$\beta$", "death rate"]
filename = savepath+"rate_parameters"
ylim = [-2,2]
#model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
#                 ylim = ylim)

param_names = ["fb_ifn", "beta"]
p_labels = [r"fb strength $\beta_0$(IFN)", "beta"]
filename = savepath+"feedback"
ylim = [-0.5,0.5]
#model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
#                 ylim)

# =============================================================================
# effect of model specific parameters
# for this, put model in unstable state by increasing beta p
# =============================================================================

param_names = ["crit_timer"]
p_labels = ["timer t0"]
filename = savepath+"crit_timer"
titles = ["Timer"]
dicts = [d_timer]
ylim = [-1,1]
#model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
#                 ylim)

param_names = ["crit_il7"]
p_labels = ["carrying capacity"]
filename = savepath+"crit_il7"
titles = ["IL7 restriction"]
dicts = [d_il7]
ylim = [-1,1]
#model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
#                 ylim)

param_names = ["crit_il2"]
p_labels = ["IL2 secretion rate"]
filename = savepath+"crit_il2"
titles = ["IL2 restriction"]
dicts = [d_il2]
ylim = [-1,1]
#model.param_scan(param_names, dicts, titles, ids, time, filename, p_labels, 
#                 ylim)

# = ============================================================================
# heatmaps
# =============================================================================
res = 30

savepath = savepath + "/heatmaps/"
label2 = r"$\beta_p$"
name2 = "beta_p"
readout_types = [0,1,4,5]
arr2 = np.linspace(10,30,res)
vmin = -5
vmax = 5

arrays = [np.linspace(0,30,res), 
         np.linspace(0,2,res), 
         np.linspace(0,1,res), 
         np.linspace(5,15,res)]

names = ["crit_timer", "beta_sad", "crit_il2", "crit_il7"]
params =  [d_timer, d_il2, d_il2, d_il7]
labels = ["timer t0", r"$\beta_e$", "crit. IL2", "Carrying Capacity"]
filenames = ["timer", "il2_beta", "il2_crit", "il7"]

def plot_heatmaps(arrays1, arr2, names, name2, time, params, labels, vmin, vmax, 
                  readout_types, savepath, filenames, label2):
    
    xlist = ["peak", "area", "decay", "tau"]
    
    for arr1, name1, param, label, filename in zip(arrays1, names, params, labels, filenames):
        for readout_type, x in zip(readout_types, xlist):
            fig = model.get_heatmap(arr1, arr2, name1, name2, time, param, 
                                    readout_type, vmin = vmin, vmax = vmax, 
                                    label1 = label, label2 = label2)
            fig.savefig(savepath+"heatmap"+"_"+x+"_"+filename+".pdf")

plot_heatmaps(arrays1 = arrays, arr2 = arr2, names = names, name2 = name2, time = time,
              params = params, labels = labels, vmin = vmin, vmax = vmax,
              readout_types = readout_types, savepath = savepath, filenames = filenames,
              label2 = label2)
            
# =============================================================================
# make heatmap to compare IL2 parameters to 
# =============================================================================
name1 = "crit_il2"
name2 = "crit_il7"
label1 = "beta e"
label2 = "C"
res = 30
arr1 = np.linspace(0,2,res)
arr2 = np.linspace(1,100,res)
readout_type = 1

time = np.arange(0,30,0.1)

# make dict that is affected by both IL7 and IL2
d_il27 = dict(d_il2)
d_il27["mode"] = "il2+"
d_il27_prolif_high = dict(d_il27)
d_il27_prolif_high["beta_p"] = 40

#model.plot_timecourse(time, dicts = [d_il27_prolif_high], labels = ["hi"])

title = r"$\beta_p=10$"
vmin = -0,3
vmax = 0.3
model.get_heatmap(arr1, arr2, name1, name2, time, d_il27, readout_type, vmin, vmax,
                  title, label1, label2)

vmin = -5
vmax = 5
title = r"$\beta_p=40$"
model.get_heatmap(arr1, arr2, name1, name2, time, d_il27_prolif_high, readout_type, vmin, vmax,
                  title, label1, label2)