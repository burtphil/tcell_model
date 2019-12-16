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

savepath = "/home/burt/Documents/projects/2019/tcell_model/results/20191212/"
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
time = np.arange(0,6,0.01)

def make_dicts(beta_p, d):
    
    d_list = [dict(d) for i in beta_p]
    for dic, rate in zip(d_list, beta_p):
        dic["beta_p"] = rate
    return d_list

# =============================================================================
# plot prolif for low rates
# =============================================================================
#prepare dicts for different prolif rates
rates_beta_p = [10, 50, 50]

d_null_prolif = make_dicts(rates_beta_p, d_null)
d_il7_prolif = make_dicts(rates_beta_p, d_il7)
d_il7_prolif[1]["crit_il7"] = 1.5
d_il7_prolif[2]["crit_il7"] = 2.

d_il2_prolif = make_dicts(rates_beta_p, d_il2)
d_il2_prolif[1]["beta_sad"] = 0.1
d_il2_prolif[1]["beta_sad"] = 10

d_timer_prolif = make_dicts(rates_beta_p, d_timer)
d_timer_prolif[1]["crit_timer"] = 1.0
d_timer_prolif[2]["crit_timer"] = 1.4

d_prolif = [d_null_prolif, d_il2_prolif, d_il7_prolif, d_timer_prolif]
labels = [r"$\beta_p< d$",r"$\beta_p>d$",r"$\beta_p=20$"]
titles = ["No Control", "IL2", "IL7", "Timer"]
colors = ["k", "tab:grey", "tab:grey"]
styles = ["-", "-", "--"]

fig, ax = plt.subplots(1,4, figsize = (18,4))

for a, dic, title in zip(ax, d_prolif, titles):
    
    all_cells = [model.get_cells2(time, d) for d in dic]
    teffs = [cells[-1] for cells in all_cells]
    
    for teff, label, c, ls in zip(teffs, labels, colors, styles):
        a.plot(time, teff, c = c, ls = ls, label = label)
    
    a.set_ylim([0,5])
  
    a.set_title(title)
    a.set_xlabel("time")
    a.set_ylabel("cell dens. norm.")
    
ax[0].legend()
    #a.set_ylim([0,1])

plt.tight_layout()
fig.savefig(savepath+"timecourse_prolif.svg")

# = ============================================================================
# heatmaps
# =============================================================================
res = 30
time = np.arange(0,30,0.01)

name2 = "beta_p"
readout_types = [0,1,4,5]
arr2 = np.linspace(10,80,res)

arrays = [np.linspace(0,3,res), 
         np.linspace(0,20,res), 
         np.linspace(0,3,res)]

names = ["crit_timer", "beta_sad", "crit_il7"]
params =  [d_timer, d_il2, d_il7]

def get_heatmap_list(arrays1, arr2, names, name2, time, params, readout_types):
    
    heatmap_list = []
    for arr1, name1, param in zip(arrays1, names, params):
        for readout_type in readout_types:
            heatmap_arr = model.get_heatmap(arr1, arr2, name1, name2, time, param, 
                                            readout_type)
            
            heatmap_list.append(heatmap_arr)
            
    return heatmap_list

heatmap_list = get_heatmap_list(arrays1 = arrays, arr2 = arr2, names = names, 
                                name2 = name2, time = time, params = params, 
                                readout_types = readout_types)
         
def get_heatmap_arr(heatmap_list, model, readout_type):
    idx_dict = {
                "timer" : 0,
                "il2_beta" : 1,    
                "il7" : 2,
                }
    readout_dict = {"area" : 0,
                    "peak" : 1,
                    "tau" : 2,
                    "decay" : 3,
                    }
    idx = idx_dict[model]
    idx1 = readout_dict[readout_type]

    heatmap_arr = heatmap_list[idx*4+idx1]
    return heatmap_arr

# =============================================================================
# plotting parameters
# =============================================================================
labels = ["timer t0", r"$\beta_e$", "crit. IL2", "Carrying Capacity"]
label2 = r"$\beta_p$"

vmin = 0
vmax = 5

model_type = "il2_beta"
readout = "area"
label1 = "rate IL2+ --> IL2-"
area = get_heatmap_arr(heatmap_list, model = model_type, readout_type = readout)   
fig = model.plot_heatmap(area, title = readout, vmin = vmin, vmax = vmax,
                         label1 = label1, label2 = label2)

fig.savefig(savepath+"heatmap_"+readout+"_"+model_type+".pdf")

vmin = 0
vmax = 5

model_type = "il7"
readout = "area"
label1 = "Carrying Capacity"
area = get_heatmap_arr(heatmap_list, model = model_type, readout_type = readout)   
fig = model.plot_heatmap(area, title = readout, vmin = vmin, vmax = vmax,
                         label1 = label1, label2 = label2)

fig.savefig(savepath+"heatmap_"+readout+"_"+model_type+".pdf")

model_type = "timer"
readout = "area"
label1 = "timer $t_0$"
area = get_heatmap_arr(heatmap_list, model = model_type, readout_type = readout)   
fig = model.plot_heatmap(area, title = readout, vmin = vmin, vmax = vmax,
                         label1 = label1, label2 = label2)

fig.savefig(savepath+"heatmap_"+readout+"_"+model_type+".pdf")

vmin = 0
vmax = 1

model_type = "il7"
readout = "tau"
label1 = "Carrying Capacity"
area = get_heatmap_arr(heatmap_list, model = model_type, readout_type = readout)   
fig = model.plot_heatmap(area, title = readout, vmin = vmin, vmax = vmax,
                         label1 = label1, label2 = label2)

fig.savefig(savepath+"heatmap_"+readout+"_"+model_type+".pdf")


model_type = "il7"
readout = "decay"
label1 = "Carrying Capacity"
area = get_heatmap_arr(heatmap_list, model = model_type, readout_type = readout)   
fig = model.plot_heatmap(area, title = readout, vmin = vmin, vmax = vmax,
                         label1 = label1, label2 = label2)

fig.savefig(savepath+"heatmap_"+readout+"_"+model_type+".pdf")

model_type = "timer"
readout = "tau"
label1 = "timer $t_0$"
area = get_heatmap_arr(heatmap_list, model = model_type, readout_type = readout)   
fig = model.plot_heatmap(area, title = readout, vmin = vmin, vmax = vmax,
                         label1 = label1, label2 = label2)

fig.savefig(savepath+"heatmap_"+readout+"_"+model_type+".pdf")

model_type = "timer"
readout = "decay"
label1 = "timer $t_0$"
area = get_heatmap_arr(heatmap_list, model = model_type, readout_type = readout)   
fig = model.plot_heatmap(area, title = readout, vmin = vmin, vmax = vmax,
                         label1 = label1, label2 = label2)

fig.savefig(savepath+"heatmap_"+readout+"_"+model_type+".pdf")