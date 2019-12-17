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
import matplotlib
from matplotlib.ticker import ScalarFormatter, NullFormatter
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
        "rate_il2" : 1000.0,
        "beta_sad" : 0,
        "fb_il2" : 1.0,
        "alpha_IL2" : 5,
        "decay_p" : 1.,
        }

d_il2 = dict(d_null)
d_timer = dict(d_null)
d_il7 = dict(d_null)

d_il2["mode"] = "il2"
d_timer["mode"] = "timer"
d_il7["mode"] = "il7"

dicts = [d_null, d_il7, d_il2, d_timer]
labels = ["Null", "IL7", "IL2", "Timer"]

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
# plot time course for multple prolif rates all models
# =============================================================================

#prepare dicts for different prolif rates
rates_beta_p = np.arange(10,40,2)


norm = matplotlib.colors.Normalize(
    vmin=np.min(rates_beta_p),
    vmax=np.max(rates_beta_p))

# choose a colormap
c_m = matplotlib.cm.Greys

# create a ScalarMappable and initialize a data structure
s_m = matplotlib.cm.ScalarMappable(cmap=c_m, norm=norm)
s_m.set_array([])


d_null_prolif = make_dicts(rates_beta_p, d_null)
d_il7_prolif = make_dicts(rates_beta_p, d_il7)
d_il2_prolif = make_dicts(rates_beta_p, d_il2)
d_timer_prolif = make_dicts(rates_beta_p, d_timer)

d_prolif = [d_null_prolif, d_il7_prolif, d_il2_prolif, d_timer_prolif]
titles = ["No Control", "IL7","IL2", "Timer"]

fig, ax = plt.subplots(1,4, figsize = (18,4))
for a, dic, title in zip(ax, d_prolif, titles):
    
    all_cells = [model.get_cells2(time, d) for d in dic]
    teffs = [cells[-1] for cells in all_cells]
    
    for teff, col in zip(teffs, rates_beta_p):
        c = s_m.to_rgba(col)
        a.plot(time, teff, c = c)
    
    a.set_ylim([0,4])
    a.set_xlim([0,time[-1]])
    a.set_title(title)
    a.set_xlabel("time")
    a.set_ylabel("cell dens. norm.")
    
#ax[0].legend()
    #a.set_ylim([0,1])
#ax[0].set_ylim([0,10])
plt.colorbar(s_m)
plt.tight_layout()

fig.savefig(savepath+"timecourse_prolif.pdf")

# =============================================================================
# plot prolif vs readout models in one plot
# =============================================================================

prolif_arr = np.logspace(1,2,50)
dicts = dicts
ids = [0]
time = np.arange(0,30,0.01)
# preo normalize
cells0 = [model.get_cells2(time, dic) for dic in dicts]
teff0 = [cell[-1] for cell in cells0]
reads0 = [model.get_readouts(teff, time) for teff in teff0]
reads0 = np.asarray(reads0)
reads0 = reads0[:, ids]

# vary param
reads = [model.vary_param(prolif_arr, "beta_p", ids, time, dic) for dic in dicts]
# normalize
reads = np.asarray(reads)
for j in range(len(dicts)):
    reads[j] = np.log2(reads[j]/reads0[j])

fig, ax = plt.subplots(figsize = (5,4))
ax.set_xscale("log")
ax.plot(prolif_arr, reads[0,:], c = "k", label = "Null")
ax.plot(prolif_arr, reads[1,:], c = "tab:grey", label = "K")
ax.plot(prolif_arr, reads[2,:], c = "tab:red", label = "IL2")
ax.plot(prolif_arr, reads[3,:], c = "tab:cyan", label = "Timer")
ax.set_ylim([0,4.5])
ax.set_xlim([15,75])
ax.set_xlabel(r"$\beta_p$")
ax.set_ylabel("log2FC")
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks(np.arange(20,80,10))
ax.set_title("peak max.")
ax.legend()
plt.tight_layout()
fig.savefig(savepath+"rate_prolif_models_peak.svg")