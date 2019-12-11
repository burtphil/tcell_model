#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 09:48:50 2019

@author: burt
"""
import os
os.chdir("/home/burt/Documents/projects/2019/tcell_model/code")

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from params_branching import d_comp, d_prec
import module_branching as m_branch
savepath = "/home/burt/Documents/projects/2019/tcell_model/results/20191211/"


sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

def update_dict(d, update_type = "RTM"):
    d = dict(d)

    x = 10 if update_type == "RTM" else 1
    d["alpha1"] = x
    d["alpha2"] = x
    d["beta1"] = float(x)
    d["beta2"] = float(x)
    return d

#======================================================================
# normal time course
# plot this as relative th1 tfh fraction?
# =============================================================================
t = np.arange(0,10,0.01)

d_prec1 = dict(d_prec)
d_comp1 = dict(d_comp)

# adjust prec one because it has one step more
d_comp1["alpha1"] = 10
d_comp1["alpha2"] = 10

d_comp2 = dict(d_comp1)
d_prec2 = dict(d_prec1)

d_comp2["beta1"] = 5.
d_prec2["beta1"] = 5.
d_prec2["beta2"] = 7.

cells_prec1 = m_branch.run_model(d_prec1, t, fun = m_branch.branch_precursor)
cells_comp1 = m_branch.run_model(d_comp1, t, fun = m_branch.branch_competetive)
cells_prec2 = m_branch.run_model(d_prec2, t, fun = m_branch.branch_precursor)
cells_comp2 = m_branch.run_model(d_comp2, t, fun = m_branch.branch_competetive)

fig, ax = plt.subplots(1,2, figsize = (10,4.5))

ax[0].plot(t, cells_prec1[0], "k")
ax[0].plot(t, cells_prec2[0], "tab:blue")
ax[0].plot(t, cells_prec2[1], "tab:red")

ax[1].plot(t, cells_comp1[0], "k")
ax[1].plot(t, cells_comp2[0], "tab:blue")
ax[1].plot(t, cells_comp2[1], "tab:red")

for a in ax:
    #a.set_xlim(0, t[-1])
    a.set_ylim([0,0.5])
    a.set_xlabel("time")
    a.set_ylabel("cell dens. norm.")
    
ax[0].set_title("precursor model")
ax[1].set_title("competition model")
plt.tight_layout()
#fig.savefig(savepath+"timecourse.svg")
plt.close()

d_fb = dict(d_prec1)
x = 1
d_fb["alpha1"] = x
d_fb["beta1"] = float(x)
d_fb["alpha2"] = x
d_fb["beta2"] = float(x)
d_fb["fb_prob1"] = 1000
cells_fb = m_branch.run_model(d_fb,t)
cells_fb = np.swapaxes(cells_fb, 0,1)
fig,ax = plt.subplots()
ax.plot(t, cells_fb)

# =============================================================================
# analyze feedback in precursor model for fb on prob fb on rate and both for all
# readouts and for diff chain lenghts
# =============================================================================

# make arr with dicts of diff chain length
alpha_arr = np.arange(1,10,1, dtype = int)
name = "SD"
dic_list_alpha = [m_branch.update_dict(val,name,d_prec1) for val in alpha_arr]

# vary fb strength and target either rate or branching prob
arr = np.linspace(0,1000,100)
name = "fb_rate1"
reads = [m_branch.vary_param(arr, name, t, dic) for dic in dic_list_alpha]
reads_norm = [m_branch.run_model(dic, t, output = "readouts") for dic in dic_list_alpha]
reads = [np.log2(read/read_norm) for read, read_norm in zip(reads, reads_norm)]
reads = np.asarray(reads)
reads = np.swapaxes(reads,0,1)
#reads has shape (alphas, fb_arr, readouts)

name2 = "fb_prob1"
reads2 = [m_branch.vary_param(arr, name2, t, dic) for dic in dic_list_alpha]
reads2 = [np.log2(read/read_norm) for read, read_norm in zip(reads2, reads_norm)]
reads2 = np.asarray(reads2)
reads2 = np.swapaxes(reads2,0,1)
#reads has shape (alphas, fb_arr, readouts)

# plot settings
ylim = [-0.6,0.8]
titles = ["rel. area", "rel. peak size", "rel. tau"]
xlabel1 = "fb --> rate"
xlabel2 = "fb --> prob"
ylabel = "log2FC"
figsize = (14,4)
xticks = np.array([0,500,1000])

# reads fb on rate plot
fig, ax = plt.subplots(1,3, figsize = figsize)
for i in range(3):
    ax[i].plot(arr, reads[:,:,i])
    ax[i].set_ylim(ylim)
    ax[i].set_xlabel(xlabel1)
    ax[i].set_ylabel(ylabel)
    ax[i].set_title(titles[i])
    ax[i].set_xticks(xticks)
plt.tight_layout()


fig.savefig(savepath+"fb_rate.pdf")
fig.savefig(savepath+"fb_rate.svg")

# make same thing for reads2
fig, ax = plt.subplots(1,3, figsize = figsize)
for i in range(3):
    ax[i].plot(arr, reads2[:,:,i])
    ax[i].set_ylim(ylim)
    ax[i].set_xlabel(xlabel2)
    ax[i].set_ylabel(ylabel)
    ax[i].set_title(titles[i])
    ax[i].set_xticks(xticks)
    
plt.tight_layout()
fig.savefig(savepath+"fb_prob.pdf")
fig.savefig(savepath+"fb_prob.svg")
