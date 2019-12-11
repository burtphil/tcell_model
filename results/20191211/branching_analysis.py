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
d_prec2["p1_def"] = 0.33

cells_prec1 = m_branch.run_model(d_prec1, t, m_branch.branch_precursor)
cells_comp1 = m_branch.run_model(d_comp1, t, m_branch.branch_competetive)
cells_prec2 = m_branch.run_model(d_prec2, t, m_branch.branch_precursor)
cells_comp2 = m_branch.run_model(d_comp2, t, m_branch.branch_competetive)

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
fig.savefig(savepath+"timecourse.svg")