#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:33:23 2020

@author: burt
"""

import numpy as np
from scipy.integrate import odeint
import pandas as pd
from module_models import th_cell_diff, branch_precursor, branch_competetive
from module_readouts import get_area, get_decay, get_peaktime, get_peak
from functools import partial
import warnings
import seaborn as sns
sns.set(style = "ticks")

def fun(state, time, beta_p, beta_p1, beta, d):
    dt_state = np.zeros_like(state)
    
    dx = dt_state[0]
    dy = dt_state[1]
    dz = dt_state[2]
    
    
    dx = -beta*state[0]
    dy = beta_p1*beta*state[0]-(beta+d)*state[1]
    dz = beta_p1*beta*state[1] + (beta_p-d)*state[2]
    
    dt_state = [dx,dy,dz]
    return dt_state

y0 = [1,0,0]

beta_p = 1
beta_p1 = 2
beta = 1
d = 1

time = np.arange(0,10,0.01)

state = odeint(fun, y0, time, args = (beta_p, beta_p1, beta, d,))

cell_names = ["x", "y", "z"]
df = pd.DataFrame(state, columns = cell_names)
df["time"] = time

id_vars = "time"

df = pd.melt(df, id_vars = id_vars, var_name = "cell")


g = sns.relplot(data = df, x = "time", y = "value", kind = "line", hue = "cell")
#g.set(ylim = (0,20))