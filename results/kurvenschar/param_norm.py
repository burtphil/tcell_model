#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 15:13:46 2020

@author: burt
"""

import sys
sys.path.append("C:/Users/Philipp/Documents/projects/tcell_model/code")
from tcell_parameters import d_il2, d_timer, d_il2_timer
import matplotlib.pyplot as plt
from test_module import norm_readout
import module_models as models
import numpy as np
import seaborn as sns
import pandas as pd

sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})
colors = ["#3498db", "#95a5a6", "#e74c3c", "#34495e"]     
sns.set_palette(colors)

time = np.arange(0,50,0.01)
pname = "crit_timer"

cond = d_timer
cond_names = ["timer"]

#out = norm_readout(pname, pmin, pmax, time, cond, cond_names)
#print(out)
arr = np.arange(30.,40,1)


def fun(arr, pname, guess_arr, cond_list, cond_names):
       
    df_list = [] 
    for cond, cond_name, guess_range in zip(cond_list, cond_names, guess_arr):
        df = fun2(arr, pname, guess_range, time, cond, cond_name)
        df_list.append(df)
        
    df_cat = pd.concat(df_list)
        
    return df_cat

def fun2(arr, pname, guess_range, time, cond, cond_name):

    out_arr = np.zeros_like(arr) 

    for i, beta_p in enumerate(arr):
        cond = dict(cond)
        cond[pname] = beta_p                   
        out = norm_readout(cond_name, guess_range, time, cond)
        print(out)
        out_arr[i] = out
        
    df = pd.DataFrame({"x" : arr, "y" : out_arr})
    df["ynorm"] = df.y / float(df.y[df.x == 35.0])
    df["ylog"] = np.log2(df.ynorm)
    df["model"] = cond_name
    
    return df

guess_arr = [(0,5), (0,2)]
pname = "beta_p"
cond_list = [d_timer, d_il2]
cond_names = ["crit_timer","rate_il2"]      

df = fun(arr, pname, guess_arr, cond_list, cond_names)