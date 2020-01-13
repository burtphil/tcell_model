
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 12:36:25 2020

@author: burt
vary beta p norm for area for homeostasis model and then plot time course
"""
import sys

#sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
sys.path.append("C:/Users/Philipp/Documents/projects/tcell_model/code")

from tcell_parameters import d_null, d_il7, d_il2, d_timer
import matplotlib.pyplot as plt
from test_module import run_exp, multi_exp, norm_readout, update_dict, update_dicts
import module_models as models
import numpy as np
import seaborn as sns
import pandas as pd
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})
import matplotlib
# =============================================================================
# define exp conditions
# =============================================================================
model = models.th_cell_diff


time = np.arange(0,200,0.01)

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

arr = np.arange(30.,40,0.5)
guess_arr = [(0,5), (0,1000), (0,5)]
pname = "beta_p"
cond_list = [d_timer, d_il2, d_il7]
cond_names = ["crit_timer","rate_il2", "rate_il7"]      

df = fun(arr, pname, guess_arr, cond_list, cond_names)
g = sns.relplot(x = "x", y = "ylog", data = df, hue = "model")

timer_arr = df.y[df.model == "crit_timer"].array
il2_arr = df.y[df.model == "rate_il2"].array
il7_arr = df.y[df.model == "rate_il7"].array

timer_dicts = [update_dict(d_timer, val, "crit_timer") for val in timer_arr]
il2_dicts = [update_dict(d_il2, val, "rate_il2") for val in il2_arr]
il7_dicts = [update_dict(d_il7, val, "rate_il7") for val in il7_arr]

cond_list = [[d1,d2,d3] for d1,d2,d3 in zip(timer_dicts, il2_dicts, il7_dicts)]
cond_list = [update_dicts(dict_list, val, "beta_p") for dict_list, val in zip(cond_list, arr)]

names = ["crit_timer", "rate_il2", "rate_il7"]
cond_names2 = ["betap"+str(val) for val in arr]

# =============================================================================
# run experiment
# =============================================================================
cond1 = [d_timer, d_il2, d_il7]
cond_names = ["Timer", "IL2", "Carr. Cap."]

time = np.arange(0,7,0.01)

exp = multi_exp(time, cond_list, cond_names, cond_names2)
exp2 = exp.loc[exp["cond2"] == "betap35.0", :]

norm = matplotlib.colors.Normalize(
    vmin=np.min(arr),
    vmax=np.max(arr))

# choose a colormap
cm = matplotlib.cm.Blues

# create a ScalarMappable and initialize a data structure
sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
sm.set_array([])

g = sns.relplot(x = "time", y = "value", kind = "line", data = exp, hue = "cond2", 
                col = "cond", palette = "Blues", height = 5,legend = False)
   
#g.set(ylim = (0, 4), xlim = (0,time[-1]))
ax = g.axes[0][0]
ax.set_ylabel("cell dens. norm.")
g.set_titles("{col_name}")
cbar = g.fig.colorbar(sm, ax = g.axes)
cbar.set_label(r"$\beta_p$")
g.savefig("proliferation_timecourse.pdf")
#g.savefig("proliferation_timecourse.svg")
