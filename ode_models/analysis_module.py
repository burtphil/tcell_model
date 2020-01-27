#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:23:41 2020

@author: burt

run ode models
"""
from ode_models import ode_model, core_model_gamma
from scipy.integrate import odeint
import numpy as np
import pandas as pd
from readout_module import get_readouts, get_area
import itertools
import matplotlib.pyplot as plt

def run_exp(time, d, model, sim_name, model_name, core_model = core_model_gamma):
    """
    provide time and sim conditions
    sim_name, model_name are strings, model is a function, d is param dict
    return data frame with simulation
    """
    d = dict(d)
    
    init = d["init"]
    init_myc = d["init_myc"]
    y0 = [init,0,0, init_myc]
    
    state = odeint(ode_model, y0, time, args = (d, model, core_model), hmax = 0.01)
    state = get_cells(state, time, sim_name, model_name)
    
    return state

def param_norm(arr_model, arr_model_name, pname, guess_range, time, d, model, norm):
    """
    arr model is model specific parameter array that is varied along x axis
    the variation in arr_model to keep response size constant   
    """
    d = dict(d)
    arr_out = np.zeros_like(arr_model)
    
    for i, arr_m in enumerate(arr_model):
        d[arr_model_name] = arr_m
        out = norm_readout(pname, guess_range, time, d, model, norm)
        arr_out[i] = out
        
    return arr_out

def multi_exp(time, d_list, model_list, sim_names, model_names, core_model = core_model_gamma):
    """
    run multiple simulations for different models (stored in model_list as list)
    and different conditions, e.g. variation of one parameter stored in d_list
    d_list is list of dicts
    """
    df_list = []
    for d, sim_name in zip(d_list, sim_names):
        for model, model_name in zip(model_list, model_names):
            df = run_exp(time, d, model, sim_name, model_name, core_model)        
            df["sim"] = sim_name
            df["model"] = model_name
            df_list.append(df)
            
    df = pd.concat(df_list)

    return df

def get_heatmap(arr1, arr2, name1, name2, time, d, model, readout_fun, norm = None,
                core_model = core_model_gamma):
    d = dict(d)
    sim_name = "dummy"
    model_name = "dummy"
    
    z  = []
    for val1, val2 in itertools.product(arr1, arr2):
        d[name1] = val1
        d[name2] = val2
        df = run_exp(time, d, model, sim_name, model_name, core_model = core_model)   
        readout = readout_fun(df.time, df.cells)
        if norm != None:
            readout = np.log2(readout/norm)
        z.append(readout)

    #print(len(z))
    z = np.asarray(z)
    z = z.reshape(len(arr1), len(arr2))
    z = z[:-1, :-1]    
    z = z.T
    heatmap_arr = [arr1, arr2, z]
    
    return heatmap_arr

def plot_heatmap(heatmap_arr, vmin = None, vmax = None, title = None, 
                 label1 = None, label2 = None, cmap = "bwr"):

    arr1, arr2, val = heatmap_arr
    fig, ax = plt.subplots(figsize = (6,4))
    color = cmap
    cmap = ax.pcolormesh(arr1, arr2, val, cmap = color, vmin = vmin, vmax = vmax,
                         rasterized = True)
    
    ax.set_xlabel(label1)
    ax.set_ylabel(label2)
    ax.set_title(title)
    cbar = plt.colorbar(cmap)
    cbar.set_label("log2FC")
        
    plt.tight_layout()

    return fig

def get_cells(state, time, sim_name, model_name):
    
    # calculate cells, leave out naive cells (first entry) and myc conc (last entry)
    cells = np.sum(state[:, 1:-1], axis = 1)

    df = pd.DataFrame(cells, columns = ["cells"])
    df["time"] = time

    return df

def generate_readouts(df):
    groups = ["sim", "model"] 
    df = df.groupby(groups).apply(get_readouts)
    
    return df

def update_dict(d, val, name):
    d = dict(d)
    d[name] = val
    
    return d

def multi_param(pnames, pnames_arr, d, time, model_list, model_names, arr = None, core_model = core_model_gamma):
    df_list = []
    
    for pname, arr in zip(pnames, pnames_arr):
        df = vary_param(arr, pname, d, time, model_list, model_names, core_model = core_model)
        df_list.append(df)
    
    df = pd.concat(df_list)
    
    return df

def vary_param(arr, pname, d, time, model_list, model_names, dic_list = None, core_model = core_model_gamma):
    """
    vary parameter provide param array and parameter name
    if dic_list is provided it exchanges d and for each arr. element there should be a d in dic_list
    """
    
    if dic_list != None:
        d_list = [update_dict(d, a, pname) for d, a in zip(dic_list, arr)]
    else:
        d_list = [update_dict(d, a, pname) for a in arr]
    
    cells = multi_exp(time, d_list, model_list, arr, model_names, core_model = core_model)
    readouts = generate_readouts(cells)
    
    df = pd.melt(readouts.reset_index(), 
                  id_vars = ["sim", "model"],
                  value_vars = ["area", "peak", "tau", "decay"], 
                  var_name = "readout",
                  value_name = "y")   
    
    df["pname"] = pname

    # get middle element of param array
    # note that this is only precise middle for odd len(sim_names)
    norm = int(len(arr) / 2)
    norm = arr[norm]
    
    # get arr where sim is middle
    df_norm = df[df["sim"] == norm]
    df_norm = df_norm.rename(columns = {"y":"ynorm"})
    df_norm = df_norm[["model", "readout", "ynorm"]]
    
    # merge to original array
    df = pd.merge(left = df, right = df_norm, how = "left", 
                   left_on = ["model", "readout"], 
                   right_on = ["model", "readout"])

    # add normalization
    logseries = df["y"]/df["ynorm"]
    logseries = logseries.astype(float)
    df["ylog"] = np.log2(logseries)
    
    return df

def arr_from_d(d, pname, res = 25):
    """
    takes dict, returns array for param scan
    with default values of dict * or / 100%
    """
    x = d[pname]
    xmin = x/10.
    xmax = x*10.
    arr1 = np.linspace(xmin, x, res, endpoint = False)
    arr2 = np.linspace(x, xmax, res)
    arr = np.concatenate((arr1, arr2))
    return arr


def norm_readout(pname, guess_range, time, d, model, norm):
    """
    guess range should be a tuple of pmin pmax, 
    cond is a dict with model params
    pname is name of parameter to be normalized
    norm cond: value of area to normalize against
    """
    sim_name = "dummy"
    model_name = "dummy"
    
    pmin, pmax = guess_range
    
    d = dict(d)
    d1 = dict(d)
    d2 = dict(d)
    d1[pname] = pmin
    d2[pname] = pmax
    
    df1 = run_exp(time, d1, model, sim_name, model_name)
    df2 = run_exp(time, d2, model, sim_name, model_name)
    
    area1 = get_area(df1.time, df1.cells)
    area2 = get_area(df2.time, df2.cells)

    #print(pname, cond, area1, area2)
    
    if not(area1 < norm < area2):
        guess = np.nan
        print("param not in guess range or nan")
        return guess
    
    else:
        guess = (pmin+pmax) / 2
        crit = False
        counter = 0
        while crit == False:
            counter = counter + 1
            if counter > 50:
                print("counter reached, stopping normalization")
                break
            
            d[pname] = guess
    
            df = run_exp(time, d, model, sim_name, model_name)
            area = get_area(df.time, df.cells)

            #print("guess...area...pmin...pmax")
            #print(guess, area, pmin, pmax)
            if area < norm:
                pmin = guess
                guess = (guess+pmax)/2
                
            else:
                pmax = guess
                guess = (guess+pmin)/2
                          
            if np.abs(area-norm) < 0.005:
                crit = True
                
        return guess