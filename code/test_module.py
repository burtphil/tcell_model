#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 16:41:01 2019

@author: burt

new module
"""

import numpy as np
from scipy.integrate import odeint
import pandas as pd
from module_models import th_cell_diff, branch_precursor, branch_competetive
from module_readouts import get_area, get_decay, get_peaktime
from functools import partial
# =============================================================================
# 
# =============================================================================
def get_cells(state, d, model):
    
    if model == th_cell_diff:
        tnaive = state[:, :-(2*d["alpha_p"])]
        teff = state[:, -(2*d["alpha_p"]):]
        
        tnaive = np.sum(tnaive, axis = 1)
        teff = np.sum(teff, axis = 1)
        cells = np.stack((tnaive, teff), axis = -1)

    if model == branch_competetive:
        th1 = state[:, 1:(d["alpha1"]+d["alpha1_p"])]
        th2 = state[:, (d["alpha1"]+d["alpha1_p"]):]
        
        th1 = th1[:,-d["alpha1_p"]:]
        th2 = th2[:,-d["alpha2_p"]:]    
        
        th1 = np.sum(th1, axis = 1)
        th2 = np.sum(th2, axis = 1)     
        cells = np.stack((th1, th2), axis = -1)
        
    if model == branch_precursor:
        th1 = state[:, 1:(d["alpha1"]+d["alpha1_p"]+1)]
        th2 = state[:, (d["alpha1"]+d["alpha1_p"]+1):]  
        
        th1 = th1[:,-d["alpha1_p"]:]
        th2 = th2[:,-d["alpha2_p"]:]      
        
        th1 = np.sum(th1, axis = 1)
        th2 = np.sum(th2, axis = 1)
        cells = np.stack((th1, th2), axis = -1)
        
    return cells

def init_model(d, model, initial_cells):
    
    if model == th_cell_diff:
        y0 = np.zeros(d["alpha"]+2*d["alpha_p"])
    if model == branch_precursor:
        y0 = np.zeros(d["alpha1"]+d["alpha1_p"]+d["alpha2"]+d["alpha2_p"]+1)
    if model == branch_competetive:
        y0 = np.zeros(d["alpha1"]+d["alpha1_p"]+d["alpha2"]+d["alpha2_p"]-1)
             
    y0[0] = initial_cells

    return y0

def run_model(time, d, model, initial_cells):
    d = dict(d)
    y0 = init_model(d, model, initial_cells)

    state = odeint(model, y0, time, args = (d,)) 
    state = get_cells(state, d, model)

    return state

def run_exp(time, cond, cond_names, model = th_cell_diff, initial_cells = 1):
    if model == th_cell_diff:
        cell_names = ["naive", "teff"]
    else:
        cell_names = ["Th1", "Tfh"]
    states = [run_model(time, d, model, initial_cells) for d in cond]
    df = df_from_exp(time, states, cond_names, cell_names)
    return df 

def df_from_exp(time, states, cond_names, cell_names):
    df_list = [sim_to_df(time, state, name, cell_names) for state, name in zip(states, cond_names)]
    df = pd.concat(df_list)
    
    return df

def sim_to_df(time, state, name, cell_names):
    df = pd.DataFrame(state, columns = cell_names)
    df["time"] = time
    df["cond"] = name
    
    id_vars = ["time", "cond"]

    df = pd.melt(df, id_vars = id_vars, var_name = "cell")
    
    return df

def generate_readouts(df, time):
    groups = ["cond", "cell"] 
    
    get_peaktime_partial = partial(get_peaktime, time)
    get_area_partial = partial(get_area, time)
    get_decay_partial = partial(get_decay, time)
    
    df = df.groupby(groups).agg(
            peak = ("value", max),
            tau = ("value", get_peaktime_partial),
            area = ("value", get_area_partial),
            decay = ("value", get_decay_partial),            
            )
    
    return df

def vary_param(param_arr, param_name, time, cond, cond_names, norm):
    # change cond dicts for each parameter
    df_arr = [get_tidy_readouts(p, param_name, time, cond, cond_names) for p in param_arr]
    df = pd.concat(df_arr)

    # get values for norm conditions and merge to original data frame
    df_norm = df[df["x"] == norm]
    # kick out x val columns
    df_norm = df_norm.iloc[:,:-1]
    df_norm = df_norm.rename(columns = {"y":"norm"})
    df = pd.merge(left = df, right = df_norm, how = "left", 
                   left_on = ["cond", "cell", "readout"], 
                   right_on = ["cond", "cell", "readout"])
    
    df["ylog"] = np.log2(df["y"]/df["norm"])
    
    return df

def update_dict(d, val, name):
    
    d = dict(d)
    d[name] = val
    
    return d
    
def get_tidy_readouts(p, param_name, time, cond, cond_names):
    cond = [update_dict(d, p, param_name) for d in cond]
    
    df = run_exp(time, cond, cond_names)
    df2 = generate_readouts(df, time)
    df2 = pd.melt(df2.reset_index(), 
                  id_vars = ["cond", "cell"],
                  value_vars = ["area", "peak", "tau", "decay"], 
                  var_name = "readout",
                  value_name = "y")
    df2["x"] = p

    return df2