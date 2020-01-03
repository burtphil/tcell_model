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

def run_exp(time, cond, cond_names, model = th_cell_diff, initial_cells = 1,
            keep_naive = False):
    
    if model == th_cell_diff:
        cell_names = ["naive", "teff"]
    else:
        cell_names = ["Th1", "Tfh"]
        
    states = [run_model(time, d, model, initial_cells) for d in cond]
    df = df_from_exp(time, states, cond_names, cell_names)
    
    if keep_naive == False:
        df = df[df["cell"] != "naive"]
        
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

def vary_param(param_arr, param_name, time, cond, cond_names, norm, model = th_cell_diff):
    """
    vary parameter values get readouts for different conditions
    cond: list of parameter dictioniaries for each condition
    cond_names: list of condition names
    param_arr: array of parameter values
    param_name: name of parameter to be varied
    norm: parameter value to normalize against
    ! note that norm value needs to be in param_arr
    returns df with readouts for each cell and condition and parameter value
    """    
    df_arr = [get_tidy_readouts(p, param_name, time, cond, cond_names, model) for p in param_arr]
    df = pd.concat(df_arr)

    # get values for norm condition (all cell types all conditions all readouts) and merge to original data frame
    # check that provided normalization value is in provided array
    assert np.amin(param_arr) <= norm <= np.amax(param_arr)
    # adjust norm to always be in param arr by changing norm to closest matching value in param_arr
    norm = param_arr[np.argmin(np.abs(norm-param_arr))]
    df_norm = df[df["x"] == norm]
    # kick out x val columns
    df_norm = df_norm.iloc[:,:-1]
    df_norm = df_norm.rename(columns = {"y":"ynorm"})

    df = pd.merge(left = df, right = df_norm, how = "left", 
                   left_on = ["cond", "cell", "readout"], 
                   right_on = ["cond", "cell", "readout"])
    
    df["ylog"] = np.log2(df["y"]/df["ynorm"])
    
    pname = convert_name(param_name)
    df["pname"] = pname
    
    return df

def get_relative_readouts(df):
    """
    take dataframe from vary_param and compute relative readouts
    """
    
    #split dataframe to have one df for each cell type
    df1 = df[df["cell"] == "Th1"]
    df2 = df[df["cell"] == "Tfh"]
   
    # use df1 as basis for new df
    df_base = df1[["cond", "readout", "pname", "x"]]
    
    df_base = df_base.reset_index(drop = True)
    df1 = df1[["y", "ynorm"]].reset_index(drop = True)
    df2 = df2[["y", "ynorm"]].reset_index(drop = True)
    
    # divide dfs of cell types with y and ynorm values and compute logFC
    # then combine with base df
    
    df3 = df1/df2 
    df3["ylog"] = np.log2(df3["y"]/df3["ynorm"])
    
    df_base = pd.concat([df_base, df3], axis = 1)
    df_base = df_base.rename(columns = {"readout":"readout(rel)"})
    
    return df_base

def multi_param(param_arrays, param_names, time, cond,
                cond_names, norm_list, model = th_cell_diff):
    
    df_arr =[vary_param(param_arr, param_name, time, cond, cond_names, norm, model) for param_arr, param_name, norm in zip(param_arrays, param_names, norm_list)]
    
    if model == branch_precursor or model == branch_competetive:
        df_arr = [get_relative_readouts(df) for df in df_arr]
    
    df = pd.concat(df_arr)
    
    return df


def update_dict(d, val, name):
    
    d = dict(d)
    d[name] = val
    
    return d

def update_dicts(dicts, val, name):
    dicts = [update_dict(dic, val, name) for dic in dicts]
    return dicts

def get_tidy_readouts(p, param_name, time, cond, cond_names, model):
    cond = [update_dict(d, p, param_name) for d in cond]
    
    df = run_exp(time, cond, cond_names, model)
    df2 = generate_readouts(df, time)
    df2 = pd.melt(df2.reset_index(), 
                  id_vars = ["cond", "cell"],
                  value_vars = ["area", "peak", "tau", "decay"], 
                  var_name = "readout",
                  value_name = "y")
    
    df2["x"] = p

    return df2

def f(x,y):
    x["cond2"] = y
    return x

def multi_exp(time, cond_list, cond_names, cond_names2, 
              model = th_cell_diff, initial_cells = 1):
    """
    cond list should be a list of lists inner list for single exp with diff conditions,
    outer list list of second conditions
    """
    assert len(cond_list) == len(cond_names2)
    
    exp_list = [run_exp(time, cond, cond_names, model, initial_cells) for cond in cond_list]

    exp_list = [f(exp, name) for exp, name in zip(exp_list, cond_names2)]
    exp = pd.concat(exp_list)
    
    return exp

def convert_name(name: str) -> str:
    """
    take input name and return as string that is nicer to read for plotting
    """
    d = {
        "beta" : r"$\beta$",
        "alpha" : r"$\alpha$",
        "d_eff" : "death rate",
        "beta_p": r"$\beta_p$",
        "beta1" : r"$\beta_1$",
        "crit_timer" : "t0",
        "rate_il2" : "rate il2"}
    
    n = d[name]
    
    return n