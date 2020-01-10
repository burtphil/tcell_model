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
from module_readouts import get_area, get_decay, get_peaktime, get_peak
from functools import partial
import warnings
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
            keep_naive = False, counter = 3, adjust_time = False):
    
    if model == th_cell_diff:
        cell_names = ["naive", "teff"]
    else:
        cell_names = ["Th1", "Tfh"]
        
    states = [run_model(time, d, model, initial_cells) for d in cond]
    df = df_from_exp(time, states, cond_names, cell_names)
    
    if keep_naive == False:
        df = df[df["cell"] != "naive"]
    
    if adjust_time == True:    
        # check that time course returns to zero for diff. cells
        # if not, increase time span and run again
        last_el = df.value.iloc[-1]
        
        if last_el < 0.01:
            return df
        elif counter > 3:
            warnings.warn("no appropriate time frame found where cells return to zero")
            return df   
        else:
            counter = counter + 1
            time = np.arange(0, 5*time[-1],0.01)
            print("adjusting time frame for simulation...")
            return run_exp(time, cond, cond_names, model, initial_cells, keep_naive, counter,
                           adjust_time)

    else:
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
    get_peak_partial = partial(get_peak, time)
    
    df = df.groupby(groups).agg(
            peak = ("value", get_peak_partial),
            tau = ("value", get_peaktime_partial),
            area = ("value", get_area_partial),
            decay = ("value", get_decay_partial),            
            )
    
    return df

def vary_param(param_arr, param_name, time, cond, cond_names, norm, model = th_cell_diff,
               convert = False, adjust_time = True):
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
    
    df_arr = [get_tidy_readouts(p, param_name, time, cond, cond_names, model = model,
                                adjust_time = adjust_time) for p in param_arr]
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
    # add xnorm column to normalise x axis for param scans
    df["xnorm"] = df["x"] / norm
    
    if convert == True:
        pname = convert_name(param_name)
    else: 
        pname = param_name
        
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

def multi_param(param_arrays, param_names, time, cond, cond_names, norm_list, 
                model = th_cell_diff, relative_readouts = False, convert = False,
                adjust_time = True):
    
    assert len(norm_list) == len(param_arrays) == len(param_names)
    df_arr =[vary_param(param_arr, 
                        param_name, 
                        time, 
                        cond, 
                        cond_names, 
                        norm, 
                        model, 
                        convert,
                        adjust_time) for param_arr, param_name, norm in zip(param_arrays, param_names, norm_list)]
    
    if (model == branch_precursor or model == branch_competetive) and relative_readouts == True:
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

def norm_to_readout1(param_arr, param_name, time, cond, cond_name, norm, norm_cond, model = th_cell_diff,
                    convert = True):
    # dummy list because vary param expects a list as input
    cond = [cond]
    cond_names = [cond_name]
    
    df = vary_param(param_arr, param_name, time, cond, cond_names, norm, model,
                    convert)   

    df = df.loc[df["readout"] == "area", ["x", "y"]]  
    
    df["crit"] = np.abs(norm_cond-df["y"])
    df = df.reset_index()

    min_idx = df["crit"].idxmin()

    assert 0 < min_idx < df["crit"].size-1, "normalization not possible in given range"
    crit_min = df["crit"].min()
 
    out = float(df.loc[min_idx, ["x"]])

    return out, crit_min
    
def norm_to_readout(param_name, time, cond, cond_name, norm_cond = 1,
                    model = th_cell_diff, convert = True, counter = 0, param_guess = 1.,
                    sample_rate = 0.1):
    """
    provide arr and parameter name
    returns value of parameter name that
    """

    assert(param_guess > 0)
    param_arr = np.linspace(0.1*param_guess, 10*param_guess, 20)
    
    if counter > 0:
        sample_min = param_guess - 1.1*sample_rate
        sample_max = param_guess + 1.1*sample_rate
        
        sample_min = sample_min if sample_min > 0 else 0
        param_arr = np.linspace(sample_min, sample_max, 20)
    
   
    norm = (param_arr[-1]+param_arr[0]) / 2

    out, crit_min = norm_to_readout1(param_arr, param_name, time, cond, cond_name, norm,
                                 norm_cond, model, convert)
    
    if counter > 5:
        return out
    
    elif crit_min > 0.01:
        
        # counter is there to stop recursive function call if things go bad
        counter = counter + 1
        norm = out
        param_guess = out
        sample_rate = np.diff(param_arr)[0]
        return norm_to_readout(param_name, time, cond, cond_name, norm_cond,
                               model, convert, counter, param_guess, sample_rate)
    else:
        return out


def norm_readout(pname, pmin, pmax, time, cond, cond_names, model = th_cell_diff,
                 norm_cond = 2):
    
    cond = dict(cond)
    cond1 = dict(cond)
    cond2 = dict(cond)
    cond1[pname] = pmin
    cond2[pname] = pmax
    cond1 = [cond1]
    cond2 = [cond2]
    cond = [cond]

    df1 = run_exp(time, cond1, cond_names, model)
    df2 = run_exp(time, cond2, cond_names, model)
    read1 = generate_readouts(df1, time)
    read2 = generate_readouts(df2, time)
    area1 = read1.area[0]
    area2 = read2.area[0]
    
    assert area1 < norm_cond < area2
    guess = (pmin+pmax) / 2
    crit = False
    
    while crit == False:

        cond[0][pname] = guess

        df = run_exp(time, cond, cond_names, model)
        read = generate_readouts(df, time)
        area = read.area[0]
       
        if area < norm_cond:
            pmin = guess
            guess = (guess+pmax)/2
            
        else:
            pmax = guess
            guess = (guess+pmin)/2
                      
        if np.abs(area-norm_cond) < 0.1:
            crit = True
            
      
    return guess


def get_tidy_readouts(p, param_name, time, cond, cond_names, model, adjust_time):

    cond = [update_dict(d, p, param_name) for d in cond]
    
    df = run_exp(time, cond, cond_names, model, adjust_time = adjust_time)
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
        "rate_il2" : "rate il2",
        "p1_def" : "$p_1$"}
    
    n = d[name]
    
    return n

def array_from_dict(d, pname, res = 25, log = True):
    """
    takes dict, returns array for param scan
    with default values of dict * or / 100%
    """
    x = d[pname]
    xmin = x/10.
    xmax = x*10.
    arr1 = np.linspace(xmin, x, res, endpoint = False)
    arr2 = np.linspace(x, xmax)
    arr = np.concatenate((arr1, arr2))
    return arr