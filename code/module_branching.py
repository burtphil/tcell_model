#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 16:18:30 2019

@author: burt

implement different branching scenarios

"""
import numpy as np
from scipy.integrate import odeint
from scipy import interpolate

def diff_effector(state, th0, alpha, beta, beta_p, p, d):
    """
    takes state vector to differentiate effector cells as linear chain
    needs alpha and beta(r) of response time distribution, probability
    and number of precursor cells
    """
    dt_state = np.zeros_like(state)
    #print(len(state))
    if alpha == 1:
        for j in range(len(state)):
            if j == 0:
                dt_state[j] = p*beta*th0+2*beta_p*state[-1]-(beta_p+d["d_eff"])*state[j]
            else:
                dt_state[j] = beta_p*state[j-1]- (beta_p+d["d_eff"])*state[j] 
    
    else:                    
        for j in range(len(state)):
            if j == 0:
                dt_state[j] = p*beta*th0 - (beta+d["d_prec"])*state[j]                
            elif j < (alpha-1):
                dt_state[j] = beta*state[j-1]-(beta+d["d_prec"])*state[j]                
            elif j == (alpha-1):
                # the problem with the 4 and 2 is that since differentiation takes 1 day it should divide twice giving 4 cells
                # however, if it has arrived in the final states if should double every half day
                dt_state[j] = beta*state[j-1]+2*beta_p*state[-1] - (d["d_eff"]+beta_p)*state[j]  

            else:
                assert j >= alpha
                dt_state[j] = beta_p*state[j-1]- (beta_p+d["d_eff"])*state[j] 
               
    return dt_state

def diff_precursor(state, th0, alpha, beta, beta_p, p_norm, d):
    """
    takes state vector to differentiate effector cells as linear chain
    needs alpha and beta(r) of response time distribution, probability
    and number of precursor cells
    """
    dt_state = np.zeros_like(state)

    for j in range(len(state)):
        if j == 0:
            dt_state[j] = p_norm*d["beta0"]*th0 - (beta+d["d_prec"])*state[j]             
        elif j < alpha:
            dt_state[j] = beta*state[j-1]- (beta+d["d_prec"])*state[j]                
        elif j == alpha:
            # the problem with the 4 and 2 is that since differentiation takes 1 day it should divide twice giving 4 cells
            # however, if it has arrived in the final states if should double every half day
            dt_state[j] = beta*state[j-1]+2*beta_p*state[-1] - (d["d_eff"]+beta_p)*state[j] 
            
        else:
            assert j > alpha        
            dt_state[j] = beta_p*state[j-1]-(beta_p+d["d_eff"])*state[j] 
                 
    return dt_state

def branch_competetive(state, time, d):
    """
    takes state vector to differentiate effector cells as linear chain
    needs alpha and beta(r) of response time distribution, probability
    and number of precursor cells
    """

    th0 = state[0]    
    th1 = state[1:(d["alpha1"]+d["alpha1_p"])]
    th2 = state[(d["alpha1"]+d["alpha1_p"]):]
    
    #print(len(state), len(th1))
    ### get all cytokine secreting cells    
    th1_all = np.sum(th1[-d["alpha1_p"]:])
    th2_all = np.sum(th2[-d["alpha2_p"]:])
    ### calculate cytokine concentrations
    cyto_1 = d["beta_cyto_1"]*th1_all + d["ifn_ext"]
    cyto_2 = d["beta_cyto_2"]*th2_all + d["il21_ext"]
    ### calculate cytokine effect on rate
    fb1 = d["fb_1"]*cyto_1**3/(cyto_1**3+d["K_1"]**3)
    fb2 = d["fb_2"]*cyto_2**3/(cyto_2**3+d["K_2"]**3)
    ### update differantiation rate
    beta1 = d["beta1"]*(1+fb1)
    beta2 = d["beta2"]*(1+fb2)
    
    ### differentiate effectors th1    
    alpha = d["alpha1"]
    p = 1.
    dt_th1 = diff_effector(th1, th0, alpha, beta1, d["beta1_p"], p, d)
    ### differentiate effectors th2
    alpha = d["alpha2"]
    p = 1.
    dt_th2 = diff_effector(th2, th0, alpha, beta2, d["beta2_p"], p, d)
    
    ### combine all cells
    dt_th0 = -(beta1+beta2)*th0
    dt_state = np.concatenate(([dt_th0], dt_th1, dt_th2))

    return dt_state

def branch_precursor(state, time, d):
    """
    takes state vector to differentiate effector cells as linear chain
    needs alpha and beta(r) of response time distribution, probability
    and number of precursor cells
    """

    th0 = state[0]
    
    dt_th0 = -d["beta0"]*th0
    
    
    th1 = state[1:(d["alpha1"]+d["alpha1_p"]+1)]
    th2 = state[(d["alpha1"]+d["alpha1_p"]+1):]
    #print(len(state), len(th1))
    ### get all cytokine secreting cells    
    th1_all = np.sum(th1[-d["alpha1_p"]:])
    th2_all = np.sum(th2[-d["alpha2_p"]:])
    ### calculate cytokine concentrations
    cyto_1 = d["beta_cyto_1"]*th1_all + d["ifn_ext"]
    cyto_2 = d["beta_cyto_2"]*th2_all + d["il21_ext"]
    
    ### calculate probability 
    p1 = d["fb_1"]*cyto_1**3/(cyto_1**3+d["K_1"]**3)
    p2 = d["fb_2"]*cyto_2**3/(cyto_2**3+d["K_2"]**3)
    # account for default probability and feedback strength
    p1 = (p1+1)*d["p1_def"]
    p2 = (p2+1)*(1-d["p1_def"])

    ### update differentiation rate
    p1_norm = p1/(p1+p2)
    p2_norm = 1-p1_norm

    alpha1 = d["alpha1"]
    beta1 = d["beta1"]
    dt_th1 = diff_precursor(th1, th0, alpha1, beta1, d["beta1_p"], p1_norm, d)

    alpha2 = d["alpha2"]
    beta2 = d["beta2"]

    dt_th2 = diff_precursor(th2, th0, alpha2, beta2, d["beta2_p"], p2_norm, d)
    
    dt_state = np.concatenate(([dt_th0], dt_th1, dt_th2))

    return dt_state

def run_model(d, t, fun, output = "cells"):
    if fun == branch_competetive:
        
        y= np.zeros(d["alpha1"]+d["alpha1_p"]+d["alpha2"]+d["alpha2_p"]-1)
        y[0] = 1
        
    if fun == branch_precursor:
        
        y= np.zeros(d["alpha1"]+d["alpha1_p"]+d["alpha2"]+d["alpha2_p"]+1)
        y[0] = 1    
    

    state = odeint(fun, y, t, args = (d,))
        
    if output == "cells":
        cells = get_cells(state, d, fun)
        return cells
    
    elif output == "readouts":
        cells = get_cells(state, d, fun)
        readouts = get_readouts(cells, t)
        return readouts
    else:
        return state

def get_cells(state, d, fun):

    if fun == branch_competetive:
        th1 = state[:, 1:(d["alpha1"]+d["alpha1_p"])]
        th2 = state[:, (d["alpha1"]+d["alpha1_p"]):]
        
    if fun == branch_precursor:
        th1 = state[:, 1:(d["alpha1"]+d["alpha1_p"]+1)]
        th2 = state[:, (d["alpha1"]+d["alpha1_p"]+1):]    
    
    th1 = th1[:,-d["alpha1_p"]:]
    th2 = th2[:,-d["alpha2_p"]:]
    
    th1 = np.sum(th1, axis = 1)
    th2 = np.sum(th2, axis = 1)
    
    cells = [th1, th2]
    
    return cells

def get_cell_rel(th1, th2):
    total_cells = th1 + th2
    total_cells[total_cells == 0] = np.nan
    th1_r = th1/total_cells
    th2_r = th2/total_cells
    
    cells_r = [th1_r, th2_r]
    
    return cells_r

def get_cell_area(th1, th2, time):
    area = np.trapz(th1, time) + np.trapz(th2, time)
    return area

def get_readouts(cells, time):
    th1 = cells[0]
    th2 = cells[1]
    
    # area
    th1_area = np.trapz(th1, time)
    th2_area = np.trapz(th2, time)
    area = th1_area+th2_area
    rel_area = th1_area/th2_area if th2_area != 0 else np.nan
    # peak
    th1_peak = np.amax(th1)
    th2_peak = np.amax(th2)
    rel_peak = th1_peak/th2_peak if th2_peak != 0 else np.nan
    
    #tau
    th1_tau = get_peaktime(th1, time)
    th2_tau = get_peaktime(th2, time)
    rel_tau = th1_tau/th2_tau if th2_tau != 0 else np.nan
    
    readouts = [rel_area, rel_peak, rel_tau, area, th1_peak, th2_peak]
    
    return readouts

def get_peaktime(cells, time):
    # only look at array where there are no nans

    cells = np.asarray(cells)    
    # test if cells are not constant
    if np.std(cells) < 0.001:
        tau = np.nan
    # get peak from global maximum
    else:
    # get max idx
        peak_idx = np.argmax(cells)
    # get max value
        peak = cells[peak_idx]
        #print(peak)
        cells = cells[:(peak_idx+1)]
        time = time[:(peak_idx+1)]
        
        # assert that peak is not at beginning
        if peak_idx <= 1:
            tau = np.nan
        else:
            f = interpolate.interp1d(cells, time)
            peak_half = peak / 2.
            tau = f(peak_half)
            tau = float(tau)
            
    return tau   

def update_dict(val, name, d):
    edge_names = ["SD", "SD_0", "SD_1", "SD_p", "chain"]
    
    d = dict(d)
    if name not in edge_names:
        d[name] = val
    
    return d

def vary_param(arr, name, time, d, fun):

    dict_list = [update_dict(val, name, d) for val in arr]
    readouts = [run_model(dic, time, fun, output = "readouts") for dic in dict_list]

    return readouts