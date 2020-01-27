#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 12:57:04 2020

@author: burt
"""
import numpy as np

def core_model_gamma(state, d, gamma):
       
    r_diff = d["r_diff"]
    n_div = d["n_div"]
    #r_death = d["r_death"]
    r_myc = d["r_myc"]
    
    dx = -r_diff*state[0]
    dy = n_div*r_diff*state[0] - r_diff*state[1]
    dz = n_div*r_diff*state[1] + gamma*state[2]
    
    d_myc = -r_myc*state[3]
    dt_state = [dx, dy, dz, d_myc]
    
    assert(len(state)==len(dt_state))
    return dt_state

def core_model_timer_il2(state, d, gamma):
    r_diff = d["r_diff"]
    n_div = d["n_div"]
    #r_death = d["r_death"]
    r_myc = d["r_myc"]
    
    dx = -r_diff*state[0]
    dy = n_div*r_diff*state[0] - r_diff*state[1]
    dz = n_div*r_diff*state[1] + gamma*state[2]

    y = state[1]
    z = state[2]    
    c_il2 = (d["r_il2"]*y) / (y+z+d["K_il2"])

    d_myc = -(r_myc/(c_il2+d["K_il2_myc"])) * state[3]
    dt_state = [dx, dy, dz, d_myc]
    
    assert(len(state)==len(dt_state))
    return dt_state    

def core_model_menten(state, d, r_d):
       
    r_diff = d["r_diff"]
    n_div = d["n_div"]
    #r_death = d["r_death"]
    r_myc = d["r_myc"]
    r_p = d["r_p"]
    
    dx = -r_diff*state[0]
    dy = n_div*r_diff*state[0] - r_diff*state[1]
    dz = n_div*r_diff*state[1] + (r_p-r_d)*state[2]
    
    d_myc = -r_myc*state[3]
    dt_state = [dx, dy, dz, d_myc]
    
    assert(len(state)==len(dt_state))
    return dt_state

def ode_model(state, time, d, homeostasis_model, core_model):
    """
    calculate gamma depending on model type
    """
    gamma = homeostasis_model(state, time, d)
    dt_state = core_model(state, d, gamma)

    return dt_state     

    
def null_model(state, time, d):
    gamma = d["gamma"]
        
    return gamma 


def null_model_menten(state, time, d):
    r_d = d["r_d"]
        
    return r_d 

def il2_model_menten(state, time, d):
    y = state[1]
    z = state[2]
    #assert y >= 0
    #assert z >= 0
    
    c_il2 = (d["r_il2"]*y) / (y+z+d["K_il2"])  
    r_d = d["r_d"] * d["K_il2"] / (c_il2+d["K_il2"])
    
    return r_d
            
def C_model_menten(state, time, d):
    z = state[2]
    #assert z >= 0
    
    c_C = (d["r_C"]*1.) / (z+d["K_C"]) 
    
    r_d = d["r_d"] * d["K_C"] / (c_C+d["K_C"])
    #if time > 10 : print(r_d)
    return r_d

def timer_model_menten(state, time, d):
    myc = state[3]
    #myc = myc if myc >= 0 else 0
    
    r_d = d["r_d"] * d["K_myc"] / (myc+d["K_myc"])

    return r_d    

#deprecated
def il2_model(state, time, d):

    gamma = d["gamma"]

    if d["crit"] == False:
        gamma = d["gamma"]-1*(time-d["t0"])
    
    else:
        y = state[1]
        z = state[2]
        conc_il2 = d["r_il2"]*y/(y+z+d["K_il2"])

        if conc_il2 < d["c_il2"] and z > 0.01:
            d["crit"] = False
            d["t0"] = time
        
    return gamma

def C_model(state, time, d):
    
    gamma = d["gamma"]
  
    if d["crit"] == False:
        gamma = d["gamma"]-1*(time-d["t0"])
    
    else:
        z = state[2]
        conc_c = d["r_C"]*1/(z+d["K_C"])
        if conc_c < d["c_C"] and z >  0.01:
            d["crit"] = False
            d["t0"] = time
        
    return gamma

def timer_model(state, time, d):

    gamma = d["gamma"]
  
    if d["crit"] == False:
        gamma = d["gamma"]-1*(time-d["t0"])
    
    else:
        myc = state[3]
        
        if myc < d["c_myc"]:
            d["crit"] = False
            d["t0"] = time
            
    return gamma   

def timer_il2_model(state, time, d):
    gamma = d["gamma"]
    
    if d["crit"] == False:
        gamma = d["gamma"]-1*(time-d["t0"])
        
    else:
        myc = state[3]
        y = state[1]
        z = state[2]
        conc_il2 = d["r_il2"]*y/(y+z+d["K_il2"])
        
        crit1 = myc < d["c_myc"]
        crit2 = conc_il2 < d["c_il2"] and z > 0.01
        
        if (crit1 or crit2):
            d["crit"] = False
            d["t0"] = time
    
    return gamma