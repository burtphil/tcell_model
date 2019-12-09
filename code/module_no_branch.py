#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 09:48:52 2019

@author: burt
module: parameter scan for cell models without branching!
two models are supported:
- model1 with direct naive to effector transition, naive cells can secrete IL2
- model2 with naive --> precursor --> effector transition, precursors can secrete IL2
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import interpolate
import itertools

sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

def th_cell_diff(th_state, time, d):
    """
    model2
    takes state vector to differentiate effector cells as linear chain
    needs alpha and beta(r) of response time distribution, probability
    and number of precursor cells
    """

    dt_state = np.zeros_like(th_state)
    tnaive = np.sum(th_state[:-(2*d["alpha_p"])])
    teff = np.sum(th_state[-(2*d["alpha_p"]):-d["alpha_p"]])
    tnoil2 = np.sum(th_state[-d["alpha_p"]:])

    #carrying capacity
    x_tot = np.sum(th_state)
    conc_il2 = d["rate_il2"]*(tnaive+teff)

    # mm kinetic feedback implementation  
    fb_ifn = 0
    if d["fb_ifn"] != 0:
        conc_ifn = d["rate_ifn"]*(teff+tnoil2)
        fb_ifn = (d["fb_ifn"]*conc_ifn**3)/(conc_ifn**3+d["K_ifn"]**3)
  
    beta = (fb_ifn+1)*d["beta"]

    beta_p = d["beta_p"]
    
    # check if crit has been updated from F to T, if not
    # check if crit should be updated, then store t0 time of update
    # if crit has been updated, update prolif rate
    #if d["mode"] in ["il7", "il2", "timer"]:
    #    beta_p = beta_p*(1-(x_tot/d["crit_il7"]))
    
    if d["mode"] in ["il7", "il2", "timer"]:
        if d["crit"] == True:
            beta_p = beta_p*np.exp(-d["decay_p"]*(time-d["t0"]))
        else:
            if d["mode"] == "timer" and time > d["crit_timer"]:
                d["t0"] = time
                d["crit"] = True
            if d["mode"] == "il2" and conc_il2 < d["crit_il2"]:
                d["t0"] = time
                d["crit"] = True
            if d["mode"] == "il7" and x_tot > d["crit_il7"]:
                d["t0"] = time
                d["crit"] = True                

    for j in range(len(th_state)):
        #print(j)
        if j == 0:
            dt_state[j] = d["b"]-beta*th_state[j] 
            
        elif j < d["alpha"]:
            dt_state[j] = beta*th_state[j-1]-(beta+d["d_prec"])*th_state[j]
            
        elif j == d["alpha"]:
            dt_state[j] = beta*th_state[j-1] + (2*beta_p*th_state[(d["alpha"]+d["alpha_p"]-1)]) - (d["d_eff"]+d["beta_sad"]+beta_p)*th_state[j]       

        elif j < (d["alpha"]+d["alpha_p"]):
            dt_state[j] = beta_p*th_state[j-1]-(beta_p+d["d_eff"]+d["beta_sad"])*th_state[j]   
        
        elif j == (d["alpha"]+d["alpha_p"]):
            dt_state[j] = d["beta_sad"]*teff + (2*beta_p*th_state[-1]) - (d["d_eff"]+beta_p)*th_state[j]       
        
        else:
            dt_state[j] = beta_p*th_state[j-1]-(beta_p+d["d_eff"])*th_state[j]
            
    return dt_state

def th_cell_diff2(th_state, time, d):
    """
    model2
    takes state vector to differentiate effector cells as linear chain
    needs alpha and beta(r) of response time distribution, probability
    and number of precursor cells
    """
    d = dict(d)
    
    dt_state = np.zeros_like(th_state)
    tnaive = np.sum(th_state[:-d["alpha_p"]])
    teff = np.sum(th_state[-d["alpha_p"]:])
    #carrying capacity
    x_tot = np.sum(th_state)
    #beta_p = d["beta_p"]*(1-(x_tot/d["C"])) 
    beta_p = d["beta_p"]
    # timer
    #if time > d["t0"]:
    #    beta_p = beta_p * np.exp(-1*(time-d["t0"]))
    
    # IL2
    #if d["mode"] == "IL2":
    #    fb_IL2 = d["fb_IL2"]*tnaive/(tnaive+teff+1)
        #print(IL2)
    #    beta_p = fb_IL2*beta_p
    
    # fb on beta_0
    #IFN = d["rate_ifn"]*teff
    #beta = d["beta"]*(d["fb_ifn"]*IFN+1)/(IFN+1)
    beta = d["beta"]
    
    for j in range(len(th_state)):
        #print(j)
        if j == 0:
            dt_state[j] = d["b"]-beta*th_state[j] 
            
        elif j < d["alpha"]:
            dt_state[j] = beta*th_state[j-1]-(beta+d["d_prec"])*th_state[j]
            
        elif j == d["alpha"]:
            dt_state[j] = beta*th_state[j-1] + (2*beta_p*th_state[-1]) - (d["d_eff"]+beta_p)*th_state[j]       

        else:
            assert j > d["alpha"] and d["alpha_p"] > 1
            dt_state[j] = beta_p*th_state[j-1]-(beta_p+d["d_eff"])*th_state[j]   
            
            
    return dt_state
    
def run_model(time, d, model = th_cell_diff):
    d = dict(d)
    if model == th_cell_diff2:
        state0 = np.zeros(d["alpha"]+d["alpha_p"])
    else:    
        state0 = np.zeros(d["alpha"]+2*d["alpha_p"]) 
        
    state0[0] = 1.
    state = odeint(model, state0, time, args = (d,))
    
    return state
    
def get_cells(state, d, model):
    if model == th_cell_diff2:
        tnaive = state[:, :-d["alpha_p"]]
        teff = state[:, -d["alpha_p"]:]
        tnaive = np.sum(tnaive, axis = 1)
        teff = np.sum(teff, axis = 1)  
        all_cells = [tnaive, teff]
    
    else:
        tnaive = state[:, :-(2*d["alpha_p"])]
        teff = state[:, -(2*d["alpha_p"]):]
        
        tnaive = np.sum(tnaive, axis = 1)
        teff = np.sum(teff, axis = 1)

    all_cells = [tnaive, teff]
    return all_cells

def get_cells2(time, d, model = th_cell_diff):
    
    state = run_model(time, d, model)
    all_cells = get_cells(state, d, model)
    return all_cells

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

def get_decay_halftime(cells, time):
    """
    get the half-time of decay
    """
    
    cellmax = np.amax(cells)
    cellmin = cells[-1]
    # check if cells are maximum in size at the end of simulation
    crit1 = np.abs(cellmax-cellmin) > 1e-3
    crit2 = cellmax > cellmin
    
    # test that there is a global max before end of array
    if crit1 and crit2:
        peak_id = np.argmax(cells)
        cells = cells[peak_id:]
        time = time[peak_id:]
        
        # make sure there are at least two values in the array
        assert len(cells) > 1
        
        # interpolate to get time unter half of diff between max and arr end is reached
        celldiff = (cellmax - cellmin) /2
        celldiff = cellmax - celldiff
        f = interpolate.interp1d(cells, time)
        #print(cellmax, cellmin, celldiff)
        tau = f(celldiff)
    else:
        tau = np.nan
    
    return tau

def get_peak(cells, time):
    return np.amax(cells)

def get_area(cells, time):
    area = np.trapz(cells, time)
    return area

def get_delta_t(cells, time, thres = 0.75):
    thres = thres*np.amax(cells)
    #print(len(cells))
    #print(len(time))
    #print(time[cells>thres])
    t = time[(cells>thres)]
    
    # get indices of first and last entry in t from original time array
    # then linear interpolation with one element before and after first or last entry respectively
    i_front = np.nonzero(time == t[0])[0][0]
    i_back = np.nonzero(time == t[-1])[0][0]
    
    # check that first entry is not at beginning of original array
    if i_front != 0:
        t1 = (time[i_front-1]+time[i_front-1])/2
    else:
        t1 = 0
    
    # make sure last entry in t is not also last entry in time array
    if i_back < len(time)-1:
        t2 = (time[i_back]+time[i_back+1])/2
    else:
        t2 = t[-1]
    delta_t = t2-t1
    
    return delta_t

def get_readouts(cells, time):
    peak = get_peak(cells, time)
    area = get_area(cells, time)
    #dt = get_delta_t(cells, time)
    dt = 10
    peaktime = get_peaktime(cells, time)
    tau_decay = get_decay_halftime(cells, time)
    
    curve = peak/area
    readouts = [peak, area, curve, dt, peaktime, tau_decay]
    return readouts

def get_readout_subset(readouts, ids):
    readouts = [readouts[i] for i in ids]
    return readouts

def get_mu(d, name):
    devs = {
       "SD": ["alpha", "beta"],
       "SD_0" : ["alpha_0", "beta_0"],
       "SD_1" : ["alpha_1", "beta_1"],
       "SD_p" : ["alpha_p", "beta_p"],
       }
    x = devs[name]
    alpha = d[x[0]]
    beta = d[x[1]]
    mu = alpha/beta

    return mu

def vary_param(arr, name, ids, time, d, model = th_cell_diff):
    d = dict(d)
    #print(arr)
    reads_all_cells = []
    
    # get average transition times
    edge_names = ["SD", "SD_0", "SD_1", "SD_p", "chain"]
    if name in edge_names and name != "chain":
        mu = get_mu(d, name)
        
    for val in arr:
        if name not in edge_names:
            d[name] = val

        elif name == "SD":
            d["alpha"] = val
            d["beta"] = val/mu

        elif name == "SD_0":
            d["alpha_0"] = val
            d["beta_0"] = val/mu

        elif name == "SD_1":
            d["alpha_1"] = val
            d["beta_1"] = val/mu

        elif name == "SD_p":
            d["alpha_p"] = val
            d["beta_p"] = val/mu
        
        elif name == "chain":
            # note that this only works for mu = 1 for both prolif and diff 
            d["alpha_p"] = val
            d["beta_p"] = val
            d["alpha"] = val
            d["beta"] = val
            
        all_cells = get_cells2(time, d, model)
        teff = all_cells[-1]
        
        readouts = get_readouts(teff, time)
        readouts = get_readout_subset(readouts, ids)
        reads_all_cells.append(readouts)    

    return reads_all_cells

    
def param_scan(param_names, dicts, titles, ids, time, filename, p_labels, ylim = [-2,2], 
               model = th_cell_diff):
    """
    make parameter scan with one row per parameter and all models contained in each subplot
    only is good for few parameters
    readouts are in columns
    """
    assert len(titles) == len(dicts)
    
    dicts = [dict(dic) for dic in dicts]
    # get default readouts
    cells0 = [get_cells2(time, dic, model) for dic in dicts]
    teff0 = [cell[-1] for cell in cells0]
    reads0 = [get_readouts(teff, time) for teff in teff0]
    reads0 = np.asarray(reads0)
    reads0 = reads0[:, ids]
    
    edge_names = ["SD", "SD_0", "SD_1", "SD_p", "chain"]
    labels = ["peak height", "area", "peak/area", "peak width", r"$\tau_{1/2}$",
              r"$\lambda_{1/2}$"]
    
    colors = ["k", "tab:gray", "tab:olive", "tab:purple", "tab:red", "tab:cyan"]

    colors = [colors[i] for i in ids]
    labels = [labels[i] for i in ids]
    
    smallsize = (5, 4)
    bigsize = (4.5*len(dicts), 3.2*len(param_names))
    
    if len(param_names) < 2 and len(dicts) < 2:
        size = smallsize
    else:
        size = bigsize
    
    fig, axes = plt.subplots(len(param_names), len(dicts), squeeze = False,
                             figsize = size)

    for i in range(len(param_names)):
        param_name = param_names[i]
        

        ax = axes[i,:]
                
        xlabel = p_labels[i]
        
        arr = generate_array(param_name, dicts[0], edge_names)            
        reads = [vary_param(arr, param_name, ids, time, dic, model) for dic in dicts]
        #reads_fb1 = [np.log2(read) for read in reads_fb1]
        reads = np.asarray(reads)
         
        for j in range(len(dicts)):
            reads[j] = np.log2(reads[j]/reads0[j])
        
        for k in range(len(dicts)):

            ax[k].set_xlabel(xlabel)
            ax[k].set_ylabel("log2FC")
            
            for j in range(len(ids)):
                
                if param_name in ["b", "fb_ifn", "d_prec", "decay_p"]:
                
                        ax[k].plot(arr, reads[k,:,j], label = labels[j], 
                        c = colors[j])
                        ax[k].set_xticks([0,5,10])                        
                        ax[k].axvline(x=0, c = "tab:grey", ls = "--", lw = 2)
    
                elif param_name == "decay_p":
                    ax[k].plot(arr, reads[k,:,j], label = labels[j],
                    c = colors[j])
                    ax[k].set_xticks([0,.5,1])

                elif param_name in ["crit_timer", "crit_il7", "crit_il2", "rate_il2"]:
                    ax[k].plot(arr, reads[k,:,j], label = labels[j],
                    c = colors[j])
                    
                elif param_name not in edge_names:
                    ax[k].set_xscale("log")
                    ax[k].plot(arr, reads[k,:,j], label = labels[j],
                    c = colors[j])
                          
                    x0 = np.amin(arr)*10
                    xticks = [np.round(arr[0],2), np.round(x0,2), np.round(arr[-1],2)]
                    ax[k].axvline(x=x0, c = "tab:grey", ls = "--", lw = 2)
                    ax[k].axhline(y=0, c = "tab:grey", ls = "--", lw = 2)
                    ax[k].set_xlim([xticks[0], xticks[2]])
                    ax[k].set_xticks(xticks)                
                      
                else:
                    x = arr/arr**2
                    ax[k].set_xscale("log")
                    ax[k].scatter(x, reads[k,:,j],alpha = 0.7, label = labels[j],
                    c = colors[j])
                    ax[k].axvline(x=0.1, c = "tab:grey", ls = "--", lw = 2)
                    ax[k].axhline(y=0, c = "tab:grey", ls = "--", lw = 2)
                    ax[k].set_xticks(np.arange(0.1,1.1,0.1))
                    ax[k].set_xlim([0.07, 1.3])                                   
         
    #axes[0][0].legend()   
    
    for i in range(len(titles)):
        axes[0][i].set_title(titles[i])
    
    axes = axes.flatten()
    for ax in axes:
        if ylim != None:
            ax.set_ylim(ylim)
             
    plt.tight_layout()
    fig.savefig(filename+".pdf", bbox_inches = "tight")
    fig.savefig(filename+".svg", bbox_inches = "tight")
    plt.close()
    return fig

def generate_array(param_name, dic, edge_names):

    d = dict(dic)
    if param_name in ["b", "fb_ifn", "d_prec", "decay_p"]:
        pmin = 0
        pmax = 10              
        arr = np.linspace(pmin, pmax, 50)
    
    elif param_name == "beta_p":
        arr = np.logspace(1,2,50)
        
    elif param_name == "crit_timer":
        arr = np.linspace(0, 10, 50)

    elif param_name == "crit_il7":
        arr = np.linspace(0, 10, 50)

    elif param_name == "crit_il2":
        arr = np.linspace(0, 1, 50)

    elif param_name == "rate_il2":
        arr = np.linspace(0, 5, 50)
        
    elif param_name not in edge_names:
        val = d[param_name]
        
        vmin = np.log10(0.1*val)
        vmax = np.log10(10*val)
        arr = np.logspace(vmin, vmax, 50)
        
    else:
        arr = np.arange(1,13,1)

    return arr

def get_readouts2(time, params, model = th_cell_diff):
    all_cells = get_cells2(time, params, model)
    teff = all_cells[-1]
    readouts = get_readouts(teff, time)
    return readouts
    
def heatmap(arr1, arr2, name1, name2, time, params, readout_type, norm):
    z = []
    d = dict(params)    
    for a1, a2 in itertools.product(arr1, arr2):
        d[name1] = a1
        d[name2] = a2
        readouts = get_readouts2(time, d)
        readout = readouts[readout_type]
        readout = np.log2(readout/norm)
        z.append(readout)  
    
    #print(len(z))
    z = np.asarray(z)
    z = z.reshape(len(arr1), len(arr2))
    z = z[:-1, :-1]    
    z = z.T
    
    return z

def get_heatmap(arr1, arr2, name1, name2, time, params, readout_type, 
                vmin = None, vmax = None, title = None, label1 = None, label2 = None):
    
    norm = get_readouts2(time, params)[readout_type]
    z = heatmap(arr1, arr2, name1, name2, time, params, readout_type, norm)
    
    plot = plot_heatmap(arr1, arr2, name1, name2, z, vmin, vmax, title,
                        label1, label2)
    return plot

def plot_heatmap(arr1, arr2, name1, name2, val, vmin = None, vmax = None, 
                 title = None, label1 = None, label2 = None):
    
    fig, ax = plt.subplots(figsize = (6,4))
    color = "bwr"
    cmap = ax.pcolormesh(arr1, arr2, val, cmap = color, vmin = vmin, vmax = vmax)
    
    ax.set_xlabel(label1)
    ax.set_ylabel(label2)
    ax.set_title(title)
    cbar = plt.colorbar(cmap)
    cbar.set_label("lfc resp. size")
        
    plt.tight_layout()

    return fig

def sensitivity_analysis(arr, name, ids, time, dicts):
    """
    d: parameter dictionary
    param_names: list of names of parameters to be varied
    param_arrays: corresponding arrays that should be varied for each parameter
    """

    # get readouts for every model
    reads_null = vary_param(arr, name, ids, time, dicts[0])
    reads_il2 = vary_param(arr, name, ids, time, dicts[1])
    reads_timer = vary_param(arr, name, ids, time, dicts[2])
    
    reads = [reads_null, reads_il2, reads_timer]
    reads = np.asarray(reads)
    
    # divide every value by the default value
    for i in range(len(dicts)):
        reads[i] = reads[i,:,:]/reads[i,1,:]

    reads = np.log2(reads)   
    # divide by default values
    #readouts = np.log2(readouts / readouts[:,1,:])
    return reads


def plot_sensitivity_analysis(reads, labels, title):
    x = np.arange(len(labels))
    width = 0.2
    
    #print(readouts.shape())
    
    fig, ax = plt.subplots(figsize = (8,6))
    ax.bar(x, reads[0,0,:], label = "Null model", width = width, color = "tab:blue")   
    ax.bar(x-width, reads[1,0,:], label = "IL2", width = width, color = "tab:orange")   
    ax.bar(x+width, reads[2,0,:], label = "timer", width = width, color = "tab:green")
    
    # check if I have pos and neg coeff, then plot neg coeff as well
    c = reads[0,:,0]

    if len(c) == 3:
        ax.bar(x, reads[0,2,:], width = width, color = "tab:blue", hatch = "/")
        ax.bar(x-width, reads[1,2,:], width = width, color = "tab:orange", hatch = "/")
        ax.bar(x+width, reads[2,2,:], width = width, color = "tab:green", hatch = "/")
    
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("log2FC")
    ax.set_title(title)
    
    plt.xticks(rotation = "vertical")
    fig.tight_layout()
    #fig.savefig("../figures/sensitivity/"+"sensitivity_"+title+".pdf")
    #fig.savefig("../figures/sensitivity/"+"sensitivity_"+title+".svg")


def run_sensitivity(arr, name, ids, time, dicts, labels, title):
    readouts = sensitivity_analysis(arr, name, ids, time, dicts)
    assert len(ids) == len(labels), "check your labels"
    plot_sensitivity_analysis(readouts, labels, title)
    

def plot_timecourse(time, dicts, labels, filename = "timecourse"):
    """
    plot time course for given number of models (provided by dicts)
    """
    all_cells = [get_cells2(time, dic) for dic in dicts]
    teffs = [cells[-1] for cells in all_cells]
    
    fig, ax = plt.subplots()
    
    for teff, label in zip(teffs, labels):
        ax.plot(time, teff, label = label)
    
    ax.set_xlabel("time")
    ax.set_ylabel("cell dens. norm.")
    ax.legend()
        
    plt.tight_layout()
  


