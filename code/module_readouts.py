# -*- coding: utf-8 -*-
import numpy as np
from scipy import interpolate
import pandas as pd

def get_readouts(df):
    
    d = {}
    time = df.time
    cells = df.value
    
    d['peak'] = get_peak(time, cells)
    d['tau'] = get_peaktime(time, cells)
    d['area'] = get_area(time, cells)
    d['decay'] = get_decay(time, cells)
    series = pd.Series(d, index=['peak', 'tau', 'area', 'decay'])
    
    return series

def get_peaktime(time, cells):

    cells = cells.array

    # only look at array where there are no nans    
    cellmax = np.amax(cells)
    cellmin = cells[-1]
    # check if cells are maximum in size at the end of simulation
    crit1 = np.abs(cellmax-cellmin) > 1e-3
    crit2 = cellmax > cellmin
    crit3 = np.std(cells) > 0.001
    # test if cells are not constant
    if crit1 and crit2 and crit3:

        peak_idx = np.argmax(cells)
    # get max value
        peak = cells[peak_idx]
        peak_half = peak / 2.
        #print(peak)
        cells = cells[:(peak_idx+1)]
        time = time[:(peak_idx+1)]
        # assert that peak is not at beginning
        if peak_idx <= 3:
            tau = np.nan
        # assert that peak half is in new cell array
        elif np.all(peak_half<cells):
            tau = np.nan
        else:
            f = interpolate.interp1d(cells, time)           
            tau = f(peak_half)
            tau = float(tau)
            
    else:
        tau = np.nan
            
    return tau    

def get_decay(time, cells):
    """
    get the half-time of decay
    """
    
    cells = cells.array
    cellmax = np.amax(cells)
    cellmin = cells[-1]
    # check if cells are maximum in size at the end of simulation
    crit1 = np.abs(cellmax-cellmin) > 1e-3
    crit2 = cellmax > cellmin
    crit3 = np.std(cells) > 0.001
    
    # test that there is a global max before end of array
    if crit1 and crit2 and crit3:
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

def get_area(time, cells):
    
    cells = cells.array 
    cellmax = np.amax(cells)
    cellmin = cells[-1]
    # check if cells are maximum in size at the end of simulation
    crit1 = np.abs(cellmax-cellmin) > 1e-3
    crit2 = cellmax > cellmin
    crit3 = np.std(cells) > 0.001
    
    if crit1 and crit2 and crit3:
        area = np.trapz(cells, time)
    else: 
        area = np.nan
        
    return area

def get_peak(time, cells):
    
    cells = cells.array 
    cellmax = np.amax(cells)
    cellmin = cells[-1]
    # check if cells are maximum in size at the end of simulation
    crit1 = np.abs(cellmax-cellmin) > 1e-3
    crit2 = cellmax > cellmin
    crit3 = np.std(cells) > 0.001
    
    if crit1 and crit2 and crit3:
        peak = cellmax
    else: 
        peak = np.nan
        
    return peak