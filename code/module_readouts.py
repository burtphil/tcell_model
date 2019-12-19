# -*- coding: utf-8 -*-
import numpy as np
from scipy import interpolate

def get_peaktime(time, series):
    # only look at array where there are no nans
    cells = series.array   
    # test if cells are not constant
    if np.std(cells) < 0.001:
        tau = np.nan
    # get peak from global maximum
    else:
    # get max idx
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
        elif np.all(peak_half<cells):
            tau = np.nan
        else:
            f = interpolate.interp1d(cells, time)           
            tau = f(peak_half)
            tau = float(tau)
            
    return tau    

def get_decay(time, series):
    """
    get the half-time of decay
    """
    cells = series.array  
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

def get_area(time, series):
    cells = series.array
    area = np.trapz(cells, time)
    return area