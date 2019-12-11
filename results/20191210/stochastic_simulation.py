# -*- coding: utf-8 -*-
"""
Spyder Editor

run a stochastic simulation for a lienar diff process
explore effect of distributions of proliferation and differentiation
"""
import numpy as np

nstates = 5

def make_cells(ncells):
    cells = np.zeros((ncells, 2))
    cells[0,0] = 1
    cells[0,1] = np.random.rand()
    # ncells are total number of cells allocated to memory incl dead cells

    # assign prolif times
    # assign diff times
    return cells

def exp_cdf(t, size, rate = 1):
    arr = np.ones(size)
    val = 1-np.exp(-t*rate)
    arr = arr*val
    return arr

def count_cells(cell_vector):
    # check how many cells of each type are there
    return cells

def add_cell(cell_vector):
    # assign new prolif time
    # assign new death time
    return 0

cells = make_cells(100)

def run_sim(time, cells):
    
    counter = np.zeros_like(time)
    
    for i in range(len(time)-1):
        t = time[i]
        t_new = time[i+1]
        
        #get indices of alive cells
        n_alive = get_alive_cells()
        exp_arr = exp_cdf(t, n_alive)
        
        # get random numbers from alive cells to compare with exp array    
        prob_div = alive_cells[:, 1]
        update_arr = exp_arr > prob_div
        
        #update arr contains information which cells divide
        # get number of cell divisions
        n_new_cells = np.sum(update_arr)
        # for divided cells, assign new division times
        
        # create number of new cells from dead cells
        counter[i] = count_cells()
    return counter    
time = np.arange(0,10,0.1)
