#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:09:30 2019

@author: burt
module for plotting functions
"""

def lineplot(ax, x, y_arrays, col = None, styles = None):
    """
    provide multiple y arrays together with colors or linestyles
    plot on some axes object
    y arrays should be a list where each entry has equal shape to x
    returns ax object
    """
    if col == None and styles == None:
        for y in y_arrays:
            ax.plot(x, y)
    elif col != None:
        assert len(col) == len(y_arrays)
        for y, c in zip(y_arrays, col):
            ax.plot(x, y, c)        
    else:
        assert len(styles) == len(y_arrays)
        for y, ls in zip(y_arrays, styles):
            ax.plot(x, y, ls = ls)
            
    return ax

