#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 10:24:33 2020

@author: burt
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

y0 = [1,0,0]

diff = 1
prolif1 = 10
prolif2 = 1
death = 1


def ode(state, time, diff, prolif1, prolif2, death):
    dt_state = np.zeros_like(state)

    dt_state[0] = -diff*state[0]
    dt_state[1] = prolif1*state[0]-diff*state[1]
    dt_state[2] = prolif1*state[1]+(prolif2-death)*state[2]
    
    return dt_state

time = np.arange(0,5,0.1)
state = odeint(ode, y0, time, args = (diff, prolif1, prolif2, death))

fig, ax = plt.subplots()
ax.plot(time, state)