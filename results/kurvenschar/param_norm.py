#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 15:13:46 2020

@author: burt
"""

import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
from tcell_parameters import d_il2, d_timer, d_il2_timer
import matplotlib.pyplot as plt
from test_module import norm_readout
import module_models as models
import numpy as np
import seaborn as sns
import pandas as pd

sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})
colors = ["#3498db", "#95a5a6", "#e74c3c", "#34495e"]     
sns.set_palette(colors)

time = np.arange(0,50,0.01)
pname = "crit_timer"
pmin = 0
pmax = 5
cond = dict(d_timer)
cond_names = ["timer"]

out = norm_readout(pname, pmin, pmax, time, cond, cond_names)

