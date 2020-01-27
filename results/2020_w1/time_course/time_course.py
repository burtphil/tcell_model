#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 10:24:11 2020

@author: burt
"""
from tcell_parameters import d_il7, d_il2, d_timer, d_null, d_il2_timer

import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
from test_module import run_exp, generate_readouts
import numpy as np
import seaborn as sns

sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})
# =============================================================================
# define exp conditions
# =============================================================================
cond = [d_timer]
cond_names = ["timer"]
#d_il7["death_mode"] = True
time = np.arange(0, 10, 0.001)

df = run_exp(time, cond, cond_names, adjust_time = True)

g = sns.relplot(x = "time", y = "value", kind = "line", data = df)
xlim = (0, None)
g.set(xlim = xlim)
test = df.value

readouts = generate_readouts(df)
print(readouts.area)