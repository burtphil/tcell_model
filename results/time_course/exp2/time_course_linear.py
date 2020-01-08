#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 10:24:11 2020

@author: burt
"""

import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
from tcell_parameters import d_il7, d_il2, d_timer, d_null, d_il2_timer
from test_module import run_exp, generate_readouts
import module_models as models
import numpy as np
import seaborn as sns

sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

colors = ["#3498db", "#95a5a6", "#e74c3c", "#34495e"]     
sns.set_palette(colors)     

# =============================================================================
# define exp conditions
# =============================================================================
cond = [d_null]
cond_names = ["Null"]
d_il7["death_mode"] = True
time = np.arange(0, 50, 0.01)
model = models.th_cell_diff

df = run_exp(time, cond, cond_names, model = model)

sns.relplot(x = "time", y = "value", kind = "line", data = df)

test = df.value

readouts = generate_readouts(df, time)
print(readouts)