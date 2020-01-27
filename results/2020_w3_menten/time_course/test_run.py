#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:34:11 2020

@author: burt
"""

import numpy as np
import seaborn as sns

from parameters import d
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/ode_models/")

from analysis_module import run_exp, multi_exp, generate_readouts
from ode_models import il2_model, timer_model, null_model, C_model
import pandas as pd
sns.set(context = "talk")
# =============================================================================
# 
# =============================================================================

# settings
time = np.arange(0,20,0.01)
model = timer_model
sim_name = "test"
model_name = "C"

# time course
df = run_exp(time, d, model, sim_name, model_name)
g = sns.relplot(data = df, x = "time", y = "cells", kind = "line")
