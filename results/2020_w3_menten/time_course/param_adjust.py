#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 12:13:37 2020

@author: burt
script to find normalization conditions
"""

import numpy as np
import seaborn as sns

from parameters_menten import d
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/ode_models/")

from analysis_module import multi_exp, generate_readouts
from ode_models import null_model_menten, il2_model_menten, timer_model_menten, C_model_menten, core_model_menten
import pandas as pd
sns.set(context = "talk")


time = np.arange(0,20,0.01)
model_list = [timer_model_menten, il2_model_menten, C_model_menten, null_model_menten]
sim_names = ["test", "a", "b"]
model_names = ["Timer", "IL2", "K", "Null"]
d_list = [d]

# time course
df = multi_exp(time, d_list, model_list, sim_names, model_names, core_model = core_model_menten)
g = sns.relplot(data = df, x = "time", y = "cells", hue = "model", kind = "line")

readouts = generate_readouts(df)
print(readouts)