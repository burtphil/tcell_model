# -*- coding: utf-8 -*-

"""
understand rate effect in IL2 model for this plot time course and vary Il2 params
"""

import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
from tcell_parameters import d_prec, d_prec_il7, d_prec_il2, d_prec_timer
import matplotlib.pyplot as plt
from test_module import run_exp
import module_models as models
import numpy as np
import seaborn as sns
sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

colors = ["#3498db", "#95a5a6", "#e74c3c", "#34495e"]     
sns.set_palette(colors)     

# =============================================================================
# define exp conditions
# =============================================================================
cond = [d_prec_il2]
cond_names = ["IL2"]
time = np.arange(0, 6, 0.01)
model = models.branch_precursor


df = run_exp(time, cond, cond_names, model = model)

sns.relplot(x = "time", y = "value", hue = "cell", kind = "line", data = df)