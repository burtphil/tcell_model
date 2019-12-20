# -*- coding: utf-8 -*-
from tcell_parameters import d_null, d_il7, d_il2, d_timer
from tcell_parameters import d_il7_death, d_il2_death, d_timer_death
import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
import matplotlib.pyplot as plt
from test_module import multi_exp
import module_models as models
import numpy as np
import seaborn as sns

sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

# =============================================================================
# define exp conditions
# =============================================================================
cond1 = [d_null, d_il7, d_il2, d_timer]
cond2 = [d_null, d_il7_death, d_il2_death, d_timer_death]

cond_list = [cond1, cond2]
cond_names2 = ["prolif", "death"]
cond_names = ["Null", "K", "IL2", "Timer"]
time = np.arange(0,6,0.01)

model = models.th_cell_diff

# =============================================================================
# run experiment
# =============================================================================
exp = multi_exp(time, cond_list, cond_names, cond_names2)

g = sns.FacetGrid(exp, hue = "cond2", col = "cond", height = 4)
g = g.map(plt.plot, "time", "value")
g.set(ylim = (0,8))
g.set_titles("{col_name}")
g.add_legend()
g.savefig("death_timer.svg")