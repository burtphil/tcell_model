# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
from tcell_parameters import d_prec

import sys
sys.path.append("/home/burt/Documents/projects/2019/tcell_model/code/")
import matplotlib.pyplot as plt
import matplotlib
from test_module import run_exp, update_dict
import module_models as models
import numpy as np
import seaborn as sns

sns.set(context = "poster", style = "ticks", rc = {"lines.linewidth": 4})

# =============================================================================
# define exp conditions
# =============================================================================
d1 = d_prec
arr1 = np.arange(0.6,0.9,0.03)
name1 = "p1_def"
conditions1 = [update_dict(d1, val ,name1) for val in arr1]
cond_names1 = [name1+"="+str(round(val,3)) for val in arr1]

d2 = d_prec
arr2 = np.arange(5,15.5,0.5)
name2 = "beta1"
conditions2 = [update_dict(d2, val ,name2) for val in arr2]
cond_names2 = [name2+"="+str(round(val,3)) for val in arr2]


time = np.arange(0,6,0.01)

model = models.branch_precursor
model_name = d_prec["mode"]
# =============================================================================
# run experiment
# =============================================================================
exp1 = run_exp(time, conditions1, cond_names1, model = model)
exp2 = run_exp(time, conditions2, cond_names2, model = model)

norm1 = matplotlib.colors.Normalize(
    vmin=np.min(arr1),
    vmax=np.max(arr1))

norm2 = matplotlib.colors.Normalize(
    vmin=np.min(arr2),
    vmax=np.max(arr2))

cm = matplotlib.cm.Blues
sm1 = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm1)
sm1.set_array([])

sm2 = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm2)
sm2.set_array([])


fig, (ax1, ax2) = plt.subplots(1,2, figsize = (14,5))

g = sns.lineplot(x = "time", y = "value", style = "cell", hue = "cond",
                 palette = "Blues", data = exp1, ax = ax1)

g.get_legend().remove()
cbar1 = fig.colorbar(sm1, ax = ax1, ticks = [0.6, 0.7, 0.8, 0.9])
cbar1.set_label("$p_1$")
g.set(ylabel = "cell dens. norm. (a.u.)")

f = sns.lineplot(x = "time", y = "value", style = "cell", hue = "cond",
                 palette = "Blues", data = exp2, ax = ax2)

f.get_legend().remove()
cbar2 = fig.colorbar(sm2, ax = ax2, ticks = [5, 10, 15])
cbar2.set_label(r"$\beta_1$")
f.set(ylabel = "cell dens. norm. (a.u.)")

plt.tight_layout()
fig.savefig("branching_timecourse_color_"+model_name+".pdf")
fig.savefig("branching_timecourse_colorbar_"+model_name+".svg")

