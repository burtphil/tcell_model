# -*- coding: utf-8 -*-

d_null = {
        "b" : 0,
        "alpha" : 10,
        "beta" : 10.,
        "alpha_p" : 10,
        "beta_p" : 10.,
        "beta_sad" : 0,
        "d_eff" : 1.5,
        "d_prec" : 0,
        "mode" : "Null",
        "rate_il2" : 100.0,
        "fb_il2" : 1.0,
        "rate_ifn" : 1.0,
        "K_ifn" : 1.0,
        "fb_ifn" : 0,
        "crit" : False,
        "crit_timer" : 1.5,
        "crit_il7" : 2.5,
        "crit_il2" : 0.5,
        "t0" : None,
        "alpha_IL2" : 4,
        "decay_p" : 1.,
        "death_mode" : False,
        "K_il2": 0.001
        }

d_il7 = dict(d_null)
d_il7["mode"] = "il7"

d_il2 = dict(d_null)
d_il2["mode"] = "il2"

d_timer = dict(d_null)
d_timer["mode"] = "timer"

d_il7_death = dict(d_il7)
d_il2_death = dict(d_il2)
d_timer_death = dict(d_timer)

for dic in [d_il7_death, d_il2_death, d_timer_death]:
    dic["death_mode"] = True