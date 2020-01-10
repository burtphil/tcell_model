# -*- coding: utf-8 -*-

d_null = {
        "b" : 0,
        "alpha" : 10,
        "beta" : 10,
        "alpha_p" : 10,
        "beta_p" : 40.,
        "beta_sad" : 0,
        "d_eff" : 1.5,
        "d_prec" : 0,
        "mode" : "Null",
        "rate_il2" : 1000000.0,
        "rate_il7" : 50.0,
        "fb_il2" : 1.0,
        "rate_ifn" : 1.0,
        "K_ifn" : 1.0,
        "fb_ifn" : 0,
        "crit" : False,
        "crit_space" : 13000,
        "crit_timer" : 0.85,
        "crit_il7" : 0.1,
        "crit_il2" : 0.1,
        "t0" : None,
        "alpha_IL2" : 8,
        "decay_p" : 1.,
        "death_mode" : False,
        "K_il2": 0.0001
        }

d_il7 = dict(d_null)
d_il7["mode"] = "il7"

d_il2 = dict(d_null)
d_il2["mode"] = "il2"
d_il2["beta_sad"] = 10

d_timer = dict(d_null)
d_timer["mode"] = "timer"

d_il2_timer = dict(d_null)
d_il2_timer["mode"] = "il2_timer"