# -*- coding: utf-8 -*-

d_null = {
        "b" : 0,
        "alpha_int" : 5,
        "alpha" : 10,
        "beta" : 10,
        "alpha_p" : 10,
        "beta_p" : 30.,
        "d_eff" : 1.5,
        "d_naive" : 0,
        "d_prec" : 0,
        "mode" : "Null",
        "rate_il2" : 100,
        "rate_il7" : 5,
        "fb_il2" : 1.0,
        "rate_ifn" : 1.0,
        "K_ifn" : 1.0,
        "fb_ifn" : 0,
        "crit" : False,
        "crit_timer" : 4,
        "crit_il7" : 0.1,
        "crit_il2" : 0.1,
        "t0" : None,
        "alpha_IL2" : 8,
        "decay_p" : 1.,
        "death_mode" : True,
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
