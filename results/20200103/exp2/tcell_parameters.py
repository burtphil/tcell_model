# -*- coding: utf-8 -*-
d_comp = {
          "alpha1" : 9, 
          "alpha2" : 9, 
          "beta1" : 10.,
          "beta2" : 10.,
          "alpha1_p" : 10,
          "alpha2_p" : 10,
          "beta1_p" : 30,
          "beta2_p" : 30,  
          "fb_rate1" : 0,
          "fb_rate2" : 0,
          "fb_prob1" : 0,
          "fb_prob2" : 0,
          "K_1" : 1.,
          "K_2" : 1.,
          "beta_cyto_1": 1.,
          "beta_cyto_2": 1.,
          "d_eff" : 1.0,
          "d_prec" : 0,
          "ifn_ext" : 0,
          "il21_ext" : 0,
          "mode" : "Null",
          "crit" : False,
          "crit_il2" : 0.5,
          "crit_il7" : 2.5,
          "crit_timer" : 1.5,
          "rate_il2" : 10000,
          "alpha_IL2" : 7,
          "K_il2" : 0.0001,
          "death_mode" : False,
          "decay_p" : 1.,
          "t0" : None
          }

# calculate corresponding beta values for precursor model
d_prec = dict(d_comp)
d_prec["p1_def"] = 0.5

d_prec_il7 = dict(d_prec)
d_prec_il7["mode"] = "il7"

d_prec_il2 = dict(d_prec)
d_prec_il2["mode"] = "il2"

d_prec_timer = dict(d_prec)
d_prec_timer["mode"] = "timer"