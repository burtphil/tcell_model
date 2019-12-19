# -*- coding: utf-8 -*-

d_null = {
        "b" : 0,
        "alpha" : 10,
        "beta" : 10.,
        "alpha_p" : 10,
        "beta_p" : 30.,
        "beta_sad" : 0,
        "d_eff" : 1.5,
        "d_prec" : 0,
        "mode" : "timer",
        "rate_il2" : 10000.0,
        "fb_il2" : 1.0,
        "rate_ifn" : 1.0,
        "K_ifn" : 1.0,
        "fb_ifn" : 0,
        "crit" : False,
        "crit_timer" : 2.,
        "crit_il7" : 10.,
        "crit_il2" : 0.5,
        "t0" : None,
        "alpha_IL2" : 4,
        "decay_p" : 1.,
        "death_mode" : False,
        }


d_comp = {
          "alpha1" : 9, 
          "alpha2" : 9, 
          "beta1" : 10.,
          "beta2" : 10.,
          "alpha1_p" : 10,
          "alpha2_p" : 10,
          "beta1_p" : 10,
          "beta2_p" : 10,  
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
          }


# calculate corresponding beta values for precursor model
d_prec = dict(d_comp)
d_prec["p1_def"] = 0.5