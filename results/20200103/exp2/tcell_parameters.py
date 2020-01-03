# -*- coding: utf-8 -*-
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
          "mode" : "Null",
          "crit" : False,
          "crit_il2" : 0.5,
          "crit_il7" : 10,
          "crit_timer" : 2,
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
