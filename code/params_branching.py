#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 09:47:20 2019

@author: burt

parameters for branching models
"""
import numpy as np

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

# =============================================================================
# functions to generate dictionaries from lognormal distributions
# =============================================================================
def convert(logmean, logsigma):
    # deprecated because logsigma is supplied but I use CV
    var = np.log(logsigma**2/logmean**2 + 1)
    sigma = np.sqrt(var)
    mean = np.log(logmean) - var/2
    return [mean, sigma]

def convert2(logmean, CV):
    # get sigma and mean from lognormal mean + CV
    logsigma = logmean*CV
    var = np.log(logsigma**2/logmean**2 + 1)
    sigma = np.sqrt(var)
    mean = np.log(logmean) - var/2
    return [mean, sigma]

def f(logmean, CV):
    # draw from lognormal distribution
    # provide mean and CV of lognormal distribution calculates mean and sigma
    # of normal distribution internally
    mean, sigma = convert2(logmean, CV)
    s = np.random.lognormal(mean, sigma)
    return s

def make_logdict(dic, CV):
    # for all provided keys take a dict with values and return new dict 
    # with vals drawn from lognormal distr. with mean of original dict
    keys = ["beta1", "beta2", "beta1_p", "beta2_p"]
    
    d2 = dict(dic)
    for key in keys:
        d2[key] = f(d2[key], CV)
    
    return d2

def gen_lognorm_dicts(dic, CV, size):
    # make a list of randomly sampled dictionaries from provided dict for given CV
    a = [make_logdict(dic, CV) for i in range(size)]
    return a

