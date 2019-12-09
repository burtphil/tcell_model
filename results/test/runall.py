#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 10:48:26 2019

@author: burt
"""

import os
os.chdir("/home/burt/Documents/projects/2019/tcell_model/code")

import module_branching as m_branch
from params_branching import d_prec

import numpy as np
t = np.arange(0,10,0.1)

cells = m_branch.run_model(d_prec, t, m_branch.branch_precursor)