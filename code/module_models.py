# -*- coding: utf-8 -*-
"""
keep all ode models here
"""

import numpy as np

# =============================================================================
# linear models
# =============================================================================
def th_cell_diff(th_state, time, d):
    """
    model2
    takes state vector to differentiate effector cells as linear chain
    needs alpha and beta(r) of response time distribution, probability
    and number of precursor cells
    """
    #tnaive = np.sum(th_state[:-(2*d["alpha_p"])])
    teff = np.sum(th_state[-(2*d["alpha_p"]):-d["alpha_p"]])
    # this is the beta sad IL2 population
    tnoil2 = np.sum(th_state[-d["alpha_p"]:])
    
    # get IL2 producers (depending on alpha IL2)
    t_il2 = np.sum(th_state[:d["alpha_IL2"]])
    t_noil2 = np.sum(th_state)-t_il2
    conc_il2 = d["rate_il2"]*t_il2/(0.001+t_noil2)

    #carrying capacity
    il7_consumers = teff+tnoil2

    # mm kinetic feedback implementation  
    fb_ifn = 0
    if d["fb_ifn"] != 0:
        conc_ifn = d["rate_ifn"]*(teff+tnoil2)
        fb_ifn = (d["fb_ifn"]*conc_ifn**3)/(conc_ifn**3+d["K_ifn"]**3)
  
    beta = (fb_ifn+1)*d["beta"]
    
    beta_p = d["beta_p"]
    rate_death = d["d_eff"]

    # check if crit has been updated from F to T, if not
    # check if crit should be updated, then store t0 time of update
    # if crit has been updated, update prolif rate
    
    if d["mode"] in ["il7", "il2", "timer", "il2_timer", "il2+", "timer+"]:

        if d["crit"] == True:
            if d["death_mode"] == True:
                rate_death = rate_death*np.exp(time-d["t0"])
            else:
                beta_p = beta_p*np.exp(-d["decay_p"]*(time-d["t0"]))
        else:
            # define criteria upon which apoptosis is induced
            crit_time = time > d["crit_timer"]
            crit_il2 = conc_il2 < d["crit_il2"]
            crit_il7 = il7_consumers > d["crit_il7"]
            
            c1 = d["mode"] == "timer" and crit_time
            c2 = d["mode"] == "il2" and crit_il2
            c3 = d["mode"] == "il7" and crit_il7
            c4 = d["mode"] == "timer+" and (crit_time or crit_il7)
            c5 = d["mode"] == "il2+" and (crit_il2 or crit_il7)
            c6 = d["mode"] == "il2_timer" and (crit_il2 or crit_time)
 
            crits = np.array([c1,c2,c3,c4,c5,c6])
            if crits.any():
                d["t0"] = time
                d["crit"] = True 
                
    dt_state = diff_effector(th_state, teff, d, beta, rate_death, beta_p)
        
 
    return dt_state

def diff_effector(th_state, teff, d, beta, rate_death, beta_p):
    dt_state = np.zeros_like(th_state)
    #print(th_state.shape)

    for j in range(len(th_state)):
        #print(j)
        if j == 0:
            dt_state[j] = d["b"]-beta*th_state[j] 
            
        elif j < d["alpha"]:
            dt_state[j] = beta*th_state[j-1]-(beta+d["d_prec"])*th_state[j]
            
        elif j == d["alpha"]:
            dt_state[j] = beta*th_state[j-1] + (2*beta_p*th_state[(d["alpha"]+d["alpha_p"]-1)]) - (rate_death+d["beta_sad"]+beta_p)*th_state[j]       

        elif j < (d["alpha"]+d["alpha_p"]):
            dt_state[j] = beta_p*th_state[j-1]-(beta_p+rate_death+d["beta_sad"])*th_state[j]   
        
        elif j == (d["alpha"]+d["alpha_p"]):
            dt_state[j] = d["beta_sad"]*teff + (2*beta_p*th_state[-1]) - (rate_death+beta_p)*th_state[j]       
        
        else:
            dt_state[j] = beta_p*th_state[j-1]-(beta_p+rate_death)*th_state[j]
        
    return dt_state

# =============================================================================
# branching models
# =============================================================================
def diff_effector2(state, th0, alpha, beta, beta_p, p, d):
    """
    takes state vector to differentiate effector cells as linear chain
    needs alpha and beta(r) of response time distribution, probability
    and number of precursor cells
    """
    dt_state = np.zeros_like(state)
    #print(len(state))
    if alpha == 1:
        for j in range(len(state)):
            if j == 0:
                dt_state[j] = p*beta*th0+2*beta_p*state[-1]-(beta_p+d["d_eff"])*state[j]
            else:
                dt_state[j] = beta_p*state[j-1]- (beta_p+d["d_eff"])*state[j] 
    
    else:                    
        for j in range(len(state)):
            if j == 0:
                dt_state[j] = p*beta*th0 - (beta+d["d_prec"])*state[j]                
            elif j < (alpha-1):
                dt_state[j] = beta*state[j-1]-(beta+d["d_prec"])*state[j]                
            elif j == (alpha-1):
                # the problem with the 4 and 2 is that since differentiation takes 1 day it should divide twice giving 4 cells
                # however, if it has arrived in the final states if should double every half day
                dt_state[j] = beta*state[j-1]+2*beta_p*state[-1] - (d["d_eff"]+beta_p)*state[j]  

            else:
                assert j >= alpha
                dt_state[j] = beta_p*state[j-1]- (beta_p+d["d_eff"])*state[j] 
               
    return dt_state

def diff_precursor(state, th0, alpha, beta, beta_p, p_adj, d):
    """
    takes state vector to differentiate effector cells as linear chain
    needs alpha and beta(r) of response time distribution, probability
    and number of precursor cells
    """
    dt_state = np.zeros_like(state)

    for j in range(len(state)):
        if j == 0:
            dt_state[j] = p_adj*beta*th0 - beta*state[j]             
        elif j < alpha:
            dt_state[j] = beta*state[j-1]- beta*state[j]                
        elif j == alpha:
            # the problem with the 4 and 2 is that since differentiation takes 1 day it should divide twice giving 4 cells
            # however, if it has arrived in the final states if should double every half day
            dt_state[j] = beta*state[j-1]+2*beta_p*state[-1] - (d["d_eff"]+beta_p)*state[j] 
            
        else:
            assert j > alpha        
            dt_state[j] = beta_p*state[j-1]-(beta_p+d["d_eff"])*state[j] 
                 
    return dt_state

def branch_competetive(state, time, d):
    """
    takes state vector to differentiate effector cells as linear chain
    needs alpha and beta(r) of response time distribution, probability
    and number of precursor cells
    """

    th0 = state[0]    
    th1 = state[1:(d["alpha1"]+d["alpha1_p"])]
    th2 = state[(d["alpha1"]+d["alpha1_p"]):]
    
    #print(len(state), len(th1))
    ### get all cytokine secreting cells    
    th1_all = np.sum(th1[-d["alpha1_p"]:])
    th2_all = np.sum(th2[-d["alpha2_p"]:])
    ### calculate cytokine concentrations
    cyto_1 = d["beta_cyto_1"]*th1_all + d["ifn_ext"]
    cyto_2 = d["beta_cyto_2"]*th2_all + d["il21_ext"]
    ### calculate cytokine effect on rate
    fb1 = d["fb_rate1"]*cyto_1**3/(cyto_1**3+d["K_1"]**3)
    fb2 = d["fb_rate2"]*cyto_2**3/(cyto_2**3+d["K_2"]**3)
    ### update differantiation rate
    beta1 = d["beta1"]*(1+fb1)
    beta2 = d["beta2"]*(1+fb2)
    
    ### differentiate effectors th1    
    alpha = d["alpha1"]
    p = 1.
    dt_th1 = diff_effector2(th1, th0, alpha, beta1, d["beta1_p"], p, d)
    ### differentiate effectors th2
    alpha = d["alpha2"]
    p = 1.
    dt_th2 = diff_effector2(th2, th0, alpha, beta2, d["beta2_p"], p, d)
    
    ### combine all cells
    dt_th0 = -(beta1+beta2)*th0
    dt_state = np.concatenate(([dt_th0], dt_th1, dt_th2))

    return dt_state

def branch_precursor(state, time, d):
    """
    takes state vector to differentiate effector cells as linear chain
    needs alpha and beta(r) of response time distribution, probability
    and number of precursor cells
    """

    th0 = state[0]
    
    th1 = state[1:(d["alpha1"]+d["alpha1_p"]+1)]
    th2 = state[(d["alpha1"]+d["alpha1_p"]+1):]
    #print(len(state), len(th1))
    ### get all cytokine secreting cells    
    th1_all = np.sum(th1[-d["alpha1_p"]:])
    th2_all = np.sum(th2[-d["alpha2_p"]:])
    ### calculate cytokine concentrations
    cyto_1 = d["beta_cyto_1"]*th1_all + d["ifn_ext"]
    cyto_2 = d["beta_cyto_2"]*th2_all + d["il21_ext"]
    
    ### calculate probability 
    p1 = d["fb_prob1"]*cyto_1**3/(cyto_1**3+d["K_1"]**3)
    p2 = d["fb_prob2"]*cyto_2**3/(cyto_2**3+d["K_2"]**3)
    # account for default probability and feedback strength

    p1 = (p1+1)*d["p1_def"]
    p2 = (p2+1)*(1-d["p1_def"])

    ### this is the effective probability after feedback integration 
    p1_norm = p1/(p1+p2)
    #p2_norm = 1-p1_norm
    p2_adj = 1.

    ### calculate fb rate effects
    fb1 = d["fb_rate1"]*cyto_1**3/(cyto_1**3+d["K_1"]**3)
    fb2 = d["fb_rate2"]*cyto_2**3/(cyto_2**3+d["K_2"]**3)
    ### update differantiation rate
    beta1 = d["beta1"]*(1+fb1)
    beta2 = d["beta2"]*(1+fb2)
      
    # adjust this parameter to effectively change branching prob because beta1 and beta2 also
    # play a role, note that fb can affect beta1,2
    p1_adj = p1_norm*beta2/(beta1*(1-p1_norm))
    #print(beta1*p1_adj/(beta1*p1_adj+beta2))
    dt_th1 = diff_precursor(th1, th0, d["alpha1"], beta1, d["beta1_p"], p1_adj, d)
    dt_th2 = diff_precursor(th2, th0, d["alpha2"], beta2, d["beta2_p"], p2_adj, d)

    dt_th0 = -(beta1*p1_adj+beta2)*th0    
    dt_state = np.concatenate(([dt_th0], dt_th1, dt_th2))

    return dt_state