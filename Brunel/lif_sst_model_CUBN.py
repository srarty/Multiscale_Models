# -*- coding: utf-8 -*-
"""
Created on Dec 12

Model parameters for Current Based LIF

@author: Artemio Soto-Breceda [artemios]


References:
    [1] BILLEH, Yazan N. et al (2020) Systematic Integration of Structural and 
        Functional Data into Multi-scale Models of Mouse Primary Visual Cortex.

    [2] CAVALLARI, Stefano et al (2014) Comparison of the dynamics of neural 
        interactions between current-based and conductance-based integrate-and-
        fire recurrent networks.
"""
from brian2 import *

def set_params(type='pyramidal'):
    '''
    Parameter space for the different populations of LIF in the modified Brunel model.

    Inputs:
            type:   Neuron type ('pyramidal', 'inhibitory').
    '''
    
    #%% Default parameters
    # Connectivity
    p_IP =  0.411 
    p_PI =  0.395 
    p_PP =  0.16 
    p_II =  0.451 
    p_SP = 0.351 # Billeh
    p_PS = 0.571
    p_SI = 0.857
    
    if type == 'pyramidal':
        #%% Pyramidal. Allen
        g_leak = 25 * nS # Leak conductance
        C = 0.5 * nF # Membrane capacitance
        
        tau_rpa = 2 * ms # Absolute refractory period
        tau_rpr = 20 * ms # Relative refractory period
        tau_rp = tau_rpa + tau_rpr # Effective refractory period
        tau_m = 20 * ms # Membrane time constant

        tau_GABA_s = 5.25 * ms
        tau_GABAb_s = 105 * ms
        tau_AMPA_s = 2.4 * ms
        tau_AMPA_s_ext = 2.4 * ms
        tau_l = 1 * ms # Latency

        # Synaptic efficacies
        j_GABA  = 97.597 * pA
        j_GABAb = 97.597 * pA # 2.3125 * pA
        j_GABAs = 43.92 * pA # 2.3125 * pA
        j_AMPA  = -230.01 * pA#-36.75 * pA
        j_AMPA_ext = -112.75 * pA #1.375 * pA # Too large
        j_AMPA_tha = -112.75 * pA # too large_
         
       
        # Delta function weight (increment with each input spike)
        # Defined experimentally with 'synaptic_functions.py'. Based on the 
        # unitary increment of the single exponential.
        weight = 1
        external_input_weight = 1.15
        
    elif type == 'inhibitory':            
        #%% Parvalbumin positive (Fast Spiking Basket cells)
        g_leak = 8.54 * nS # (<- allen institute) #20 * nS # Leak conductance
        C = 61.63 * pF # (<- allen institute) # 0.2 * nF # Membrane capacitance

        tau_rpa = 1 * ms # Absolute refractory period
        tau_rpr = 12.5 * ms # Relative refractory period
        tau_rp = tau_rpa + tau_rpr # Effective refractory period
        tau_m = C/g_leak #10 * ms # Membrane time constant

        tau_GABA_s = 2.6 * ms # Galarreta 2002 # 5.25 * ms
        tau_GABAb_s = 0 * ms
        tau_AMPA_s = 1.2 * ms
        tau_AMPA_s_ext = 1.2*ms 
        tau_l = 1 * ms # Latency

        # Synaptic efficacies
        j_GABA  = 28.4 * pA #45.24 * pA
        j_GABAb  = 0 * pA
        j_GABAs  = 21.29 * pA
        j_AMPA  = -94.4625 * pA
        j_AMPA_ext = -1.9 * pA
        j_AMPA_tha = -19 * pA
        
        # Delta function weight (increment with each input spike)
        weight = 1
        external_input_weight = 8.2
        
    elif type == 'sst':
        #%% SST (Non-martinoti)
        g_leak = 1/(132 * Mohm) # Leak conductance
        C = 62.12 * pF # Membrane capacitance

        tau_rpa = 2 * ms # Absolute refractory period
        tau_rpr = 20 * ms # 100 * ms # Relative refractory period
        tau_rp = tau_rpa + tau_rpr # Effective refractory period
        tau_m = C/g_leak #10 * ms # Membrane time constant

        tau_GABA_s = 5.25 * ms # To Do
        tau_GABAb_s = 0 * ms # To Do
        tau_AMPA_s = 1.2 * ms # To Do
        tau_AMPA_s_ext = 1.2*ms  # To Do
        tau_l = 1 * ms # Latency

        # Synaptic efficacies
        j_GABA  = 0 * pA # 22.62 * pA #45.24 * pA # To Do
        j_GABAb  = 0 * pA
        j_GABAs  = 0 * pA
        j_AMPA  = -36.7125 * pA
        j_AMPA_ext = 0 * pA
        j_AMPA_tha = 0 * pA
        
        # Delta function weight (increment with each input spike)
        weight = 1
        external_input_weight = 8.2
        
    else:
        return 0


    #%% Parse output
    params = {
        "g_leak":       g_leak,
        "C":            C,
        "tau_rp":       tau_rp,
        "tau_m":        tau_m,
        "tau_GABA_s":   tau_GABA_s,
        "tau_GABAb_s":   tau_GABAb_s,
        "tau_AMPA_s":   tau_AMPA_s,
        "tau_AMPA_s_ext":   tau_AMPA_s_ext,
        "tau_l":        tau_l,
        "j_GABA":       j_GABA,
        "j_GABAb":      j_GABAb,
        "j_GABAs":       j_GABAs,
        "j_AMPA":       j_AMPA,
        "j_AMPA_ext":   j_AMPA_ext,
        "j_AMPA_tha":   j_AMPA_tha,
        "p_IP":         p_IP,
        "p_PI":         p_PI,
        "p_PP":         p_PP,
        "p_II":         p_II,
        "p_PS":         p_PS,
        "p_SP":         p_SP,
        "p_SI":         p_SI,
        "weight":       weight,
        "external_input_weight":               external_input_weight
    }
    
    return params
    
#%% Equations
def get_equations(type = 'pyramidal'):
    if type == 'pyramidal':
        eqs = '''
            dv / dt = (-v + V_leak - (I_tot/g_m_P)) / tau_m_P : volt (unless refractory)
            
            dv_pi /dt = (-v_pi - ( I_GABA_rec / g_m_P)) / tau_m_P : volt (unless refractory)
            dv_pp /dt = (-v_pp - ( I_AMPA_rec / g_m_P)) / tau_m_P : volt (unless refractory)
            dv_pb /dt = (-v_pb - ( I_GABAb    / g_m_P)) / tau_m_P : volt (unless refractory)
            dv_pu /dt = (-v_pu - ( I_AMPA_cor / g_m_P)) / tau_m_P : volt (unless refractory)
            dv_ps /dt = (-v_ps - ( I_GABAs    / g_m_P)) / tau_m_P : volt (unless refractory)
        
            I_tot = I_exc + I_inh + I_injected : amp
            I_exc = I_AMPA_cor + I_AMPA_tha + I_AMPA_rec : amp
            I_inh = I_GABAb + I_GABAs + I_GABA_rec : amp
            
            I_AMPA_cor = j_AMPA_cor_P * s_AMPA_cor : amp
            ds_AMPA_cor / dt = -s_AMPA_cor / tau_s_AMPA_P : 1    
            
            I_AMPA_tha = j_AMPA_tha_P * s_AMPA_tha : amp
            ds_AMPA_tha / dt = -s_AMPA_tha / tau_s_AMPA_P : 1
            
            I_GABA_rec = j_GABA_P * s_GABA : amp
            ds_GABA / dt = -s_GABA / tau_s_GABA_P : 1

            I_GABAb = j_GABAb_P * s_GABAb: amp
            ds_GABAb / dt = -s_GABAb / tau_s_GABAb_P : 1
            
            I_GABAs = j_GABAs_P * s_GABAs: amp
            ds_GABAs / dt = -s_GABAs / tau_s_GABA_P : 1
            
            I_AMPA_rec = j_AMPA_rec_P * s_AMPA : amp
            ds_AMPA / dt = -s_AMPA / tau_s_AMPA_P: 1    
            
            ref : second
            v_th : volt
        '''
        
        
    elif type == 'inhibitory':
        eqs = '''
            dv / dt = (-v + V_leak - (I_tot/g_m_I)) / tau_m_I : volt (unless refractory)
        
            dv_ip /dt = (-v_ip -(I_AMPA_rec / g_m_I)) / tau_m_I : volt (unless refractory)
            dv_iu /dt = (-v_iu -(I_AMPA_cor / g_m_I)) / tau_m_I : volt (unless refractory)
            dv_ii /dt = (-v_ii -(I_GABA_rec / g_m_I)) / tau_m_I : volt (unless refractory)
            dv_is /dt = (-v_is -(I_GABAs    / g_m_I)) / tau_m_I : volt (unless refractory)
            
            I_tot = I_exc + I_inh + I_injected_I: amp
            I_exc = I_AMPA_cor + I_AMPA_tha + I_AMPA_rec : amp
            I_inh = I_GABAs + I_GABA_rec : amp
            
            I_AMPA_cor = j_AMPA_cor_I * s_AMPA_cor : amp
            ds_AMPA_cor / dt = - s_AMPA_cor / tau_s_AMPA_I : 1
            
            I_AMPA_tha = j_AMPA_tha_I * s_AMPA_tha : amp
            ds_AMPA_tha / dt = - s_AMPA_tha / tau_s_AMPA_I_ext : 1
            
            I_GABA_rec = j_GABA_I * s_GABA : amp
            ds_GABA / dt = -s_GABA / tau_s_GABA_I : 1
            
            I_GABAs = j_GABA_I * s_GABAs : amp
            ds_GABAs / dt = -s_GABAs / tau_s_GABAs_I : 1
            
            I_AMPA_rec = j_AMPA_rec_I * s_AMPA : amp
            ds_AMPA / dt = -s_AMPA / tau_s_AMPA_I : 1
            
            ref : second
            v_th : volt
        '''
        
    elif type == 'sst':
        eqs = '''
            dv / dt = (-v + V_leak - (I_tot/g_m_I)) / tau_m_I : volt (unless refractory)
        
            dv_sp /dt = (-v_sp -(I_AMPA_rec / g_m_I)) / tau_m_I : volt (unless refractory)
            
            I_tot = I_exc + I_inh: amp
            I_exc = I_AMPA_rec : amp
            I_inh = 0 * amp : amp
            
            I_AMPA_rec = j_AMPA_rec_S * s_AMPA : amp
            ds_AMPA / dt = -s_AMPA / tau_s_AMPA_S : 1
            
            ref : second
            v_th : volt
        '''
        
    else:
        raise ValueError(format('The option type = %s is not a valid one.' %(type)))
    return eqs