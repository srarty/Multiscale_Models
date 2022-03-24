# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 12:17:10 2022


@author: Artemio Soto-Breceda [artemios]


References:
    [1] BILLEH, Yazan N. et al (2020) Systematic Integration of Structural and 
        Functional Data into Multi-scale Models of Mouse Primary Visual Cortex.

    [2] CAVALLARI, Stefano et al (2014) Comparison of the dynamics of neural 
        interactions between current-based and conductance-based integrate-and-
        fire recurrent networks.
"""
from brian2 import *

def set_params(type='pyramidal', source='brunel'):
    '''
    Parameter space for the different populations of LIF in the modified Brunel model.

    Inputs:
            type:   Neuron type ('pyramidal', 'inhibitory' or 'spiny').
            source: Reference for the parameters ('brunel'-from Callavari 2014, 'allen'-from billeh 2020 )
    '''
    
    if type == 'pyramidal':
        if source == 'brunel':
            g_leak = 25 * nS # Leak conductance
            C = 0.5 * nF # Membrane capacitance
    
            tau_rp = 2 * ms # Absolute refractory period
            tau_m = 20 * ms # Membrane time constant
    
            tau_GABA_r = 0.25 * ms
            tau_GABA_d = 5 * ms
            tau_AMPA_r = 0.4 * ms
            tau_AMPA_d = 2 * ms
            tau_l = 1 * ms # Latency
    
            # Synaptic efficacies
            j_GABA  = 42.5 * pA 
            j_AMPA  = -10.5 * pA
            j_AMPA_ext = -13.75 * pA
            
            # Delta function weight (increment with each input spike)
            # Defined experimentally with 'synaptic_functions.py'. Based on the 
            # unitary increment of the single exponential.
            single_exp_weight = 1 # Inverse: 1/1.52 = 0.6578947368421053
            alpha_simple_weight_AMPA = 12.5
            alpha_simple_weight_GABA = 2.6
            alpha_simple_weight_AMPA_ext = 12.5
            single_exponential_weight = 8.2
            delayed_exp_weight_AMPA = 6.0995 # Inverse 1.52/6.0995 = 0.24920075416017706
            delayed_exp_weight_GABA = 6.0995 # TODO
            
        else: # if 'allen'
            g_leak = 25 * nS # Leak conductance
            C = 0.5 * nF # Membrane capacitance
    
            tau_rp = 2 * ms # Absolute refractory period
            tau_m = 20 * ms # Membrane time constant
    
            tau_GABA_r = 0.25 * ms
            tau_GABA_d = 5 * ms
            tau_AMPA_r = 0.4 * ms
            tau_AMPA_d = 2 * ms
            tau_l = 1 * ms # Latency
    
            # Synaptic efficacies
            j_GABA  = 49.7 * pA 
            j_AMPA  = -147 * pA
            j_AMPA_ext = -13.75 * pA
            
            # Delta function weight (increment with each input spike)
            # Defined experimentally with 'synaptic_functions.py'. Based on the 
            # unitary increment of the single exponential.
            single_exp_weight = 1 # TODO
            alpha_simple_weight_AMPA = 1
            alpha_simple_weight_GABA = 1
            alpha_simple_weight_AMPA_ext = 12.5 # TODO
            single_exponential_weight = 8.2 # TODO
            delayed_exp_weight_AMPA = 6.0995 # TODO
            delayed_exp_weight_GABA = 6.0995 # TODO
        
    elif type == 'inhibitory':
        if source == 'brunel':
            g_leak = 20 * nS # Leak conductance
            C = 0.2 * nF # Membrane capacitance
    
            tau_rp = 1 * ms # Absolute refractory period
            tau_m = 10 * ms # Membrane time constant
    
            tau_GABA_r = 0.25 * ms
            tau_GABA_d = 5 * ms
            tau_AMPA_r = 0.2 * ms
            tau_AMPA_d = 1 * ms
            tau_l = 1 * ms # Latency
    
            # Synaptic efficacies
            j_GABA  = 54 * pA 
            j_AMPA  = -14 * pA
            j_AMPA_ext = -19 * pA
            
            # Delta function weight (increment with each input spike)
            single_exp_weight = 1
            alpha_simple_weight_AMPA = 13 # 2.85
            alpha_simple_weight_GABA = 1.25 # 0.69
            alpha_simple_weight_AMPA_ext = 23
            single_exponential_weight = 8.2
            delayed_exp_weight_AMPA = 27.96 
            delayed_exp_weight_GABA = 39
        else:
            g_leak = 20 * nS # Leak conductance
            C = 0.2 * nF # Membrane capacitance
    
            tau_rp = 1 * ms # Absolute refractory period
            tau_m = 10 * ms # Membrane time constant
    
            tau_GABA_r = 0.25 * ms
            tau_GABA_d = 5 * ms
            tau_AMPA_r = 0.2 * ms
            tau_AMPA_d = 1 * ms
            tau_l = 1 * ms # Latency
    
            # Synaptic efficacies
            j_GABA  = 35.1 * pA 
            j_AMPA  = -896 * pA
            j_AMPA_ext = -19 * pA
            
            # Delta function weight (increment with each input spike)
            single_exp_weight = 1
            alpha_simple_weight_AMPA = 1
            alpha_simple_weight_GABA = 1
            alpha_simple_weight_AMPA_ext = 23
            single_exponential_weight = 8.2
            delayed_exp_weight_AMPA = 27.96 
            delayed_exp_weight_GABA = 39
        
    elif type == 'spiny':
        # TODO
        g_leak = 25 * nS # Leak conductance
        C = 0.5 * nF # Membrane capacitance

        tau_rp = 2 * ms # Absolute refractory period
        tau_m = 20 * ms # Membrane time constant

        tau_GABA_r = 0.25 * ms
        tau_GABA_d = 5 * ms
        tau_AMPA_r = 0.4 * ms
        tau_AMPA_d = 2 * ms
        tau_l = 1 * ms # Latency

        # Synaptic efficacies
        j_GABA  = 42.5 * pA 
        j_AMPA  = -10.5 * pA
        j_AMPA_ext = -13.75 * pA

        # Delta function weight (increment with each input spike)
        single_exp_weight = 1 # Inverse: 1/1.52 = 0.6578947368421053
        alpha_simple_weight_AMPA = 1.52
        alpha_simple_weight_AMPA_ext = 1.52
        single_exponential_weight = 8.2 # TODO
        alpha_simple_weight_GABA = 1.52 
        delayed_exp_weight_AMPA = 6.0995 # Inverse 1.52/6.0995 = 0.24920075416017706
        delayed_exp_weight_GABA = 6.0995 # Inverse 1.52/6.0995 = 0.24920075416017706

    else:
        return 0

    #%% Connectivity
    if source == 'brunel':
        p_IP =  0.2
        p_PI =  0.2
        p_PE =  0.2 
        p_EP =  0.2 
        p_PP =  0.2
        p_II =  0.2
    else:
        p_IP =  0.411
        p_PI =  0.395
        p_PE =  0.2 
        p_EP =  0.2 
        p_PP =  0.16
        p_II =  0.451

    #%% Parse output
    params = {
        "g_leak":       g_leak,
        "C":            C,
        "tau_rp":       tau_rp,
        "tau_m":        tau_m,
        "tau_GABA_r":   tau_GABA_r,
        "tau_GABA_d":   tau_GABA_d,
        "tau_AMPA_r":   tau_AMPA_r,
        "tau_AMPA_d":   tau_AMPA_d,
        "tau_l":        tau_l,
        "j_GABA":       j_GABA,
        "j_AMPA":       j_AMPA,
        "j_AMPA_ext":   j_AMPA_ext,
        "p_IP":         p_IP,
        "p_PI":         p_PI,
        "p_PE":         p_PE,
        "p_EP":         p_EP,
        "p_PP":         p_PP,
        "p_II":         p_II,
        "alpha_weight_AMPA":        alpha_simple_weight_AMPA,
        "alpha_weight_AMPA_ext":    alpha_simple_weight_AMPA_ext,
        "alpha_weight_GABA":        alpha_simple_weight_GABA,
        "single_exp":               single_exponential_weight
    }
    
    return params
    
def set_equations():
    # TODO
    return 0