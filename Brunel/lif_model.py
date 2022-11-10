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
    
    #%% Default parameters
    # Connectivity
    p_IP =  0.411 #* 0.25
    p_PI =  0.395 #* 0.25
    p_PP =  0.16 #* 0.25
    p_II =  0.451 #* 0.25
    
    if type == 'pyramidal':
        #%% Pyramidal. Allen
        g_leak = 25 * nS # Leak conductance
        C = 0.5 * nF # Membrane capacitance
        
        tau_rpa = 2 * ms # Absolute refractory period
        tau_rpr = 20 * ms #100 * ms # Relative refractory period
        tau_rp = tau_rpa + tau_rpr # Effective refractory period
        tau_m = 20 * ms # Membrane time constant

        tau_GABA_s = 5.25 * ms
        tau_GABAb_s = 105 * ms
        tau_AMPA_s = 2.4 * ms
        tau_AMPA_s_ext = 2.4 * ms
        tau_l = 1 * ms # Latency

        # Synaptic efficacies
        j_GABA  = 37 * pA # 37 * pA <- (normal parameters) # 18.3 * pA <- corresponds to alpha_i=0.26 # (default = 37 | seizure = j(alpha=-0.3) = 21.0666) 
        j_GABAb  = 0 * pA
        j_AMPA  = -73.5 * pA
        j_AMPA_ext = -1.375 * pA
        j_AMPA_tha = -13.75 * pA
        
        # Delta function weight (increment with each input spike)
        # Defined experimentally with 'synaptic_functions.py'. Based on the 
        # unitary increment of the single exponential.
        weight = 1
        external_input_weight = 8.2
        
    elif type == 'inhibitory':            
        #%% Inhibitory. Allen
        g_leak = 20 * nS # Leak conductance
        C = 0.2 * nF # Membrane capacitance

        tau_rpa = 1 * ms # Absolute refractory period
        tau_rpr = 12.5 * ms # 100 * ms # Relative refractory period
        tau_rp = tau_rpa + tau_rpr # Effective refractory period
        tau_m = C/g_leak #10 * ms # Membrane time constant

        tau_GABA_s = 5.25 * ms
        tau_GABAb_s = 0 * ms
        tau_AMPA_s = 1.2 * ms
        tau_AMPA_s_ext = 1.2*ms 
        tau_l = 1 * ms # Latency

        # Synaptic efficacies
        j_GABA  = 45.24 * pA 
        j_GABAb  = 0 * pA
        j_AMPA  = -165 * pA 
        j_AMPA_ext = -1.9 * pA
        j_AMPA_tha = -19 * pA
        
        # Delta function weight (increment with each input spike)
        weight = 1
        external_input_weight = 8.2
        
    elif type == 'gabab':            
        #%% GABAb. Allen
        g_leak = 20 * nS # Leak conductance
        C = 0.2 * nF # Membrane capacitance

        tau_rpa = 1 * ms # Absolute refractory period
        tau_rpr = 12.5 * ms # 100 * ms # Relative refractory period
        tau_rp = tau_rpa + tau_rpr # Effective refractory period
        tau_m = C/g_leak #10 * ms # Membrane time constant

        tau_GABA_s = 5.25 * ms
        tau_GABAb_s = 0 * ms
        tau_AMPA_s = 1.2 * ms
        tau_AMPA_s_ext = 1.2*ms 
        tau_l = 1 * ms # Latency

        # Synaptic efficacies
        j_GABA  = 45.24 * pA 
        j_GABAb  = 0 * pA
        j_AMPA  = -165 * pA
        j_AMPA_ext = -1.9 * pA
        j_AMPA_tha = -19 * pA
        
        # Delta function weight (increment with each input spike)
        weight = 1
        external_input_weight = 8.2
        
    else:
        return 0
    
    if source == 'allen':
        #%% Allen parameters
        # Same as default, so: do nothing
        print('allen parameters loaded')
    elif source == 'three_pop':
        #%% Three population parameters
        if type == 'pyramidal':            
            j_GABA  = 18.5 * pA
            j_GABAb = 2.3125 * pA
            j_AMPA  = -36.75 * pA
            
            j_AMPA_ext = -112.75 * pA #1.375 * pA # Too large
            j_AMPA_tha = -112.75 * pA # too large_
            external_input_weight = 1.15
            
            
        elif type == 'inhibitory':   
            j_AMPA  = -41.25 * pA
            j_AMPA_tha = -19 * pA
            
        elif type == 'gabab':
            j_GABA  = 4.524 * pA
            j_AMPA  = -41.25 * pA
            j_AMPA_tha = -15.2 * pA
            
            
    elif source == 'brunel':
        #%% Brunel parameters
        # Connectivity
        p_IP =  0.2
        p_PI =  0.2
        p_PP =  0.2
        p_II =  0.2
        
        if type == 'pyramidal':
            #%% Pyramidal. Brunel
            g_leak = 25 * nS # Leak conductance
            C = 0.5 * nF # Membrane capacitance
    
            tau_rp = 2 * ms # Absolute refractory period
            tau_m = 20 * ms # Membrane time constant
    
            tau_GABA_s = 5.25 * ms
            tau_AMPA_s = 2.4 * ms
            tau_AMPA_s_ext = 2.4 * ms
            tau_l = 1 * ms # Latency
    
            # Synaptic efficacies
            j_GABA  = 42.5 * pA 
            j_AMPA  = -10.5 * pA
            j_AMPA_ext = -13.75 * pA
            j_AMPA_tha = -13.75 * pA
            
            # Delta function weight (increment with each input spike)
            # Defined experimentally with 'synaptic_functions.py'. Based on the 
            # unitary increment of the single exponential.
            alpha_simple_weight_AMPA = 12.5
            alpha_simple_weight_GABA = 1
            alpha_simple_weight_AMPA_ext = 12.5
            single_exponential_weight = 8.2
            
        elif type == 'inhibitory':    
            #%% Inhibitory. Brunel
            g_leak = 20 * nS # Leak conductance
            C = 0.2 * nF # Membrane capacitance
    
            tau_rp = 1 * ms # Absolute refractory period
            tau_m = 10 * ms # Membrane time constant
    
            tau_GABA_s = 5.25 * ms
            tau_AMPA_s = 1.2 * ms
            tau_AMPA_s_ext = 1.2 * ms
            tau_l = 1 * ms # Latency
    
            # Synaptic efficacies
            j_GABA  = 54 * pA 
            j_AMPA  = -14 * pA
            j_AMPA_ext = -19 * pA
            j_AMPA_tha = -19 * pA
            
            # Delta function weight (increment with each input spike)
            alpha_simple_weight_AMPA = 13 # 2.85
            alpha_simple_weight_GABA = 1.25 # 0.69
            alpha_simple_weight_AMPA_ext = 23
            single_exponential_weight = 8.2
            
        else:
            return 0
        
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
        "j_AMPA":       j_AMPA,
        "j_AMPA_ext":   j_AMPA_ext,
        "j_AMPA_tha":   j_AMPA_tha,
        "p_IP":         p_IP,
        "p_PI":         p_PI,
        "p_PP":         p_PP,
        "p_II":         p_II,
        "weight":       weight,
        "external_input_weight":               external_input_weight
    }
    
    return params
    
#%% Equations
def get_equations(type = 'pyramidal'):
    if type == 'pyramidal':
        eqs = '''
            dv / dt = (-v + V_leak - (I_tot/g_m_P)) / tau_m_P : volt (unless refractory)
            
            dv_pe /dt = (-v_pe - ((I_AMPA_spi + I_AMPA_cor) / g_m_P)) / tau_m_P : volt (unless refractory)
            dv_pi /dt = (-v_pi - ( I_GABA_rec               / g_m_P)) / tau_m_P : volt (unless refractory)
            dv_pp /dt = (-v_pp - ( I_AMPA_rec               / g_m_P)) / tau_m_P : volt (unless refractory)
            dv_pb /dt = (-v_pb - ( I_GABAb               / g_m_P)) / tau_m_P : volt (unless refractory)
            dv_pu /dt = (-v_pu - ( I_AMPA_cor               / g_m_P)) / tau_m_P : volt (unless refractory)
        
            I_tot = I_AMPA_cor + I_AMPA_tha + I_AMPA_rec + I_AMPA_spi + I_GABAb + I_GABA_rec + I_injected : amp
            
            I_AMPA_cor = j_AMPA_cor_P * s_AMPA_cor : amp
            ds_AMPA_cor / dt = -s_AMPA_cor / tau_s_AMPA_P : 1    
            
            I_AMPA_tha = j_AMPA_tha_P * s_AMPA_tha : amp
            ds_AMPA_tha / dt = -s_AMPA_tha / tau_s_AMPA_P : 1
            
            I_GABA_rec = agonist * j_GABA_P * s_GABA : amp
            ds_GABA / dt = -s_GABA / tau_s_GABA_P : 1

            I_GABAb = j_GABAb_P * s_GABAb: amp
            ds_GABAb / dt = -s_GABAb / tau_s_GABAb_P : 1
            
            I_AMPA_rec = j_AMPA_rec_P * s_AMPA : amp
            ds_AMPA / dt = -s_AMPA / tau_s_AMPA_P: 1    
            
            I_AMPA_spi = j_AMPA_rec_P * s_AMPA_spi : amp
            ds_AMPA_spi / dt = -s_AMPA_spi / tau_s_AMPA_P : 1
            
            ref : second
            v_th : volt
        '''
        
        
    elif type == 'inhibitory':
        eqs = '''
            dv / dt = (-v + V_leak - (I_tot/g_m_I)) / tau_m_I : volt (unless refractory)
        
            dv_ip /dt = (-v_ip -(I_AMPA_rec / g_m_I)) / tau_m_I : volt (unless refractory)
            dv_ii /dt = (-v_ii -(I_GABA_rec / g_m_I)) / tau_m_I : volt (unless refractory)
            
            I_tot = I_AMPA_cor + I_AMPA_tha + I_AMPA_rec + I_GABA_rec + I_injected_I : amp
            
            I_AMPA_cor = j_AMPA_cor_I * s_AMPA_cor : amp
            ds_AMPA_cor / dt = - s_AMPA_cor / tau_s_AMPA_I : 1
            
            I_AMPA_tha = j_AMPA_tha_I * s_AMPA_tha : amp
            ds_AMPA_tha / dt = - s_AMPA_tha / tau_s_AMPA_I_ext : 1
            
            I_GABA_rec = j_GABA_I * s_GABA : amp
            ds_GABA / dt = -s_GABA / tau_s_GABA_I : 1
            
            I_AMPA_rec = j_AMPA_rec_I * s_AMPA : amp
            ds_AMPA / dt = -s_AMPA / tau_s_AMPA_I : 1
            
            ref : second
            v_th : volt
        '''
        
        
    elif type == 'gabab':
        eqs = '''
            dv / dt = (-v + V_leak - (I_tot/g_m_I)) / tau_m_I : volt (unless refractory)

            dv_bi /dt = (-v_bi -(I_GABA / g_m_I)) / tau_m_I : volt (unless refractory)
            dv_bp /dt = (-v_bp -(I_AMPA / g_m_I)) / tau_m_I : volt (unless refractory)
            
            I_tot = I_AMPA + I_GABA + I_AMPA_tha + I_injected_I: amp
            
            I_GABA = agonist * j_GABA_B * s_GABA : amp
            ds_GABA / dt = -s_GABA / tau_s_GABA_I : 1
                        
            I_AMPA = j_AMPA_B * s_AMPA : amp
            ds_AMPA / dt = -s_AMPA / tau_s_AMPA_I : 1
            
            I_AMPA_tha = j_AMPA_tha_I * s_AMPA_tha : amp
            ds_AMPA_tha / dt = - s_AMPA_tha / tau_s_AMPA_I_ext : 1
            
            ref : second
            v_th : volt
        '''
        
        
    else:
        raise ValueError(format('The option type = %s is not a valid one.' %(type)))
    return eqs
