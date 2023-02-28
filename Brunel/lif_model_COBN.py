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
    [3] DESTEXHE, et al (1996) }\ Reversal
    [4] Connors, et al (1988)  }/ potentials
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
    p_IP = 0.437 #* 0.6 #L2/3:0.411 #* 0.25
    p_BP = 0.437 #* 0.4
    p_PI = 0.43  #L2/3:0.395 #* 0.25
    p_PP = 0.243 #L2/3:0.16 #* 0.25
    p_II = 0.451 #L2/3:0.451 #* 0.25
    
    if type == 'pyramidal':
        #%% Pyramidal. Allen
        g_leak = 25 * nS # Leak conductance
        C = 0.5 * nF # Membrane capacitance
        
        tau_rpa = 2 * ms # Absolute refractory period
        tau_rpr = 20 * ms #100 * ms # Relative refractory period
        tau_rp = tau_rpa + tau_rpr # Effective refractory period
        tau_m = 20 * ms # Membrane time constant

        tau_GABA_s = 5.25 * ms
        tau_GABAb_s = 52.5 * ms
        tau_AMPA_s = 2.4 * ms
        tau_AMPA_s_ext = 2.4 * ms
        tau_l = 1 * ms # Latency

        # Synaptic efficacies
        g_GABA  = 18.5 * 0.04371582 * nS
        g_GABAb = 2.3125 * 0.03041824 * nS
        g_AMPA  = -36.75 * -0.01750547 * nS
        g_AMPA_ext = -112.75 * -0.01750547 * nS
        g_AMPA_tha = -112.75 * -0.01750547 * nS
         
       
        # Delta function weight (increment with each input spike)
        # Defined experimentally with 'synaptic_functions.py'. Based on the 
        # unitary increment of the single exponential.
        weight = 1
        external_input_weight = 1
        
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
        g_GABA  = 22.62 * 0.04030332 * nS
        g_GABAb  = 0 * nS
        g_AMPA  = -41.25 * -0.01811983 * nS
        g_AMPA_ext = -1.9  * -0.01811983 * nS
        g_AMPA_tha = -19  * -0.01811983 * nS
        
        # Delta function weight (increment with each input spike)
        weight = 1
        external_input_weight = 1
        
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
        "g_GABA":       g_GABA,
        "g_GABAb":      g_GABAb,
        "g_AMPA":       g_AMPA,
        "g_AMPA_ext":   g_AMPA_ext,
        "g_AMPA_tha":   g_AMPA_tha,
        "p_IP":         p_IP,
        "p_BP":         p_BP,
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
            
            dv_pi /dt = (-v_pi - ( I_GABA_rec / g_m_P)) / tau_m_P : volt (unless refractory)
            dv_pp /dt = (-v_pp - ( I_AMPA_rec / g_m_P)) / tau_m_P : volt (unless refractory)
            dv_pb /dt = (-v_pb - ( I_GABAb    / g_m_P)) / tau_m_P : volt (unless refractory)
            dv_pu /dt = (-v_pu - ( I_AMPA_cor / g_m_P)) / tau_m_P : volt (unless refractory)
        
            I_tot = I_exc + I_inh + I_injected : amp
            I_exc = I_AMPA_cor + I_AMPA_tha + I_AMPA_rec : amp
            I_inh = I_GABAb + I_GABA_rec : amp
            
            I_AMPA_cor = g_AMPA_cor_P * s_AMPA_cor * (v - v_AMPA): amp
            ds_AMPA_cor / dt = -s_AMPA_cor / tau_s_AMPA_P : 1    
            
            I_AMPA_tha = g_AMPA_tha_P * s_AMPA_tha * (v - v_AMPA) : amp
            ds_AMPA_tha / dt = -s_AMPA_tha / tau_s_AMPA_P : 1
            
            I_GABA_rec = g_GABA_P * s_GABA * (v - v_GABA): amp
            ds_GABA / dt = -s_GABA / tau_s_GABA_P : 1

            I_GABAb = g_GABAb_P * s_GABAb * (v - v_GABAb): amp
            ds_GABAb / dt = -s_GABAb / tau_s_GABAb_P : 1
            
            I_AMPA_rec = g_AMPA_rec_P * s_AMPA * (v - v_AMPA): amp
            ds_AMPA / dt = -s_AMPA / tau_s_AMPA_P: 1    
                        
            ref : second
            v_th : volt
        '''
        
        
    elif type == 'inhibitory':
        eqs = '''
            dv / dt = (-v + V_leak - (I_tot/g_m_I)) / tau_m_I : volt (unless refractory)
        
            dv_ip /dt = (-v_ip -(I_AMPA_rec / g_m_I)) / tau_m_I : volt (unless refractory)
            dv_ii /dt = (-v_ii -(I_GABA_rec / g_m_I)) / tau_m_I : volt (unless refractory)
           
            I_tot = I_exc + I_inh + I_injected_I: amp
            I_exc = I_AMPA_cor + I_AMPA_tha + I_AMPA_rec : amp
            I_inh = I_GABA_rec : amp
            
            I_AMPA_cor = g_AMPA_cor_I * s_AMPA_cor * (v - v_AMPA) : amp
            ds_AMPA_cor / dt = - s_AMPA_cor / tau_s_AMPA_I : 1
            
            I_AMPA_tha = g_AMPA_tha_I * s_AMPA_tha * (v - v_AMPA) : amp
            ds_AMPA_tha / dt = - s_AMPA_tha / tau_s_AMPA_I_ext : 1
            
            I_GABA_rec = g_GABA_I * s_GABA * (v - v_GABA) : amp
            ds_GABA / dt = -s_GABA / tau_s_GABA_I : 1
            
            I_AMPA_rec = g_AMPA_rec_I * s_AMPA * (v - v_AMPA): amp
            ds_AMPA / dt = -s_AMPA / tau_s_AMPA_I : 1
            
            ref : second
            v_th : volt
        '''
        
    else:
        raise ValueError(format('The option type = %s is not a valid one.' %(type)))
    return eqs
