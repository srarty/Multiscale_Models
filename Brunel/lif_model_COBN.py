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
        g_GABA  = 18.5 * nS 
        g_GABAb = 2.3125 * nS 
        g_AMPA  = -36.75 * nS 
        g_AMPA_ext = -112.75 * nS 
        g_AMPA_tha = -112.75 * nS 
         
       
        # Delta function weight (increment with each input spike)
        # Defined experimentally with 'synaptic_functions.py'. Based on the 
        # unitary increment of the single exponential.
        weight = 1
        external_input_weight = 1.15
        
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
        g_GABA  = 22.62 * nS * 0.337967
        g_GABAb  = 0 * nS
        g_AMPA  = -41.25 * nS * 0.012980
        g_AMPA_ext = -1.9 * nS * 0.012980
        g_AMPA_tha = -19 * nS * 0.012980
        
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
        "g_GABA":       g_GABA,
        "g_GABAb":      g_GABAb,
        "g_AMPA":       g_AMPA,
        "g_AMPA_ext":   g_AMPA_ext,
        "g_AMPA_tha":   g_AMPA_tha,
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
            dv_pb /dt = (-v_pb - ( I_GABAb                  / g_m_P)) / tau_m_P : volt (unless refractory)
            dv_pu /dt = (-v_pu - ( I_AMPA_cor               / g_m_P)) / tau_m_P : volt (unless refractory)
        
            I_tot = I_AMPA_cor + I_AMPA_tha + I_AMPA_rec + I_AMPA_spi + I_GABAb + I_GABA_rec + I_injected : amp
            
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
            
            I_AMPA_spi = g_AMPA_rec_P * s_AMPA_spi * (v - v_AMPA): amp
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



#%% Plot 
def plot_results(sp_P, sp_I, Py_monitor, In_monitor, r_P, r_I, lfp_v):
    # spike rates
    window_size = 10*ms#%100.1*ms # Size of the window for the smooth spike rate # 100.1 instead of 100 to avoid an annoying warning at the end of the simulation
    # window_size2 = 0.1*ms
    
    r_P_rate = r_P.smooth_rate(window='gaussian', width=window_size)
    if shape(r_P_rate) != shape(r_P.t):
        r_P_rate = r_P_rate[5:]
    
    r_I_rate = r_I.smooth_rate(window='gaussian', width=window_size)
    if shape(r_I_rate) != shape(r_I.t):
        r_I_rate = r_I_rate[5:]
    
    
    f, axs = plt.subplots(4, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    
    # colors
    c_inter = 'C6'      # pink
    c_py    = 'C9'      # light blue
    c_b    = 'C0'      # blue
    c_Cor   = 'C1'      # orange
    c_gray  = '#e0e0e0' # grey
    
    # raster
    axs[0].set_title('Raster ({} neurons/pop) N = {}, u = {} spikes/ms'.format(N_activity_plot, N, input_spike_rate))
    axs[0].set_ylabel('Neuron')
    axs[0].set_yticks([])
    axs[0].spines["top"].set_visible(False)
    axs[0].spines["right"].set_visible(False)
    
    axs[0].plot(sp_P.t / ms, sp_P.i + N_I + N_B, '.', markersize=2, label='Pyramidal', c=c_py)
    axs[0].plot(sp_I.t / ms, sp_I.i + N_B, '.', markersize=2, label='GABAa', c=c_inter)
    # axs[0].plot(sp_B.t / ms, sp_B.i, '.', markersize=2, label='GABAb', c=c_b)
    axs[0].legend(loc=1)
            
    axs[1].set_title('Population rates, moving average')
    axs[1].set_ylabel('Spike rate (Hz)')
    axs[1].spines["top"].set_visible(False)
    axs[1].spines["right"].set_visible(False)   
        
    axs[1].plot(r_P.t / ms, r_P_rate / Hz, label='Pyramidal', c=c_py)
    axs[1].plot(r_I.t / ms, r_I_rate / Hz, label='GABAa', c=c_inter)
    # axs[1].plot(r_B.t / ms, r_B_rate / Hz, label='GABAb', c=c_b)
    # axs[1].plot(r_Cor.t / ms, r_Cor_rate / Hz, label='Cortico-cortical (OU)', c=c_Cor)
    # axs[1].plot(T_u / ms, u / Hz)
    axs[1].legend(loc=1)
    
    # synaptic currents
    axs[2].set_title('Synaptic currents')
    axs[2].set_ylabel('Amplitude (unitless)')
    axs[2].set_xlabel('Time (ms)')
    axs[2].spines["top"].set_visible(False)
    axs[2].spines["right"].set_visible(False)
    # Others
    # axs[2].plot(T*1000, np.array(st_GABA_I.s_GABA).transpose(), lw=0.5, c=c_gray) # , label='GABA (Inter)'
    # axs[2].plot(T*1000, np.array(st_AMPA_P.s_AMPA).transpose(), lw=0.5, c=c_gray) # , label='AMPA (Py)'
    # axs[2].plot(T*1000, np.array(st_AMPA_cor_I.s_AMPA_cor).transpose(), lw=0.5, c=c_gray) # , label='AMPA_cor (Inter)'
    # # alphas
    # axs[2].plot(T*1000, np.array(st_GABA_P.s_GABA).transpose(), label='GABA (Py)', c=c_py)
    # axs[2].plot(T*1000, np.array(st_AMPA_I.s_AMPA).transpose(), label='AMPA (In)', c=c_inter)
    # axs[2].plot(T*1000, np.array(st_AMPA_cor_P.s_AMPA_cor).transpose(), lw=0.5, c=c_Cor , label='AMPA_cor (Cortical)')
    # axs[2].legend(loc=1)
    
    # LFP
    axs[3].set_title('LFP')
    axs[3].set_xlabel('Time (ms)')
    axs[3].set_ylabel('mV')
    axs[3].spines["top"].set_visible(False)
    axs[3].spines["right"].set_visible(False)
    axs[3].plot(T*1000, np.transpose(lfp_v), lw=0.5, label='LFP_V')
    # axs[3].plot(T*1000, np.transpose(lfp_), lw=0.5, label='LFP_I')
    axs[3].legend(loc=1)

    f.tight_layout() # Fixes the positions of subplots and labels

    # Second figure. PSP
    f2, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    
    axs[0].set_title('Pyramidal Vm (selected cells)')
    axs[0].set_xlabel('Time (ms)')
    axs[0].set_ylabel('mV')
    axs[0].plot(T*1000, np.transpose(Py_monitor.v[0:5])*1e3, lw=0.5, c=c_py)

    axs[1].set_title('Interneurons Vm (selected cells)')
    axs[1].set_xlabel('Time (ms)')
    axs[1].set_ylabel('mV')
    axs[1].plot(T*1000, np.transpose(In_monitor.v[0:5])*1e3, lw=0.5, c=c_inter)
    
    # f3.tight_layout()
   
    # f4, axs = plt.subplots(1, 1, figsize=(6,6))
    # axs.plot(np.transpose(v_pi) * 1000, np.transpose(v_ip) * 1000)
    # axs.set_xlabel('x1 = V_pi')
    # axs.set_ylabel('x3 = V_ip')
    
    plt.show()    
    
    #%% Plot ISI distance and synchronization 
    stp = ast.get_spike_trains(sp_P.t/second, sp_P.i) # Spike trains of Pyramidal neurons
    sti = ast.get_spike_trains(sp_I.t/second, sp_I.i) # Spike trains of inhibitory neurons
    
    # merge spike trains
    merged_p = spk.merge_spike_trains(stp)
    merged_i = spk.merge_spike_trains(sti)
    
    # Define the same interval for both merged spike trains
    if merged_p.t_end > merged_i.t_end:
        merged_p.t_end = merged_i.t_end
    else:
        merged_i.t_end = merged_p.t_end
        
    # Get the distance between Py and In
    d = ast.distance(merged_p, merged_i)
    
    # PSTH
    f = plt.figure()
    ast.plot_psth(stp, fig=f)
    ast.plot_psth(sti, fig=f)
    
    # Spike syncrhony
    ast.spike_si(stp, 'Py')
    ast.spike_si(sti, 'In')
    ast.spike_si([merged_i, merged_p])