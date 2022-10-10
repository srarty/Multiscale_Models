# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 12:01:40 2022

@author: Artemio Soto-Breceda [artemios]
"""

import os
import scipy.io
import numpy as np
from brian2 import *
from scipy import signal
from termcolor import colored  # Coloured text in the terminal
import matplotlib.pyplot as plt
prefs.codegen.target = 'numpy'  # use the Python fallback instead of C compilation
devices.device.shape = []       # This and the following line remove an annoying warning when brian2 is imported and loaded into RAM
devices.device.size = []

from lif_model import set_params

# # Save commands:
# # Pyramidal:
# save = {'ipsp': Py_monitor.v1*1e3}
# scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/pyramidal_ipsp.mat', mdict=save)

# # Interneurons:
# save = {'epsp': Py_monitor.v1*1e3}
# scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/inhibitory_epsp.mat', mdict=save)

# RECURSIVE Pyramidal:
# save = {'epsp': Py_monitor.v1*1e3}
# scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/pyramidal_epsp.mat', mdict=save)

# # RECURSIVE Interneurons:
# save = {'ipsp': Py_monitor.v1*1e3}
# scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/inhibitory_ipsp.mat', mdict=save)

# # External input -> Pyramidal:
# save = {'epsp': Py_monitor.v1*1e3}
# scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/pyramidal_externalEPSP.mat', mdict=save)

# # External input -> Inhibitory:
# save = {'epsp': Py_monitor.v1*1e3}
# scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/inhibitory_externalEPSP.mat', mdict=save)

def synaptic_functions_exploration(alpha_ei='',alpha_ie='',alpha_ee='',alpha_ii='', PLOT=True):
    # plt.close('all')
    #%% options  --------------------------------------------------------------
    
    source          = 'three_pop'       # 'brunel', 'allen'  or 'three_pop'
    synaptic_type   = 'GABA'        # AMPA (excitatory), GABA (inhibitory) or GABAb
    neuron_type     = 'pyramidal'  # pyramidal, inhibitory or gabab
    external        = False         # When AMPA, synapsis can be external or recurrent (local)
    input_spike_rate = 0            # spikes/ms/cell 
    simulation_time = 0.3 * second
    
        
    #%% parameters  -----------------------------------------------------------
    dt_ = 100 * usecond
    T = linspace(0, simulation_time, round(simulation_time/dt_)) # Time vector for plots (in seconds)
    
    # populations
    N_P = 1
    
    # voltage
    V_leak = -70. * mV      # Resting membrane potential
    V_thr = -50 * mV        # Threshold
    V_reset = -59 * mV      # Reset voltage. Equal to V_leak-> To use Burkitt's, 2006 Eq. (12)
    
    # parse inputs to synaptic_functions_exploration -------------------------
    if alpha_ee != '':
        neuron_type = 'pyramidal'
        synaptic_type = 'AMPA'
        params = set_params(neuron_type, source)
        params["j_AMPA"] = -1 * alpha_ee * pA
        
    elif alpha_ii != '':        
        neuron_type = 'inhibitory'
        synaptic_type = 'GABA'
        params = set_params(neuron_type, source)
        params["j_GABA"] = -1 * alpha_ii * pA
        
    elif alpha_ei != '':        
        neuron_type = 'inhibitory'
        synaptic_type = 'AMPA'
        params = set_params(neuron_type, source)
        params["j_AMPA"] = -1 * alpha_ei * pA
        
    elif alpha_ie != '':
        neuron_type = 'pyramidal'
        synaptic_type = 'GABA'
        params = set_params(neuron_type, source)
        params["j_GABA"] = -1 * alpha_ie * pA
    
    else:
        params = set_params(neuron_type, source)
    
    # membrane capacitance
    C_m = params["C"]
    
    # membrane leak
    g_m = params["g_leak"]
    
    # membrane time constants
    tau_m = C_m / g_m
    
    # connectivity time constants
    if synaptic_type == 'AMPA':
        tau_s =  params["tau_AMPA_s"]
    else:
        tau_s =  params["tau_GABA_s"]             
        
    
    tau_l = 0.0 * ms
    tau_rp =  params["tau_rp"] # refractory period
    
    # Cortical input
    num_inputs = 1                    # Both thalamo-cortical and cortico-cortical 
    I_input = -300 * pA
    
    # Synaptic efficacies
    # AMPA (external)
    if synaptic_type == 'AMPA':
        if external:
            j =  params["j_AMPA_ext"]
            alpha_weight = params["external_input_weight"]
        else:
            j =  params["j_AMPA"]        
            alpha_weight = params["weight"]
    else:
        j =  params["j_GABA"]
        alpha_weight = params["weight"]
    
    
    # Alpha function's parameter (and double exponential) to fix the units in ds/dt
    k = 1 / ms # Dimmensionless?, check Nicola and Campbell 2013
    
    #%%
    # model equations
    # Input
    input_eqs = '''
    dv / dt = (-v + V_leak - (I_input/g_m)) / tau_m : volt (unless refractory)
    '''
    
    # Population
    eqs = '''
    dv / dt = (-v + V_leak - (I_tot/g_m)) / tau_m : volt (unless refractory)
    
    dv1 / dt = (-v1 - (I_AMPA1/g_m)) / tau_m : volt (unless refractory)
    dv6 / dt = (-v6 - (I_AMPA6/g_m)) / tau_m : volt (unless refractory)
    
    I_tot = I_AMPA1 + I_AMPA6: amp
    
    I_AMPA1 = j * s_AMPA1 : amp
    ds_AMPA1 / dt = -s_AMPA1 / (tau_s) : 1    
    
    I_AMPA6 = (-896 * pA) * s_AMPA6 : amp
    s_AMPA6 : 1
    '''
    
    eqs_pre_ampa1 = '''
    s_AMPA1 += alpha_weight
    ''' 
    
    eqs_ampa6 ='''
    s_AMPA6_post = s_AMPA6_syn : 1 (summed)
    ds_AMPA6_syn / dt = - s_AMPA6_syn / (tau_s) + k * x6 : 1 (clock-driven)
    dx6 / dt = - x6 / (tau_s) : 1 (clock-driven)
    '''
    
    eqs_pre_ampa6 = '''
    x6 += 1
    ''' 
    
    Pyramidal = NeuronGroup(N_P, eqs, threshold='v > V_thr', reset='''v = V_reset
                                                    v1 = V_reset-V_leak
                                                    v6 = V_reset-V_leak
                                                    ''',refractory=tau_rp, method='rk4', dt=dt_, name='PyramidalPop') # Pyramidal population
    Pyramidal.v = V_leak
    
  
    if num_inputs > 1:
        Input = PoissonGroup(num_inputs, rates=input_spike_rate * Hz, dt=dt_)
    else:
        Input = NeuronGroup(num_inputs, input_eqs, threshold='v > V_thr', reset='v = V_reset', refractory=tau_rp, method='rk4', dt=dt_) # Pyramidal population
    
    # if external:
    #     AMPA1_synapses = Synapses(Input, Pyramidal, on_pre=eqs_pre_ampa1, method='rk4', clock=Input.clock)
    #     AMPA1_synapses.connect(p = 1)
    
    AMPA1_synapses = Synapses(Input, Pyramidal, on_pre=eqs_pre_ampa1, method='rk4', delay=tau_l)#, dt=10*usecond)
    AMPA1_synapses.connect(p = 1)
        
    AMPA6_synapses = Synapses(Input, Pyramidal, model=eqs_ampa6, on_pre=eqs_pre_ampa6, method='rk4')#, dt=10*usecond)
    AMPA6_synapses.connect(p = 1)
    AMPA6_synapses.active = False
    
    Py_monitor = StateMonitor(Pyramidal, ['v', 'v1', 'v6', 'I_AMPA1', 'I_AMPA6', 'I_tot'], record = True) # Monitoring the AMPA and GABA currents in the Pyramidal population

    #%% Run
    net = Network(collect())
    net.run(simulation_time, report='stdout')
    
    
    #%% Plot
    if PLOT:
        f, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
            
        axs[0].set_title('Neuron type: {} | Synapses: {} | j: {} pA'.format(neuron_type, synaptic_type, j/pA))
        axs[0].set_ylabel('PSP (mV)')
        axs[0].plot(T * 1e3, (np.transpose(Py_monitor.v1) * 1e3), lw=1, label='v1 (single exp)')
        # axs[0].plot(T * 1e3, (np.transpose(Py_monitor.v6) * 1e3), lw=1, label='v6 (alpha)', linestyle='dashed')
        # axs[0].plot(T * 1e3, ((np.transpose(Py_monitor.v) - V_leak) * 1e3), lw=1, label='Vm', linestyle='solid')
        # axs[0].plot(T * 1e3, (np.transpose(Py_monitor.v1 + Py_monitor.v6) * 1e3), lw=1, label='Sum', linestyle='dashed')
        # if external:
            # axs[0].plot(T * 1e3, (np.transpose(Py_monitor.v1) * 1e3), lw=1, label='single exp')
        axs[0].legend()
            
        axs[1].set_xlabel('Time (ms)')
        axs[1].set_ylabel('PSC (pA)')
        axs[1].plot(T * 1e3, np.transpose(Py_monitor.I_AMPA1), lw=1)
        # axs[1].plot(T * 1e3, np.transpose(Py_monitor.I_AMPA6), lw=1, linestyle='dashed')
        # axs[1].plot(T * 1e3, np.transpose(Py_monitor.I_tot), lw=1)
        # if external:
            # axs[1].plot(T * 1e3, (np.transpose(Py_monitor.s_AMPA1)*j/pA), lw=1)
        
        # f.tight_layout()
        
        plt.show()
    
    return Py_monitor

Py_monitor = synaptic_functions_exploration()