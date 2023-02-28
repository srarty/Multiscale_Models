# -*- coding: utf-8 -*-
"""
Created on 8/Dec/2022

@author: Artemio Soto-Breceda [artemios]

Estimate g in COBN from j in CUBN

References: 
    
    [1] CAVALLARI, Stefano et al (2014)
"""

import os
import scipy.io
import numpy as np
from brian2 import *
from scipy import signal
from termcolor import colored  # Coloured text in the terminal
import matplotlib.pyplot as plt
import pyspike as spk
prefs.codegen.target = 'numpy'  # use the Python fallback instead of C compilation
# prefs.codegen.target = 'cython'  # use the Python fallback instead of C compilation
devices.device.shape = []       # This and the following line remove an annoying warning when brian2 is imported and loaded into RAM
devices.device.size = []

import lif_model_CUBN as cubn
import lif_model_COBN as cobn
from lif_plot import plot_results
import analyse_spike_trains as ast

def lif(MODEL='cobn', k_AMPA_P=0.01/mV, k_GABA_P=0.05/mV, k_GABAb_P=0.05/mV, k_AMPA_I=0.01/mV, k_GABA_I=0.1/mV, PLOT = False):    
    #%% Options:
    RECURRENT_PYRAMIDAL = True      # Self excitation 
    RECURRENT_INHIBITORY = True     # Self inhibition
    ACTIVE_INTERNEURONS = True      # Inhibitory population
    ACTIVE_GABAb = True            # Second inhibitory (slow) population (Wendling-like model)
    INHIBIT_INPUT = True           # Excitatory cortical input to inhibitory population
    GAUSSIAN_REFRACTORY = True      # If true, the refractory period of each cell is taken from a gaussian distribution, otherwise it is the same for all
    GAUSSIAN_THRESHOLD = True       # If true, the refractory period of each cell is taken from a gaussian distribution, otherwise it is the same for all
    
    corriente = 0#50
    # Balanced-rate network (?) with input currents: Py = 500.01 pA, In = 398 pA
    input_current = corriente  # 437.5 # 500.01       # Injected current to Pyramidal population # Use this to calculate the nonlinearity (Vm -> Spike_rate sigmoid) on the disconnected model
    input_current_I = corriente # 350 # 398 # 400.01     # Inhibitory interneurons
    
    input_spike_rate = [0]#[0, 2.5, 5, 7.5]#[0, 1, 3, 5] #[u] #[5] #  [0, 2.5, 5] # spikes/ms/cell (driving input)
    input_spike_rate_thalamic = 2 # spikes/ms/cell (spontaneous activity)
    input_spike_rate_thalamic_in = 1.5 # spikes/ms/cell (spontaneous activity)
    
    #%% parameters  --------------------------------------------------------------
    simulation_time = 1 * second
    dt_ = 100 * usecond
    T = np.linspace(0, simulation_time, round(simulation_time/dt_)) # Time vector for plots (in seconds)
       
    # populations
    N = 2000
    N_P = int(N * 4)  # pyramidal neurons
    N_I = int(N * 1)  # interneurons
    
    # set populations parameters
    cubn_py = cubn.set_params('pyramidal')
    cubn_in = cubn.set_params('inhibitory')
    
    cobn_py = cobn.set_params('pyramidal')
    cobn_in = cobn.set_params('inhibitory')
    
    # Probability of connection
    p_IP = cubn_py.get('p_IP')/2
    p_PI = cubn_py.get('p_PI')
    p_PP = cubn_py.get('p_PP')
    p_II = cubn_py.get('p_II')/2
    
    # voltage
    V_leak = -70. * mV      # Resting membrane potential
    V_reset = -59 * mV #    # Reset voltage. Equal to V_leak-> To use Burkitt's, 2006 Eq. (12)
    V_thr = -50 * mV        # Threshold
    v_AMPA = 0 * mV # AMPA reversal potential
    v_GABA = -80 * mV # GABA reversal potential
    v_GABAb = -90 * mV # GABA reversal potential
    
    # membrane capacitance
    C_m_P = cubn_py.get('C')
    C_m_I = cubn_in.get('C')
    
    # membrane leak conductance
    g_m_P = cubn_py.get('g_leak')
    g_m_I = cubn_in.get('g_leak')
    
    # membrane time constants
    tau_m_P = cubn_py.get('tau_m')
    tau_m_I = cubn_in.get('tau_m')
    
    # synaptic time constants
    # Pyramidal
    tau_s_AMPA_P = cubn_py.get('tau_AMPA_s') # Decay time constant (From Brunel and Wang 2001)
    tau_s_AMPA_P_ext = cubn_py.get('tau_AMPA_s_ext') # Decay time constant (From Brunel and Wang 2001)
    tau_s_GABA_P = cubn_py.get('tau_GABA_s')
    tau_s_GABAb_P = cubn_py.get('tau_GABAb_s')
    
    # Inhibitory Interneurons
    tau_s_AMPA_I =  cubn_in.get('tau_AMPA_s')  
    tau_s_AMPA_I_ext = cubn_in.get('tau_AMPA_s_ext')  
    tau_s_GABA_I = cubn_in.get('tau_GABA_s')
    tau_s_GABAb_I = cubn_in.get('tau_GABAb_s')
    
    # refractory period
    tau_rp_P = cubn_py.get('tau_rp')
    tau_rp_I = cubn_in.get('tau_rp')
    
    # Synaptic delay
    delay = 0.0 * ms 
    
    # Cortical input
    num_inputs = 800                    # Both thalamo-cortical and cortico-cortical 
    
    
    
    # Weight constants. Amplitude of the synaptic input
    # Pyramidal 
    increment_AMPA_P        = cubn_py.get('weight')
    increment_AMPA_ext_P    = cubn_py.get('external_input_weight')
    increment_GABA_P        = cubn_py.get('weight')
    
    # Inhibitory interneurons
    increment_AMPA_I        = cubn_in.get('weight')
    increment_AMPA_ext_I    = cubn_in.get('external_input_weight')
    increment_GABA_I        = cubn_in.get('weight')
    
    
    
    
    
    
    # Synaptic efficacies (CUBN)
    # AMPA (excitatory)
    j_AMPA_rec_P = cubn_py.get('j_AMPA') * 2000/N 
    j_AMPA_rec_I = cubn_in.get('j_AMPA') * 2000/N
    
    j_AMPA_cor_P = cubn_py.get('j_AMPA_ext')
    j_AMPA_cor_I = cubn_in.get('j_AMPA_ext')
    
    j_AMPA_tha_P = cubn_py.get('j_AMPA_tha')
    j_AMPA_tha_I = cubn_in.get('j_AMPA_tha')
    
    # GABAergic (inhibitory)
    j_GABA_P    = cubn_py.get('j_GABA') * 2000/N 
    j_GABAb_P   = cubn_py.get('j_GABAb') * 2000/N
    j_GABA_I    = cubn_in.get('j_GABA') * 2000/N 
    
    
    
    
    
    
    # Synaptic conductances (COBN)
    # AMPA (excitatory)
    g_AMPA_rec_P = k_AMPA_P * j_AMPA_rec_P
    g_AMPA_rec_I = k_AMPA_I * j_AMPA_rec_I
    
    g_AMPA_cor_P = k_AMPA_P * j_AMPA_cor_P
    g_AMPA_cor_I = k_AMPA_I * j_AMPA_cor_I
    
    g_AMPA_tha_P = k_AMPA_P * j_AMPA_tha_P
    g_AMPA_tha_I = k_AMPA_I * j_AMPA_tha_I
    
    # GABAergic (inhibitory)
    g_GABA_P    = k_GABA_P * j_GABA_P
    g_GABAb_P   = k_GABAb_P * j_GABAb_P
    g_GABA_I    = k_GABA_I * j_GABA_I
    
    
    
    
    
    
    # Injected current
    I_injected      = -input_current * pA # Input current to Pyramidal population. Sets a baseline spiking rate
    I_injected_I    = -input_current_I * pA # Input current to Pyramidal population. Sets a baseline spiking rate
    
    
    #%% modeling  ----------------------------------------------------------------
    # Model equations
    if MODEL == 'cobn':
        eqs_P = cobn.get_equations('pyramidal')
        eqs_I = cobn.get_equations('inhibitory')
    else:
        eqs_P = cubn.get_equations('pyramidal')
        eqs_I = cubn.get_equations('inhibitory')
    
    # Neuron groups
    Py_Pop = NeuronGroup(N_P, eqs_P, threshold='v > v_th', reset='''v = V_reset
                                                                    v_pe = V_reset-V_leak
                                                                    v_pi = V_reset-V_leak
                                                                    v_pb = V_reset-V_leak
                                                                    v_pp = V_reset-V_leak
                                                                    ''', refractory='ref', method='rk4', dt=dt_, name='PyramidalPop') # Pyramidal population
    Py_Pop.v = V_leak
    
    
    
    In_Pop = NeuronGroup(N_I, eqs_I, threshold='v > V_thr', reset='''v = V_reset
                                                                    v_ip = V_reset-V_leak
                                                                    v_ii = V_reset-V_leak
                                                                    ''', refractory='ref', method='rk4', dt=dt_, name='InhibitoryPop') # Interneuron population
    In_Pop.v = V_leak
    
    # Refractoriness
    if GAUSSIAN_REFRACTORY:
        Py_Pop.ref  = tau_rp_P + (3*ms * randn(N_P,))
        In_Pop.ref  = tau_rp_I + (3*ms * randn(N_I,))
    else:
        Py_Pop.ref  = tau_rp_P
        In_Pop.ref  = tau_rp_I    
        
    # Thresholds
    if GAUSSIAN_THRESHOLD:
        Py_Pop.v_th = V_thr + (3*mV * randn(N_P,))
        In_Pop.v_th = V_thr + (3*mV * randn(N_I,))
    else:
        Py_Pop.v_th = V_thr
        In_Pop.v_th = V_thr
        
    
    #%% synaptic equations -------------------------------------------------------
    eqs_pre_glut_P = '''
    s_AMPA += increment_AMPA_P
    '''
    
    eqs_pre_gaba_P = '''
    s_GABA += increment_GABA_P
    '''
    
    eqs_pre_gabab_P = '''
    s_GABAb += increment_GABA_P
    '''
    
    eqs_pre_glut_I = '''
    s_AMPA += increment_AMPA_I
    '''
    
    eqs_pre_gaba_I = '''
    s_GABA += increment_GABA_I
    '''
       
    # Synapses
    # P to P
    C_P_P = Synapses(Py_Pop, Py_Pop, on_pre=eqs_pre_glut_P, method='rk4', dt=dt_, delay=delay, name='synapses_pp')
    C_P_P.connect('i != j', p = p_PP)
    C_P_P.active = RECURRENT_PYRAMIDAL    # Testing no recursive connections to match NMM
    
    # P to I
    C_P_I = Synapses(Py_Pop, In_Pop, on_pre=eqs_pre_glut_I, method='rk4', dt=dt_, delay=delay, name='synapses_ip')
    C_P_I.connect(p = p_PI)     
    C_P_I.active = ACTIVE_INTERNEURONS
    
    # I to I
    C_I_I = Synapses(In_Pop, In_Pop, on_pre=eqs_pre_gaba_I, method='rk4', dt=dt_, delay=delay, name='synapses_ii')
    C_I_I.connect('i != j', p = p_II)
    C_I_I.active = RECURRENT_INHIBITORY
    
    # I to P
    C_I_P = Synapses(In_Pop, Py_Pop, on_pre=eqs_pre_gaba_P, method='rk4', dt=dt_, delay=delay, name='synapses_pi')
    C_I_P.connect(p = p_IP)    
    C_I_P.active = ACTIVE_INTERNEURONS
    
    # GABAb to P
    C_B_P = Synapses(In_Pop, Py_Pop, on_pre=eqs_pre_gabab_P, method='rk4', dt=dt_, delay=delay, name='synapses_pb')
    C_B_P.connect(p = p_IP)
    C_B_P.active = ACTIVE_INTERNEURONS
    
    
    # External inputs
    # Poisson input (Cortico-cortical)
    input1 =  PoissonInput(Py_Pop, 's_AMPA_cor', num_inputs, (input_spike_rate[0] * 1000/num_inputs) * Hz, increment_AMPA_ext_P)
    if np.size(input_spike_rate) > 1:
        input1.active = False
        input2 =  PoissonInput(Py_Pop, 's_AMPA_cor', num_inputs, (input_spike_rate[1] *1000/num_inputs) * Hz, increment_AMPA_ext_P)
        input2.active = False
    if np.size(input_spike_rate) > 2:
        input3 =  PoissonInput(Py_Pop, 's_AMPA_cor', num_inputs, (input_spike_rate[2] *1000/num_inputs) * Hz, increment_AMPA_ext_P)
        input3.active = False
    if np.size(input_spike_rate) > 3:
        input4 =  PoissonInput(Py_Pop, 's_AMPA_cor', num_inputs, (input_spike_rate[3] *1000/num_inputs) * Hz, increment_AMPA_ext_P)
        input4.active = False
        
    # Poisson input (Thalamic, baseline spike rate)
    C_Tha_P = PoissonInput(Py_Pop, 's_AMPA_tha', num_inputs, (input_spike_rate_thalamic*1000/num_inputs) * Hz, increment_AMPA_ext_P)
    C_Tha_I = PoissonInput(In_Pop, 's_AMPA_tha', num_inputs, (input_spike_rate_thalamic_in*1000/num_inputs) * Hz, increment_AMPA_ext_I)
    
    
    #%% monitors  -----------------------------------------------------------------
    N_activity_plot = 30 # How many neurons in the raster plots (too large takes longer to monitor and plot)
    sp_P = SpikeMonitor(Py_Pop[:]) #[:N_activity_plot])
    sp_I = SpikeMonitor(In_Pop[:]) #[:N_activity_plot])
    # sp_B = SpikeMonitor(B_Pop[:]) #[:N_activity_plot])
    
    
    r_P = PopulationRateMonitor(Py_Pop)
    r_I = PopulationRateMonitor(In_Pop)
    
    st_AMPA_P = StateMonitor(Py_Pop, ('s_AMPA'), record = 0)
    st_GABA_P = StateMonitor(Py_Pop, 's_GABA', record = 0)
    st_AMPA_cor_P = StateMonitor(Py_Pop, 's_AMPA_cor', record = 0)
    
    st_AMPA_I = StateMonitor(In_Pop, 's_AMPA', record = 0)
    st_GABA_I = StateMonitor(In_Pop, 's_GABA', record = 0)
    st_AMPA_cor_I = StateMonitor(In_Pop, 's_AMPA_cor', record = 0)
    
    Py_monitor = StateMonitor(Py_Pop, ['v',  'v_pi', 'I_tot'], record = True) # Monitoring the AMPA and GABA currents in the Pyramidal population
    In_monitor = StateMonitor(In_Pop, ['v', 'v_ip', 'I_tot'], record = True)
    
    
    
    
    #%% simulate  -----------------------------------------------------------------
    net = Network(collect())
    
    input1.active = True
    net.run(simulation_time/size(input_spike_rate), report='stdout') # Run first segment, if running more segments, run for a fraction of simulation_time
    
    if np.size(input_spike_rate) > 1:
        input1.active = False
        input2.active = True
        net.run(simulation_time/size(input_spike_rate), report='stdout') # Run first segment, if running more segments, run for a fraction of simulation_time
    
    if np.size(input_spike_rate) > 2:    
        input2.active = False
        input3.active = True
        net.run(simulation_time/size(input_spike_rate), report='stdout') # Run first segment, if running more segments, run for a fraction of simulation_time
           
    if np.size(input_spike_rate) > 3:
        input3.active = False
        input4.active = True
        net.run(simulation_time/size(input_spike_rate), report='stdout') # Run first segment, if running more segments, run for a fraction of simulation_time
    
        
    #%% analysis ------------------------------------------------------------------
    # Generate LFP
    # current based;
    # lfp = sum((Py_monitor.I_GABAb + Py_monitor.I_GABA_rec),0) - sum((Py_monitor.I_AMPA_cor + Py_monitor.I_AMPA_rec + Py_monitor.I_AMPA_spi),0) # Difference of currents
    # lfp_ = lfp / g_m_P # Sum across all Pyramidal neurons and divide by the leak conductance to get volts
    
    # voltage based:
    mean_v_Py = np.transpose(mean(Py_monitor.v,0) - V_leak) * 1e3
    lfp_v = mean_v_Py/volt 
    
    if PLOT:
        plot_results(sp_P, sp_I, Py_monitor, In_monitor, r_P, r_I, lfp_v)
    
    vPy = mean(mean(Py_monitor.v[: , 2500:], 0))
    vIn = mean(mean(In_monitor.v[: , 2500:], 0))
    
    return vPy, vIn
    


#%% Run Cavallari's approach
v_AMPA = 0 * mV
v_GABA = -80 * mV
v_GABAb = -90 * mV

# First time:
VguessPy, VguessIn = lif(MODEL='cubn')    

#%%
k_AMPA_P = []
k_GABA_P = []
k_GABAb_P = []
k_AMPA_I = []
k_GABA_I = []
Py_diff = []
In_diff = []

k_AMPA_P.append(1/(VguessPy - v_AMPA))
k_GABA_P.append(1/(VguessPy - v_GABA))
k_GABAb_P.append(1/(VguessPy - v_GABAb))
k_AMPA_I.append(1/(VguessIn - v_AMPA))
k_GABA_I.append(1/(VguessIn - v_GABA))

print('k_AMPA_P = %f' %k_AMPA_P[0])
print('k_GABA_P = %f' %k_GABA_P[0])
print('k_GABAb_P = %f' %k_GABAb_P[0])
print('k_AMPA_I = %f' %k_AMPA_I[0])
print('k_GABA_I = %f' %k_GABA_I[0])

# Loop
VsimulPy = 0 * mV
VsimulIn = 0 * mV
iteration = 0

converged = False
while not(converged):
    VsimulPy, VsimulIn = lif(MODEL='cobn', k_AMPA_P=k_AMPA_P[iteration], k_GABA_P=k_GABA_P[iteration], k_GABAb_P=k_GABAb_P[iteration], k_AMPA_I=k_AMPA_I[iteration], k_GABA_I=k_GABA_I[iteration])
    
    Py_diff.append(abs(VsimulPy - VguessPy))
    In_diff.append(abs(VsimulIn - VguessIn))
    print("\niteration = %s | Py_diff = %s | In_diff = %s \nk_AMPA_P=%s | k_GABA_P=%s | k_GABAb_P=%s | k_AMPA_I=%s | k_GABA_I=%s" %(iteration, Py_diff[iteration], In_diff[iteration], k_AMPA_P[iteration], k_GABA_P[iteration], k_GABAb_P[iteration], k_AMPA_I[iteration], k_GABA_I[iteration]))
    
    # converged = (Py_diff <= 0.01*mV) & (In_diff <= 0.01*mV)
    converged = (Py_diff[iteration] <= 0.01*mV) & (In_diff[iteration] <= 0.01*mV)
    
    iteration +=1

    VguessPy = VsimulPy
    VguessIn = VsimulIn

    k_AMPA_P.append(1/(VguessPy - v_AMPA))
    k_GABA_P.append(1/(VguessPy - v_GABA))
    k_GABAb_P.append(1/(VguessPy - v_GABAb))
    k_AMPA_I.append(1/(VguessIn - v_AMPA))
    k_GABA_I.append(1/(VguessIn - v_GABA))
    
    
    
print("Success. k_AMPA_P=%s, k_GABA_P=%s, k_GABAb_P=%s, k_AMPA_I=%s, k_GABA_I=%s" %(k_AMPA_P[-1], k_GABA_P[-1], k_GABAb_P[-1], k_AMPA_I[-1], k_GABA_I[-1]))
    
# If on Spartan, save results:
save_dictionary={'k_AMPA_P': k_AMPA_P,
                 'k_GABA_P': k_GABA_P,
                 'k_GABAb_P': k_GABAb_P,
                 'k_AMPA_I': k_AMPA_I,
                 'k_GABA_I': k_GABA_I,
                 } 

# folder_path = '/data/gpfs/projects/punim0643/artemios/simulations/'
# scipy.io.savemat(folder_path + 'CUBN_to_COBN_results.mat', mdict = save_dictionary)




















