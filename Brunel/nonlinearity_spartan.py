# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 13:20:17 2022

@author: Artemio Soto-Breceda [artemios]

Current based
3 populations (like JR, Spiny Stellate cells optional)

References: 
    
    [1] BRUNEL, Nicolas et WANG, Xiao-Jing. Effects of neuromodulation in a 
        cortical network model of object working memory dominated by recurrent 
        inhibition. Journal of computational neuroscience, 2001, vol. 11, no 1,
        p. 63-85.
        
    [2] MARTINEZ-CANADA, Pablo et al. Computation of the electroencephalogram 
        (EEG) from network models of point neurons. PLoS Computational Biology,
        2021, vol. 17, no 4.
    
    [3] MAIMON (2009) Experimental firing rates. Mean firing rate = 16.6 Hz,max firing rate ~= 140 spikes/s (monkey, visual and motor cortex)
    
    [4] ZIBURKUS (2006) Same as above. Different firing rates for Inhibitory and Pyramidal neurons
    
    [5] CAVALLARI, Stefano et al (2014)
"""

import os
import scipy.io
import numpy as np
from brian2 import *
from scipy import signal
#from termcolor import colored  # Coloured text in the terminal
import matplotlib.pyplot as plt
prefs.codegen.target = 'cython'  # use the Python fallback instead of C compilation
devices.device.shape = []       # This and the following line remove an annoying warning when brian2 is imported and loaded into RAM
devices.device.size = []
from brian2.__init__ import clear_cache

from lif_model import set_params, get_equations

#def brunel(p_II = 0.125):
#def brunel(p_PP = 0.16, p_II = 0.125):
#def brunel(p_PP = 0.16, p_II = 0.451, u = 20):
#def brunel(alpha_ii=0, u=1):
def brunel(corriente = 0):
    
    # arguments from sbatch
    # u = float(sys.argv[1])
    	    
        
    
    #%% Options:
    RECURRENT_PYRAMIDAL = False    # Self excitation 
    RECURRENT_INHIBITORY = False   # Self inhibition
    INHIBIT_INPUT = False         # Excitatory cortical input to inhibitory population
    PARAMS_SOURCE = 'three_pop'       # 'brunel' or 'allen' or 'three_pop'
    ACTIVE_INTERNEURONS = True    # Inhibitory population
    ACTIVE_GABAb = True           # Second inhibitory (slow) population (Wendling-like model)
    GAUSSIAN_REFRACTORY = True    # If true, the refractory period of each cell is taken from a gaussian distribution, otherwise it is the same for all
    GAUSSIAN_THRESHOLD = True     # If true, the refractory period of each cell is taken from a gaussian distribution, otherwise it is the same for all
    SAVE = False                   # Save ground truth data
    PLOT = True                   # Plot results (main Figure)
    PLOT_EXTRA = True             # Plot extra things.
    PSP_FR = 0                    # Presynaptic firing rate for TEST_PSP (TEST_PSP needs to be diff to none for this to take effect)                               
    TEST_PSP = 'none'             # Testing the post synaptic potential of given synapses to a specified input firing rate. Options: 'pu', 'pp', 'pi', 'ii', 'ip', 'bp', 'bi', 'pb', 'none'. To prevent neurons spiking, make V_thr large.
    GABA_A_MULTIPLIER = 1         # GABA_A Agonist applied for the whole duration
    MIDWAY_MULTIPLIER = 1         # GABA_A Agonist applied midsimulation
    RUNTYPE = 'normal'            # Simulation: 'normal', 'current_pulse' or 'gaba_agonist'
    
    corriente = 0
    input_current = corriente # Injected current to Pyramidal population # Use this to calculate the nonlinearity (Vm -> Spike_rate sigmoid) on the disconnected model
    input_current_I = corriente # Inhibitory interneurons
    
    input_spike_rate = [0, 0.5, 1]#[0, 2.5, 5, 7.5]#[0, 1, 3, 5] #[u] #[5] #  [0, 2.5, 5] # spikes/ms/cell (driving input)
    input_spike_rate_thalamic = 1.5 # spikes/ms/cell (spontaneous activity)
    
    #%% parameters  --------------------------------------------------------------
    simulation_time = 1 * second
    dt_ = 100 * usecond
    T = np.linspace(0, simulation_time, round(simulation_time/dt_)) # Time vector for plots (in seconds)
       
    # populations
    N = 2000
    N_P = int(N * 4)  # pyramidal neurons
    N_I = int(N * (0.5 ** ACTIVE_GABAb))    # interneurons
    N_B = int(N * 0.5)  # GABAb
    
    # set populations parameters
    params_py = set_params('pyramidal', PARAMS_SOURCE)
    params_in = set_params('inhibitory', PARAMS_SOURCE)
    params_b = set_params('gabab', PARAMS_SOURCE)
    
    # Probability of connection
    p_IP = params_py.get('p_IP') # * np.sqrt(1000/N) #0.2 #* 100/N # Inhibitory to Pyramidal
    p_PI = params_py.get('p_PI') # * np.sqrt(1000/N) #0.2 #* 100/N # Pyramidal to Inhibitory
    p_PP = params_py.get('p_PP') # * np.sqrt(1000/N)  #0.2 #* 100/N # recurrent excitation (pyramidal) # Generally less than PI, IP connectivity (Bryson et al., 2021)
    p_II = params_py.get('p_II') # * np.sqrt(1000/N)  #0.2 #* 100/N # recurrent inhibition
    
    # voltage
    V_leak = -70. * mV      # Resting membrane potential
    V_reset = -59 * mV #    # Reset voltage. Equal to V_leak-> To use Burkitt's, 2006 Eq. (12)
    if TEST_PSP == 'none':
        V_thr = -50 * mV        # Threshold
    else:
        V_thr = 1000 * mV        # Threshold
    
    # membrane capacitance
    C_m_P = params_py.get('C')
    C_m_I = params_in.get('C')
    
    # membrane leak conductance
    g_m_P = params_py.get('g_leak')
    g_m_I = params_in.get('g_leak')
    
    # membrane time constants
    tau_m_P = params_py.get('tau_m')
    tau_m_I = params_in.get('tau_m')
    
    # synaptic time constants
    # Pyramidal
    tau_s_AMPA_P = params_py.get('tau_AMPA_s') # Decay time constant (From Brunel and Wang 2001)
    tau_s_AMPA_P_ext = params_py.get('tau_AMPA_s_ext') # Decay time constant (From Brunel and Wang 2001)
    tau_s_GABA_P = params_py.get('tau_GABA_s')
    tau_s_GABAb_P = params_py.get('tau_GABAb_s')
    
    # Inhibitory Interneurons
    tau_s_AMPA_I =  params_in.get('tau_AMPA_s')  
    tau_s_AMPA_I_ext = params_in.get('tau_AMPA_s_ext')  
    tau_s_GABA_I = params_in.get('tau_GABA_s')
    
    # refractory period
    tau_rp_P = params_py.get('tau_rp')
    tau_rp_I = params_in.get('tau_rp')
    
    # Synaptic delay
    delay = 0.0 * ms 
    
    # Cortical input
    num_inputs = 800                    # Both thalamo-cortical and cortico-cortical 
    
    
    # Synaptic efficacies
    # AMPA (excitatory)
    j_AMPA_rec_P = params_py.get('j_AMPA') * 2000/N # * np.sqrt(1000/N)
    j_AMPA_rec_I = params_in.get('j_AMPA') * 2000/N # * np.sqrt(1000/N)
    j_AMPA_B = params_b.get('j_AMPA') * 2000/N # * np.sqrt(1000/N)
        
    j_AMPA_cor_P = params_py.get('j_AMPA_ext')
    j_AMPA_cor_I = params_in.get('j_AMPA_ext')
    
    j_AMPA_tha_P = params_py.get('j_AMPA_tha')
    j_AMPA_tha_I = params_in.get('j_AMPA_tha')
    
    # GABAergic (inhibitory)
    j_GABA_P = GABA_A_MULTIPLIER * params_py.get('j_GABA') * 2000/N # * np.sqrt(1000/N)
    j_GABAb_P = params_py.get('j_GABAb') * 2000/N # * np.sqrt(1000/N)
    j_GABA_I = GABA_A_MULTIPLIER * params_in.get('j_GABA') * 2000/N # * np.sqrt(1000/N)
    j_GABA_B = GABA_A_MULTIPLIER * params_b.get('j_GABA') * 2000/N # * np.sqrt(1000/N)
    
    # Weight constants. Amplitude of the synaptic input
    # Pyramidal 
    increment_AMPA_ext_P = params_py.get('external_input_weight')
    increment_AMPA_ext_I = params_in.get('external_input_weight')
    
    
    # Injected current
    I_injected = -input_current * pA # Input current to Pyramidal population. Sets a baseline spiking rate
    I_injected_I = -input_current_I * pA # Input current to Pyramidal population. Sets a baseline spiking rate
    
    #%% modeling  ----------------------------------------------------------------
    # model equations
    eqs_P = get_equations('pyramidal')
    eqs_I = get_equations('inhibitory')
    eqs_B = get_equations('gabab')
    
    
    # Neuron groups
    Py_Pop = NeuronGroup(N_P, eqs_P, threshold='v > v_th', reset='''v = V_reset
                                                                    v_pu = V_reset-V_leak
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
    
    B_Pop = NeuronGroup(N_I, eqs_B, threshold='v > V_thr', reset='''v = V_reset
                                                                    v_bi = V_reset-V_leak
                                                                    v_bp = V_reset-V_leak
                                                                    ''', refractory='ref', method='rk4', dt=dt_, name='GABAbPop') # Interneuron population
    B_Pop.v = V_leak
    
    # Refractoriness
    if GAUSSIAN_REFRACTORY:
        Py_Pop.ref = tau_rp_P + (3*ms * randn(N_P,))
        In_Pop.ref = tau_rp_I + (3*ms * randn(N_I,))
    else:
        Py_Pop.ref = tau_rp_P
        In_Pop.ref = tau_rp_I
        
        
    # Thresholds
    if GAUSSIAN_THRESHOLD:
        Py_Pop.v_th = V_thr + (3*mV * randn(N_P,))
        In_Pop.v_th = V_thr + (3*mV * randn(N_I,))
    else:
        Py_Pop.v_th = V_thr
        In_Pop.v_th = V_thr
        
        
    
    # Custom Poisson population to test PSP on a given synapse
    if (TEST_PSP == 'ip') | (TEST_PSP == 'pp') | (TEST_PSP == 'bp'):
        N_PSP_Test = N_P
    elif (TEST_PSP == 'pi') | (TEST_PSP == 'ii') | (TEST_PSP == 'bi'):
        N_PSP_Test = N_I
    elif (TEST_PSP == 'ib') | (TEST_PSP == 'pb'):
        N_PSP_Test = N_B
    else:
        N_PSP_Test = 0    
    Pop_PSP_Test = PoissonGroup(N_PSP_Test, rates = PSP_FR * Hz, dt=dt_) # poisson input
    
    #%% synaptic equations -------------------------------------------------------
    eqs_pre_glut = '''
    s_AMPA += 1
    '''
    
    eqs_pre_gaba = '''
    s_GABA += 1
    '''
    
    # GABAb
    eqs_pre_gabab = '''
    s_GABAb += 1
    '''
    
       
    # Synapses
    # P to P
    C_P_P = Synapses(Py_Pop, Py_Pop, on_pre=eqs_pre_glut, method='rk4', dt=dt_, delay=delay, name='synapses_pp')
    C_P_P.connect('i != j', p = p_PP)
    C_P_P.active = RECURRENT_PYRAMIDAL    # Testing no recursive connections to match NMM
    
    # P to I
    C_P_I = Synapses(Py_Pop, In_Pop, on_pre=eqs_pre_glut, method='rk4', dt=dt_, delay=delay, name='synapses_ip')
    C_P_I.connect(p = p_PI)     
    C_P_I.active = ACTIVE_INTERNEURONS
    
    # I to I
    C_I_I = Synapses(In_Pop, In_Pop, on_pre=eqs_pre_gaba, method='rk4', dt=dt_, delay=delay, name='synapses_ii')
    C_I_I.connect('i != j', p = p_II)
    C_I_I.active = RECURRENT_INHIBITORY & (not ACTIVE_GABAb)
    
    # I to P
    C_I_P = Synapses(In_Pop, Py_Pop, on_pre=eqs_pre_gaba, method='rk4', dt=dt_, delay=delay, name='synapses_pi')
    C_I_P.connect(p = p_IP)    
    C_I_P.active = ACTIVE_INTERNEURONS
    
    # P to GABAb
    C_P_B = Synapses(Py_Pop, B_Pop, on_pre=eqs_pre_glut, method='rk4', dt=dt_, delay=delay, name='synapses_bp')
    C_P_B.connect(p = p_PI)     
    C_P_B.active = ACTIVE_GABAb
    
    # I to GABAb
    C_I_B = Synapses(B_Pop, In_Pop, on_pre=eqs_pre_gaba, method='rk4', dt=dt_, delay=delay, name='synapses_bi')
    C_I_B.connect(p = p_II)
    C_I_B.active = ACTIVE_GABAb
    
    # GABAb to P
    C_B_P = Synapses(B_Pop, Py_Pop, on_pre=eqs_pre_gabab, method='rk4', dt=dt_, delay=delay, name='synapses_pb')
    C_B_P.connect(p = p_IP)
    C_B_P.active = ACTIVE_GABAb
    
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
    C_Tha_I = PoissonInput(In_Pop, 's_AMPA_tha', num_inputs, (input_spike_rate_thalamic*1000/num_inputs) * Hz, increment_AMPA_ext_I)
    
    # Testing PSP on chosen synapses
    ACTIVE_TEST = True
    if (TEST_PSP == 'pi'):
        C_Test_PSP = Synapses(Pop_PSP_Test, Py_Pop, on_pre=eqs_pre_gaba, method='rk4', dt=dt_, delay=delay, name='synapses_pi_test')
        connection_probability = p_IP
    elif (TEST_PSP == 'pp'):
        C_Test_PSP = Synapses(Pop_PSP_Test, Py_Pop, on_pre=eqs_pre_glut, method='rk4', dt=dt_, delay=delay, name='synapses_pp_test')
        connection_probability = p_PP
    elif (TEST_PSP == 'ip'):
        C_Test_PSP = Synapses(Pop_PSP_Test, In_Pop, on_pre=eqs_pre_glut, method='rk4', dt=dt_, delay=delay, name='synapses_ip_test')
        connection_probability = p_PI
    elif (TEST_PSP == 'ii'):
        C_Test_PSP = Synapses(Pop_PSP_Test, In_Pop, on_pre=eqs_pre_gaba, method='rk4', dt=dt_, delay=delay, name='synapses_ii_test')
        connection_probability = p_II
    elif (TEST_PSP == 'bi'):
        C_Test_PSP = Synapses(Pop_PSP_Test, B_Pop, on_pre=eqs_pre_gaba, method='rk4', dt=dt_, delay=delay, name='synapses_bi_test')
        connection_probability = p_II
    elif (TEST_PSP == 'bp'):
        C_Test_PSP = Synapses(Pop_PSP_Test, B_Pop, on_pre=eqs_pre_glut, method='rk4', dt=dt_, delay=delay, name='synapses_bp_test')
        connection_probability = p_PI
    elif (TEST_PSP == 'pb'):
        C_Test_PSP = Synapses(Pop_PSP_Test, Py_Pop, on_pre=eqs_pre_gabab, method='rk4', dt=dt_, delay=delay, name='synapses_pb_test')
        connection_probability = p_IP
    else:
        C_Test_PSP = Synapses(Pop_PSP_Test, Pop_PSP_Test, method='rk4', name='none_synapse') # This won't do anything. It's to prevent runtime errors
        connection_probability = 1
        ACTIVE_TEST = False
    
    C_Test_PSP.connect('i != j', p = connection_probability)
    C_Test_PSP.active = ACTIVE_TEST 
    
    
    
    #%% monitors  -----------------------------------------------------------------
    N_activity_plot = 30 # How many neurons in the raster plots (too large takes longer to monitor and plot)
    sp_P = SpikeMonitor(Py_Pop[:]) #[:N_activity_plot])
    sp_I = SpikeMonitor(In_Pop[:]) #[:N_activity_plot])
    sp_B = SpikeMonitor(B_Pop[:]) #[:N_activity_plot])
    
    
    r_P = PopulationRateMonitor(Py_Pop)#[0:N_activity_plot])
    r_I = PopulationRateMonitor(In_Pop)
    r_B = PopulationRateMonitor(B_Pop)
    
    st_AMPA_P = StateMonitor(Py_Pop, ('s_AMPA'), record = 0)
    st_GABA_P = StateMonitor(Py_Pop, 's_GABA', record = 0)
    st_AMPA_cor_P = StateMonitor(Py_Pop, 's_AMPA_cor', record = 0)
    
    st_AMPA_I = StateMonitor(In_Pop, 's_AMPA', record = 0)
    st_GABA_I = StateMonitor(In_Pop, 's_GABA', record = 0)
    
    Py_monitor = StateMonitor(Py_Pop, ['I_AMPA_cor', 'I_AMPA_rec', 'I_GABA_rec', 'I_GABAb', 'v',  'v_pi', 'I_tot'], record = True) # Monitoring the AMPA and GABA currents in the Pyramidal population
    In_monitor = StateMonitor(In_Pop, ['v', 'v_ip', 'I_tot', 'v_ii'], record = True)
    B_monitor = StateMonitor(B_Pop, ['v', 'I_tot'], record = True)
    
    
    
    
    #%% simulate  -----------------------------------------------------------------
    net = Network(collect())
    agonist = 1 # Default simulation should run with agonist = 1 so the GABA_A synapse is not boosted.
    
    if RUNTYPE == 'normal':
        ## NORMAL:
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
    
    elif RUNTYPE == 'current_pulse':
        # INJECTED CURRENT PULSE
        I_injected = 0 * pA
        I_injected_I = 0 * pA
        net.run(1.49 * second, report='stdout') # Run first segment, if running more segments, run for a fraction of simulation_time
        
        I_injected = -input_current * pA # Input current to Pyramidal population. Sets a baseline spiking rate
        I_injected_I = -input_current_I * pA # Input current to Pyramidal population. Sets a baseline spiking rate
        net.run(0.01 * second, report='stdout') # Run first segment, if running more segments, run for a fraction of simulation_time
        
        I_injected = 0 * pA
        I_injected_I = 0 * pA
        net.run(0.5 * second, report='stdout') # Run first segment, if running more segments, run for a fraction of simulation_time
    
    elif RUNTYPE == 'gaba_agonist':
        ## GABA_A AGONIST
        
        net.run(simulation_time/2, report='stdout')
        agonist = MIDWAY_MULTIPLIER
        net.run(simulation_time/2, report='stdout')
                
    else:
        print(colored('The simulation did not run. RUNTYPE most be one of the three options.', 'yellow'))
        
    #%% analysis ------------------------------------------------------------------
    # spike rates
    window_size = 10*ms#%100.1*ms # Size of the window for the smooth spike rate # 100.1 instead of 100 to avoid an annoying warning at the end of the simulation
    # window_size2 = 0.1*ms
    
    r_P_rate = r_P.smooth_rate(window='gaussian', width=window_size)
    if shape(r_P_rate) != shape(r_P.t):
        r_P_rate = r_P_rate[5:]
    
    r_I_rate = r_I.smooth_rate(window='gaussian', width=window_size)
    if shape(r_I_rate) != shape(r_I.t):
        r_I_rate = r_I_rate[5:]
    
    
    r_B_rate = r_B.smooth_rate(window='gaussian', width=window_size)
    if shape(r_B_rate) != shape(r_B.t):
        r_B_rate = r_B_rate[5:]
    
    # r_P_rate2 = r_P.smooth_rate(window='flat', width=window_size2)
    # if shape(r_P_rate2) != shape(r_P.t):
    #     r_P_rate2 = r_P_rate2[5:]
    
    # r_I_rate2 = r_I.smooth_rate(window='flat', width=window_size2)
    # if shape(r_I_rate2) != shape(r_I.t):
    #     r_I_rate2 = r_I_rate2[5:]
        
    # Calculate mean PSP (NMM states)
    v_pi = mean(Py_monitor.v_pi, 0)
    v_ip = mean(In_monitor.v_ip, 0)
    
    # Generate LFP
    # current based;
    lfp = sum((Py_monitor.I_GABAb + Py_monitor.I_GABA_rec),0) - sum((Py_monitor.I_AMPA_cor + Py_monitor.I_AMPA_rec),0) # Difference of currents
    lfp_ = lfp / g_m_P # Sum across all Pyramidal neurons and divide by the leak conductance to get volts
    # voltage based:
    mean_v_Py = np.transpose(mean(Py_monitor.v,0) - V_leak) * 1e3
    lfp_v = mean_v_Py/volt 
           
    #%% plotting  -----------------------------------------------------------------
    if PLOT:
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
        axs[0].plot(sp_B.t / ms, sp_B.i, '.', markersize=2, label='GABAb', c=c_b)
        axs[0].legend(loc=1)
                
        axs[1].set_title('Population rates, moving average')
        axs[1].set_ylabel('Spike rate (Hz)')
        axs[1].spines["top"].set_visible(False)
        axs[1].spines["right"].set_visible(False)   
            
        axs[1].plot(r_P.t / ms, r_P_rate / Hz, label='Pyramidal', c=c_py)
        axs[1].plot(r_I.t / ms, r_I_rate / Hz, label='GABAa', c=c_inter)
        axs[1].plot(r_B.t / ms, r_B_rate / Hz, label='GABAb', c=c_b)
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
        axs[2].plot(T*1000, np.array(st_GABA_I.s_GABA).transpose(), lw=0.5, c=c_gray) # , label='GABA (Inter)'
        axs[2].plot(T*1000, np.array(st_AMPA_P.s_AMPA).transpose(), lw=0.5, c=c_gray) # , label='AMPA (Py)'
        axs[2].plot(T*1000, np.array(st_AMPA_cor_I.s_AMPA_cor).transpose(), lw=0.5, c=c_gray) # , label='AMPA_cor (Inter)'
        # alphas
        axs[2].plot(T*1000, np.array(st_GABA_P.s_GABA).transpose(), label='GABA (Py)', c=c_py)
        axs[2].plot(T*1000, np.array(st_AMPA_I.s_AMPA).transpose(), label='AMPA (In)', c=c_inter)
        axs[2].plot(T*1000, np.array(st_AMPA_cor_P.s_AMPA_cor).transpose(), lw=0.5, c=c_Cor , label='AMPA_cor (Cortical)')
        axs[2].legend(loc=1)
        
        # LFP
        axs[3].set_title('LFP')
        axs[3].set_xlabel('Time (ms)')
        axs[3].set_ylabel('mV')
        axs[3].spines["top"].set_visible(False)
        axs[3].spines["right"].set_visible(False)
        axs[3].plot(T*1000, np.transpose(lfp_v), lw=0.5, label='LFP_V')
        axs[3].plot(T*1000, np.transpose(lfp_), lw=0.5, label='LFP_I')
        axs[3].legend(loc=1)
    
        f.tight_layout() # Fixes the positions of subplots and labels
    
    if PLOT_EXTRA:
        # Second figure. PSP
        f2, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
        
        axs[0].set_title('Pyramidal population (v_pi)')
        axs[0].set_xlabel('Time (ms)')
        axs[0].set_ylabel('IPSP (mV)')
        axs[0].plot(T*1000, np.transpose(v_pi)*1000)
        
        axs[1].set_title('Inhibitory population (v_ip)')
        axs[1].set_xlabel('Time (ms)')
        axs[1].set_ylabel('EPSP (mV)')
        axs[1].plot(T*1000, np.transpose(v_ip)*1000)
    
        axs[0].set_title('LFP (from voltage)')
        axs[0].set_xlabel('Time (ms)')
        axs[0].set_ylabel('mV')
        axs[0].plot(T*1000, np.transpose(lfp_v), label='LFP_V')
        axs[0].legend()
    
        axs[1].set_title('LFP (from current)')
        axs[1].set_xlabel('Time (ms)')
        axs[1].set_ylabel('mV')
        axs[1].plot(T*1000, np.transpose(lfp_), label='LFP_I')
        axs[1].legend()
        
        # Third figure. Membrane potential
        f3, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
        
        axs[0].set_title('Pyramidal Vm (selected cells)')
        axs[0].set_xlabel('Time (ms)')
        axs[0].set_ylabel('mV')
        axs[0].plot(T*1000, np.transpose(Py_monitor.v[0:5])*1e3, lw=0.5, c=c_py)
    
        axs[1].set_title('Interneurons Vm (selected cells)')
        axs[1].set_xlabel('Time (ms)')
        axs[1].set_ylabel('mV')
        axs[1].plot(T*1000, np.transpose(In_monitor.v[0:5])*1e3, lw=0.5, c=c_inter)
        
        f3.tight_layout()
        
        f4, axs = plt.subplots(1, 1, figsize=(6,6))
        axs.plot(np.transpose(v_pi) * 1000, np.transpose(v_ip) * 1000)
        axs.set_xlabel('x1 = V_pi')
        axs.set_ylabel('x3 = V_ip')
        
        plt.show()    
        
        
    #%% Find Coefficient of Variation of ISI
    isi_P = np.zeros(N_activity_plot)
    isi_I = np.zeros(N_activity_plot)
    isis_P = np.zeros(N_activity_plot)
    isis_I = np.zeros(N_activity_plot)
    
    for i in range(N_activity_plot):
        isi_P[i] = mean(diff(sp_P.t[sp_P.i[sp_P.t/second >= 0.5] == i]/second))
        isis_P[i] = std(diff(sp_P.t[sp_P.i[sp_P.t/second >= 0.5] == i]/second))
        
        isi_I[i] = mean(diff(sp_I.t[sp_I.i[sp_I.t/second >= 0.5] == i]/second))
        isis_I[i] = std(diff(sp_I.t[sp_I.i[sp_I.t/second >= 0.5] == i]/second))
        
         
    #%% Save simulation  ------------------------------------------------------------
    if SAVE:    
        P_ = np.array(list(sp_P.spike_trains().values()))
        I_ = np.array(list(sp_I.spike_trains().values()))
        for i in range(0,shape(P_)[0]):
            P_[i] = P_[i]/second
            
        for i in range(0,shape(I_)[0]):
            I_[i] = I_[i]/second
        
        save_dictionary={'LFP': lfp_,
                        'LFP_V': lfp_v,
                        'lfp_dt' : dt_,
                        'v_rest': V_leak,
                        'v_p': mean(Py_monitor.v,0),
                        'v_i': mean(In_monitor.v,0),
                        'v_b': mean(B_monitor.v,0),
                        'v_pi': mean(Py_monitor.v_pi,0),
                        'v_ip': mean(In_monitor.v_ip,0),
                        'R_py': r_P_rate, # 1/diff(np.array(sp_P.t)).mean(),
                        'R_in': r_I_rate,
                        'R_b': r_B_rate,
                        'I_in': I_in,
                        'I_py': I_py,
                        'I_b': I_b,
                        'I_in_tha': I_in_tha,
                        'I_py_tha': I_py_tha,
                        'I_b_tha': I_b_tha,
                        'RECURRENT_PYRAMIDAL': RECURRENT_PYRAMIDAL,
                        'RECURRENT_INHIBITORY': RECURRENT_INHIBITORY,
                        'INHIBIT_INPUT': INHIBIT_INPUT,
                        'ACTIVE_INTERNEURONS': ACTIVE_INTERNEURONS,
                        'input_spike_rate': input_spike_rate,
                        'input_spike_rate_thalamic': input_spike_rate_thalamic,
                        'input_current': input_current}   
    	               
     
        # Save as lfp_last
        # scipy.io.savemat('/data/gpfs/projects/punim0643/artemios/simulations/lfp_last.mat',
        #                  mdict = save_dictionary)
        
       # i = 0
       # while os.path.exists('/data/gpfs/projects/punim0643/artemios/simulations/lfp_%s.mat' % i):
       #     i += 1
       # save_str = format('lfp_%s.png' %(i))
       # scipy.io.savemat('/data/gpfs/projects/punim0643/artemios/simulations/lfp_%s.mat' %(i),
       #                  mdict = save_dictionary)
    
        save_str = format('sweep/lfp_current_%s.png' %(corriente))
        scipy.io.savemat('/data/gpfs/projects/punim0643/artemios/simulations/nonlinearity_three_pop/lfp_current_%s.mat' %(corriente),
                         mdict = save_dictionary)
        
        save_dictionary = None
        I_in = None
        I_py = None
        
    else:
        print('Attention! Results of simulation were not saved. SAVE = False')
    
    ##%% plot  ------------------------------------------------------------
    #f, axs = plt.subplots(4, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    #
    ## colors
    #c_inter = 'C6'  # pink
    #c_py = 'C9'     # light blue
    #c_ex = 'C0'    # blue
    #c_Cor = 'C1'    # orange
    #c_gray = '#e0e0e0' # grey
    #
    ## raster
    ##axs[0].set_title('Raster ({} neurons/pop) N = {}, u = {} spikes/ms'.format(N_activity_plot, N, input_spike_rate))
    #axs[0].set_title('Raster N = {}, u = {} spikes/ms'.format(N, input_spike_rate))
    #axs[0].set_ylabel('Neuron')
    #axs[0].set_yticks([])
    #axs[0].spines["top"].set_visible(False)
    #axs[0].spines["right"].set_visible(False)
    #
    #axs[0].plot(sp_P.t / ms, sp_P.i + 1 * N_I, '.', markersize=2, label='Pyramidal', c=c_py)
    #axs[0].plot(sp_I.t / ms, sp_I.i, '.', markersize=2, label='Inhibitory', c=c_inter)
    #axs[0].legend(loc=1)
    #		
    #axs[1].set_title('Population rates, moving average')
    #axs[1].set_ylabel('Spike rate (Hz)')
    #axs[1].spines["top"].set_visible(False)
    #axs[1].spines["right"].set_visible(False)   
    #	
    #axs[1].plot(r_P.t / ms, r_P_rate / Hz, label='Pyramidal', c=c_py)
    #axs[1].plot(r_I.t / ms, r_I_rate / Hz, label='Interneuron', c=c_inter)
    #axs[1].legend(loc=1)
    #
    ## synaptic currents
    #axs[2].set_title('Synaptic currents')
    #axs[2].set_ylabel('Amplitude (unitless)')
    #axs[2].set_xlabel('Time (ms)')
    #axs[2].spines["top"].set_visible(False)
    #axs[2].spines["right"].set_visible(False)
    ## Others
    #axs[2].plot(T*1000, np.array(st_GABA_I.s_GABA).transpose(), lw=0.5, c=c_gray) # , label='GABA (Inter)'
    #axs[2].plot(T*1000, np.array(st_AMPA_P.s_AMPA).transpose(), lw=0.5, c=c_gray) # , label='AMPA (Py)'
    #axs[2].plot(T*1000, np.array(st_AMPA_cor_I.s_AMPA_cor).transpose(), lw=0.5, c=c_gray) # , label='AMPA_cor (Inter)'
    ## alphas
    #axs[2].plot(T*1000, np.array(st_GABA_P.s_GABA).transpose(), label='GABA (Py)', c=c_py)
    #axs[2].plot(T*1000, np.array(st_AMPA_I.s_AMPA).transpose(), label='AMPA (In)', c=c_inter)
    #axs[2].plot(T*1000, np.array(st_AMPA_cor_P.s_AMPA_cor).transpose(), lw=0.5, c=c_Cor , label='AMPA_cor (Cortical)')
    #axs[2].legend(loc=1)
    #
    ## LFP
    #axs[3].set_title('LFP')
    #axs[3].set_xlabel('Time (ms)')
    #axs[3].set_ylabel('mV')
    #axs[3].spines["top"].set_visible(False)
    #axs[3].spines["right"].set_visible(False)
    #axs[3].plot(T*1000, np.transpose(lfp_v), lw=0.5, label='LFP_V')
    #axs[3].plot(T*1000, np.transpose(lfp_), lw=0.5, label='LFP_I')
    #axs[3].legend(loc=1)
    #
    #f.tight_layout() # Fixes the positions of subplots and labels
    #
    #plt.savefig('/data/gpfs/projects/punim0643/artemios/simulations/' + save_str)
    #print('Results saved as:' + save_str)
    #plt.close('all')

ranges = np.arange(-300, 2500, 25) 
for iterations in ranges:
    brunel(corriente = iterations)

#ranges = np.arange(0, 5.1, 0.1)
#ranges2 = np.arange(-61, 0.1, 1.2) # alpha_ii, 51 elements
#for iterations in ranges:
#    for iterations2 in ranges2:
#        brunel(u = iterations, alpha_ii=iterations2)
#        #clear_cache('cython')
