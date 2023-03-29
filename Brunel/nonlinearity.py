# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 19:51:25 2022

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


from lif_model import set_params, get_equations
import lif_model_CUBN as cubn
import lif_model_COBN as cobn

def brunel(corriente=0, PLOT=False, SAVE=True):
    plt.close('all')
        
            
    #%% Options:
    MODEL           = 'cubn'        # cubn vs cobn
    PARAMS_SOURCE   = 'three_pop'   # 'brunel' or 'allen' or 'three_pop' or''
    
    RECURRENT_PYRAMIDAL     = True  # Self excitation 
    RECURRENT_INHIBITORY    = True  # Self inhibition
    ACTIVE_INTERNEURONS     = True  # Inhibitory population
    ACTIVE_GABAb            = True # Second inhibitory (slow) population (Wendling-like model)
    INHIBIT_INPUT           = True # Excitatory cortical input to inhibitory population
    
    GABA_A_MULTIPLIER = 1           # GABA_A Agonist applied for the whole duration
    MIDWAY_MULTIPLIER = 1           # GABA_A Agonist applied midsimulation
    
    GAUSSIAN_REFRACTORY = False      # If true, the refractory period of each cell is taken from a gaussian distribution, otherwise it is the same for all
    GAUSSIAN_THRESHOLD  = False      # If true, the refractory period of each cell is taken from a gaussian distribution, otherwise it is the same for all
    
    # SAVE = False                    # Save ground truth data
    # PLOT = True                     # Plot results 
    STATS = False                    # Calculate spike statistics (ISI distance, CV, etc)
    
    PSP_FR   = 0                    # Presynaptic firing rate for TEST_PSP (TEST_PSP needs to be diff to none for this to take effect)                               
    TEST_PSP = 'none'               # Testing the post synaptic potential of given synapses to a specified input firing rate. Options: 'pu', 'pp', 'pi', 'ii', 'ip', 'bp', 'bi', 'pb', 'none'. To prevent neurons spiking, make V_thr large.
    RUNTYPE  = 'normal'             # Simulation: 'normal', 'current_pulse' or 'gaba_agonist'
    
    #corriente = 0#500 #50
    # Balanced-rate network (?) with input currents: 
    input_current   = corriente    # Injected current to Pyramidal population # Use this to calculate the nonlinearity (Vm -> Spike_rate sigmoid) on the disconnected model
    input_current_I = corriente  # Inhibitory interneurons
    
    input_spike_rate = [0]#[0, 0.25]#[0, 1, 3, 5] #[u] #[5] #  [0, 2.5, 5] # spikes/ms/cell (driving input)
    input_spike_rate_thalamic = 1 #1.5 # spikes/ms/cell (spontaneous activity)
    input_spike_rate_thalamic_in = 1 #1.5 # spikes/ms/cell (spontaneous activity)
    
    #%% parameters  --------------------------------------------------------------
    simulation_time = 1 * second
    dt_ = 100 * usecond
    T = np.linspace(0, simulation_time, round(simulation_time/dt_)) # Time vector for plots (in seconds)
       
    # populations
    N = 2000
    N_P = 1#int(N * 4)  # pyramidal neurons
    # N_I = int(N * (0.5 ** ACTIVE_GABAb))    # interneurons
    N_I = 1#int(N)    # interneurons
    
    # set populations parameters
    if MODEL == 'cubn':
        params_py = cubn.set_params('pyramidal')
        params_in = cubn.set_params('inhibitory')
        
    elif MODEL == 'cobn':
        params_py = cobn.set_params('pyramidal')
        params_in = cobn.set_params('inhibitory')
    else:
        raise Exception('Model %s does not exist (the options are cubn and cobn)' %MODEL)
        
    if not(ACTIVE_GABAb):
        params_py = set_params('pyramidal', 'allen')
        params_in = set_params('inhibitory', 'allen')
    
    # # Modify parameter according to the arguments ('parameter' and 'value_')
    # if parameter != '':
    #     if pop_=='py':
    #         params_py[parameter] = params_py.get(parameter) * value_
    #     elif pop_=='in':
    #         params_in[parameter] = params_in.get(parameter) * value_
        
    
    # Probability of connection
    p_IP = params_py.get('p_IP')
    p_BP = params_py.get('p_BP')
    p_PI = params_py.get('p_PI')
    p_PP = params_py.get('p_PP')
    p_II = params_py.get('p_II')
    
    # voltage
    V_leak = -70. * mV      # Resting membrane potential
    V_reset = -59 * mV      # Reset voltage. Equal to V_leak-> To use Burkitt's, 2006 Eq. (12)
    v_AMPA = 0 * mV         # AMPA reversal potential
    v_GABA = -80 * mV       # GABA_A reversal potential
    v_GABAb = -90 * mV      # GABA_B reversal potential
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
    if MODEL == 'cubn':
        # AMPA (excitatory)
        j_AMPA_rec_P = params_py.get('j_AMPA') * 2000/N 
        j_AMPA_rec_I = params_in.get('j_AMPA') * 2000/N 
            
        j_AMPA_cor_P = params_py.get('j_AMPA_ext')
        j_AMPA_cor_I = params_in.get('j_AMPA_ext')
        
        j_AMPA_tha_P = params_py.get('j_AMPA_tha')
        j_AMPA_tha_I = params_in.get('j_AMPA_tha')
        
        # GABAergic (inhibitory)
        j_GABA_P    = GABA_A_MULTIPLIER * params_py.get('j_GABA') * 2000/N
        j_GABA_I    = GABA_A_MULTIPLIER * params_in.get('j_GABA') * 2000/N 
        j_GABAb_P   = params_py.get('j_GABAb') * 2000/N 
        
    elif MODEL == 'cobn': 
        # AMPA (excitatory)
        g_AMPA_rec_P = params_py.get('g_AMPA') * 2000/N
        g_AMPA_rec_I = params_in.get('g_AMPA') * 2000/N 
            
        g_AMPA_cor_P = params_py.get('g_AMPA_ext')
        g_AMPA_cor_I = params_in.get('g_AMPA_ext')
        
        g_AMPA_tha_P = params_py.get('g_AMPA_tha')
        g_AMPA_tha_I = params_in.get('g_AMPA_tha')
        
        # GABAergic (inhibitory)
        g_GABA_P    = GABA_A_MULTIPLIER * params_py.get('g_GABA') * 2000/N 
        g_GABA_I    = GABA_A_MULTIPLIER * params_in.get('g_GABA') * 2000/N
        g_GABAb_P   = params_py.get('g_GABAb') * 2000/N
        
    else:
        raise Exception('Model %s does not exist' %MODEL)
    
    
    
    
    
    # Weight constants. Amplitude of the synaptic input
    # Pyramidal 
    increment_AMPA_P        = params_py.get('weight')
    increment_AMPA_ext_P    = params_py.get('external_input_weight')
    increment_GABA_P        = params_py.get('weight')
    
    # Inhibitory interneurons
    increment_AMPA_I        = params_in.get('weight')
    increment_AMPA_ext_I    = params_in.get('external_input_weight')
    increment_GABA_I        = params_in.get('weight')
    
    
    # Injected current
    I_injected      = -input_current * pA # Input current to Pyramidal population. Sets a baseline spiking rate
    I_injected_I    = -input_current_I * pA # Input current to Pyramidal population. Sets a baseline spiking rate
    
    # I_injected = np.random.normal(loc=-input_current, scale=10, size=int(l))
    # I_injected[:] = I_injected
    # I_injected = TimedArray(I_injected * pA, dt=dt_)
    
    #%% modeling  ----------------------------------------------------------------
    # model equations
    if MODEL == 'cubn':
        eqs_P = cubn.get_equations('pyramidal')
        eqs_I = cubn.get_equations('inhibitory')
    elif MODEL == 'cobn':
        eqs_P = cobn.get_equations('pyramidal')
        eqs_I = cobn.get_equations('inhibitory')
    else: 
        raise Exception('Model %s does not exist' %MODEL)
    
    
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
                                                                    v_iu = V_reset-V_leak
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
        
        
    
    # Custom Poisson population to test PSP on a given synapse
    if (TEST_PSP == 'ip') | (TEST_PSP == 'pp') | (TEST_PSP == 'bp'):
        N_PSP_Test = N_P
    elif (TEST_PSP == 'pi') | (TEST_PSP == 'ii') | (TEST_PSP == 'bi'):
        N_PSP_Test = N_I
    elif (TEST_PSP == 'ib') | (TEST_PSP == 'pb'):
        N_PSP_Test = N_I
    else:
        N_PSP_Test = 0    
    Pop_PSP_Test = PoissonGroup(N_PSP_Test, rates = PSP_FR * Hz, dt=dt_) # poisson input
        
    # Ornstein-Uhlenbeck process:
    # sigma_sq_noise = 0.16 * volt * volt   # Ornstein-Uhlenbeck process (cortico-cortical) units need to be volt^2 for units in the population's equation to be consistent. Makes sense, because sigma should be in volts, hence sigma^2 in volts^2. Check cavalleri 2014
    # tau_noise = 16 * ms                   # Ornstein-Uhlenbeck 
    # input_cortical = NeuronGroup(num_inputs, 'dv/dt = -v/tau_noise + np.sqrt(2 * sigma_sq_noise/tau_noise) * xi : volt', threshold = 'v >= 0.5 * volt', reset = 'v = 0 * volt', dt=dt_, method='euler')
    
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
    
    eqs_pre_gaba_B = '''
    s_GABA += increment_GABA_I
    '''
    
    eqs_pre_glut_B = '''
    s_AMPA += increment_AMPA_I
    '''
    
       
    # Synapses
    # P to P
    C_P_P = Synapses(Py_Pop, Py_Pop, on_pre=eqs_pre_glut_P, method='rk4', dt=dt_, delay=delay, name='synapses_pp')
    C_P_P.connect('i != j', p = p_PP)
    C_P_P.active = False    # Testing no recursive connections to match NMM
    
    # P to I
    C_P_I = Synapses(Py_Pop, In_Pop, on_pre=eqs_pre_glut_I, method='rk4', dt=dt_, delay=delay, name='synapses_ip')
    C_P_I.connect(p = p_PI)     
    C_P_I.active = False
    
    # GABAa to I
    C_I_I = Synapses(In_Pop, In_Pop, on_pre=eqs_pre_gaba_I, method='rk4', dt=dt_, delay=delay, name='synapses_ii')
    C_I_I.connect('i != j', p = p_II)
    C_I_I.active = False
    
    # GABAa to P
    C_I_P = Synapses(In_Pop, Py_Pop, on_pre=eqs_pre_gaba_P, method='rk4', dt=dt_, delay=delay, name='synapses_pi')
    C_I_P.connect(p = p_IP)    
    C_I_P.active = False
    
    # GABAb to P
    C_B_P = Synapses(In_Pop, Py_Pop, on_pre=eqs_pre_gabab_P, method='rk4', dt=dt_, delay=delay, name='synapses_pb')
    C_B_P.connect(p = p_BP)
    C_B_P.active = False
    
    
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
        
    # Poisson input (Cortico-cortical) input to inhibitory interneurons. Controlled by INHIBIT_INPUT
    C_Cor_I = PoissonInput(In_Pop, 's_AMPA_cor', num_inputs, (input_spike_rate[0]*1000/num_inputs) * Hz, increment_AMPA_ext_I)
    C_Cor_I.active = INHIBIT_INPUT # Innactive cortico-cortical -> interneuron
    
    # Poisson input (Thalamic, baseline spike rate)
    C_Tha_P = PoissonInput(Py_Pop, 's_AMPA_tha', num_inputs, (input_spike_rate_thalamic*1000/num_inputs) * Hz, increment_AMPA_ext_P)
    C_Tha_I = PoissonInput(In_Pop, 's_AMPA_tha', num_inputs, (input_spike_rate_thalamic_in*1000/num_inputs) * Hz, increment_AMPA_ext_I)
    
    
    
    # Testing PSP on chosen synapses
    ACTIVE_TEST = True
    if (TEST_PSP == 'pi'):
        C_Test_PSP = Synapses(Pop_PSP_Test, Py_Pop, on_pre=eqs_pre_gaba_P, method='rk4', dt=dt_, delay=delay, name='synapses_pi_test')
        connection_probability = p_IP
    
    elif (TEST_PSP == 'pp'):
        C_Test_PSP = Synapses(Pop_PSP_Test, Py_Pop, on_pre=eqs_pre_glut_P, method='rk4', dt=dt_, delay=delay, name='synapses_pp_test')
        connection_probability = p_PP
    
    elif (TEST_PSP == 'ip'):
        C_Test_PSP = Synapses(Pop_PSP_Test, In_Pop, on_pre=eqs_pre_glut_I, method='rk4', dt=dt_, delay=delay, name='synapses_ip_test')
        connection_probability = p_PI
    
    elif (TEST_PSP == 'ii'):
        C_Test_PSP = Synapses(Pop_PSP_Test, In_Pop, on_pre=eqs_pre_gaba_I, method='rk4', dt=dt_, delay=delay, name='synapses_ii_test')
        connection_probability = p_II
    
    elif (TEST_PSP == 'bp'):
        C_Test_PSP = Synapses(Pop_PSP_Test, B_Pop, on_pre=eqs_pre_glut_B, method='rk4', dt=dt_, delay=delay, name='synapses_bp_test')
        connection_probability = p_PI
    
    elif (TEST_PSP == 'pb'):
        C_Test_PSP = Synapses(Pop_PSP_Test, Py_Pop, on_pre=eqs_pre_gabab_P, method='rk4', dt=dt_, delay=delay, name='synapses_pb_test')
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
    
    r_P = PopulationRateMonitor(Py_Pop)#[0:N_activity_plot])
    r_I = PopulationRateMonitor(In_Pop)
    
    # st_AMPA_P = StateMonitor(Py_Pop, ('s_AMPA'), record = 0)
    # st_GABA_P = StateMonitor(Py_Pop, 's_GABA', record = 0)
    # st_AMPA_cor_P = StateMonitor(Py_Pop, 's_AMPA_cor', record = 0)
    
    # st_AMPA_I = StateMonitor(In_Pop, 's_AMPA', record = 0)
    # st_GABA_I = StateMonitor(In_Pop, 's_GABA', record = 0)
    # st_AMPA_cor_I = StateMonitor(In_Pop, 's_AMPA_cor', record = 0)
    
    Py_monitor = StateMonitor(Py_Pop, ['I_GABA_rec', 'I_GABAb', 'I_AMPA_tha', 'v',  'v_pi', 'I_tot'], record = True) # Monitoring the AMPA and GABA currents in the Pyramidal population
    In_monitor = StateMonitor(In_Pop, ['I_AMPA_rec', 'I_AMPA_tha', 'v', 'v_ip', 'I_tot'], record = True)
    B_monitor = []
    
    # If TEST_PSP is on, create a new monitor just to observe that PSP
    if TEST_PSP[0]=='p':
        PSP_monitor = StateMonitor(Py_Pop, ['v_'+TEST_PSP], record = True)
    elif TEST_PSP[0]=='i':
        PSP_monitor = StateMonitor(In_Pop, ['v_'+TEST_PSP], record = True)
    
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
        net.run(0.49 * second, report='stdout') # Run first segment, if running more segments, run for a fraction of simulation_time
        
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
        # print(colored('The simulation did not run. RUNTYPE most be one of the three options.', 'yellow'))
        raise Exception("The simulation did not run. RUNTYPE most be one of the three options.")
        
        
        
        
        
        
        
        
        
        
        
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
    
        
    # Calculate mean PSP (NMM states)
    v_pi = mean(Py_monitor.v_pi, 0)
    v_ip = mean(In_monitor.v_ip, 0)
    
    # Mean currents
    i_pi = np.transpose(mean(Py_monitor.I_GABA_rec,0)) / amp
    i_pb = np.transpose(mean(Py_monitor.I_GABAb,0)) / amp
    i_ip = np.transpose(mean(In_monitor.I_AMPA_rec,0)) / amp
    
    # Local Field Potential
    # current based;
    # lfp = sum((Py_monitor.I_GABAb + Py_monitor.I_GABA_rec),0) - sum((Py_monitor.I_AMPA_cor + Py_monitor.I_AMPA_rec),0) # Difference of currents
    # lfp_ = lfp / g_m_P # Sum across all Pyramidal neurons and divide by the leak conductance to get volts
    
    # voltage based:
    mean_v_Py = np.transpose(mean(Py_monitor.v,0) - V_leak) * 1e3
    lfp_v = mean_v_Py/volt 
    
    v_p = np.transpose(Py_monitor.v[0:5])*1e3
    v_i = np.transpose(In_monitor.v[0:5])*1e3
           

    #%% plotting  -----------------------------------------------------------------
    if PLOT:
        f, axs = plt.subplots(4, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
        
        # colors
        c_inter = 'C6'  # pink
        c_py = 'C9'     # light blue
        c_ex = 'C0'    # blue
        c_Cor = 'C1'    # orange
        c_gray = '#e0e0e0' # grey
        
        # raster
        axs[0].set_title('Raster ({} neurons/pop) N = {}, u = {} spikes/ms'.format(N_activity_plot, N, input_spike_rate))
        axs[0].set_ylabel('Neuron')
        axs[0].set_yticks([])
        axs[0].spines["top"].set_visible(False)
        axs[0].spines["right"].set_visible(False)
        
        # axs[0].plot(sp_E.t / ms, sp_E.i + 2 * N_activity_plot, '.', markersize=2, label='Spiny', c=c_ex)
        axs[0].plot(sp_P.t / ms, sp_P.i + 1 * N_I, '.', markersize=2, label='Pyramidal', c=c_py)
        axs[0].plot(sp_I.t / ms, sp_I.i, '.', markersize=2, label='Inhibitory', c=c_inter)
        axs[0].legend(loc=1)
                
        axs[1].set_title('Population rates, moving average')
        axs[1].set_ylabel('Spike rate (Hz)')
        axs[1].spines["top"].set_visible(False)
        axs[1].spines["right"].set_visible(False)   
            
        axs[1].plot(r_P.t / ms, r_P_rate / Hz, label='Pyramidal', c=c_py)
        axs[1].plot(r_I.t / ms, r_I_rate / Hz, label='Interneuron', c=c_inter)
        # axs[1].plot(r_Cor.t / ms, r_Cor_rate / Hz, label='Cortico-cortical (OU)', c=c_Cor)
        # axs[1].plot(T_u / ms, u / Hz)
        axs[1].legend(loc=1)
        
        # synaptic currents
        # axs[2].set_title('Synaptic currents')
        # axs[2].set_ylabel('Amplitude (unitless)')
        # axs[2].set_xlabel('Time (ms)')
        # axs[2].spines["top"].set_visible(False)
        # axs[2].spines["right"].set_visible(False)
        # # Others
        # axs[2].plot(T*1000, np.array(st_GABA_I.s_GABA).transpose(), lw=0.5, c=c_gray) # , label='GABA (Inter)'
        # axs[2].plot(T*1000, np.array(st_AMPA_P.s_AMPA).transpose(), lw=0.5, c=c_gray) # , label='AMPA (Py)'
        # axs[2].plot(T*1000, np.array(st_AMPA_cor_I.s_AMPA_cor).transpose(), lw=0.5, c=c_gray) # , label='AMPA_cor (Inter)'
        # # alphas
        # axs[2].plot(T*1000, np.array(st_GABA_P.s_GABA).transpose(), label='GABA (Py)', c=c_py)
        # axs[2].plot(T*1000, np.array(st_AMPA_I.s_AMPA).transpose(), label='AMPA (In)', c=c_inter)
        # axs[2].plot(T*1000, np.array(st_AMPA_cor_P.s_AMPA_cor).transpose(), lw=0.5, c=c_Cor , label='AMPA_cor (Cortical)')
        # axs[2].legend(loc=1)
        
        # # LFP
        # axs[3].set_title('LFP')
        # axs[3].set_xlabel('Time (ms)')
        # axs[3].set_ylabel('mV')
        # axs[3].spines["top"].set_visible(False)
        # axs[3].spines["right"].set_visible(False)
        # axs[3].plot(T*1000, np.transpose(lfp_v), lw=0.5, label='LFP_V')
        # axs[3].plot(T*1000, np.transpose(lfp_), lw=0.5, label='LFP_I')
        # axs[3].legend(loc=1)
    
        f.tight_layout() # Fixes the positions of subplots and labels
    
    
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
    
        # axs[0].set_title('LFP (from voltage)')
        # axs[0].set_xlabel('Time (ms)')
        # axs[0].set_ylabel('mV')
        # axs[0].plot(T*1000, np.transpose(lfp_v), label='LFP_V')
        # axs[0].legend()
    
        # axs[1].set_title('LFP (from current)')
        # axs[1].set_xlabel('Time (ms)')
        # axs[1].set_ylabel('mV')
        # axs[1].plot(T*1000, np.transpose(lfp_), label='LFP_I')
        # axs[1].legend()
        
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
    
        save_dictionary={'LFP_V': lfp_v,
                        'lfp_dt' : dt_,
                        'Vm': -(I_injected/g_m_P), # To calculate the nonlinearity, need to simulate single cell disconnected network 
                        'Vm_interneurons': -(I_injected_I/g_m_I), # To calculate the nonlinearity, need to simulate single cell disconnected network 
                        'R_py': r_P_rate, # 1/diff(np.array(sp_P.t)).mean(),
                        'R_in': r_I_rate,
                        'I_py': Py_monitor.I_tot,
                        'I_in': In_monitor.I_tot,
                        'I_py_tha': Py_monitor.I_AMPA_tha,
                        'I_in_tha': In_monitor.I_AMPA_tha,
                        'RECURRENT_PYRAMIDAL': RECURRENT_PYRAMIDAL,
                        'RECURRENT_INHIBITORY': RECURRENT_INHIBITORY,
                        'INHIBIT_INPUT': INHIBIT_INPUT,
                        'ACTIVE_INTERNEURONS': ACTIVE_INTERNEURONS,
                        'input_spike_rate': input_spike_rate,
                        'input_spike_rate_thalamic': input_spike_rate_thalamic,
                        'input_spike_rate_thalamic_in': input_spike_rate_thalamic_in,
                        'input_current': input_current}   
            
        # Save as lfp_last
        scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/lfp_last.mat',
                         mdict = save_dictionary)
        
        i = 0
        while os.path.exists('C://Users/artemios/Documents/Multiscale_Models_Data/2023/nonlinearity/one_tha/lfp_inputCurrent_%s_%s.mat' % (floor(input_current),i)):
            i += 1
        
        scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/2023/nonlinearity/one_tha/lfp_inputCurrent_%s_%s.mat' % (floor(input_current),i),
                          mdict=save_dictionary)
    
        # scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/lfp_%s.mat' % i,
        #                   mdict = save_dictionary)
        
        # # Save at the end to keep all recordings
        # if RECURRENT_PYRAMIDAL:
        #     i = 0
        #     while os.path.exists('simulations/CUBN/recurrent_excitation/inhibitory_input_lfp_inputRate_%s_%s.mat' % (floor(timed_rate.values[0]),i)):
        #         i += 1
            
        #     scipy.io.savemat('simulations/CUBN/recurrent_excitation/inhibitory_input_lfp_inputRate_%s_%s.mat' % (floor(timed_rate.values[0]),i),
        #                       mdict=save_dictionary)
            
        # elif RECURRENT_INHIBITORY:
        #     i = 0
        #     while os.path.exists('simulations/CUBN/recurrent_inhibition/inhibitory_input_lfp_inputRate_%s_%s.mat' % (floor(timed_rate.values[0]),i)):
        #         i += 1
            
        #     scipy.io.savemat('simulations/CUBN/recurrent_inhibition/inhibitory_input_lfp_inputRate_%s_%s.mat' % (floor(timed_rate.values[0]),i),
        #                       mdict=save_dictionary)
        # elif ACTIVE_INTERNEURONS:
        #     if INHIBIT_INPUT:
        #         i = 0
        #         while os.path.exists('simulations/CUBN/no_recurrent_connections/inhibitory_input/inhibitory_input_lfp_inputRate_%s_%s.mat' % (floor(timed_rate.values[0]),i)):
        #             i += 1
                
        #         scipy.io.savemat('simulations/CUBN/no_recurrent_connections/inhibitory_input/inhibitory_input_lfp_inputRate_%s_%s.mat' % (floor(timed_rate.values[0]),i),
        #                           mdict=save_dictionary)
        #     else:                
        #         i = 0
        #         while os.path.exists('simulations/CUBN/no_recurrent_connections/iterneurons_lfp_inputRate_%s_%s.mat' % (floor(timed_rate.values[0]), i)):
        #             i += 1
                
        #         scipy.io.savemat('simulations/CUBN/no_recurrent_connections/iterneurons_lfp_inputRate_%s_%s.mat' % (floor(timed_rate.values[0]),i),
        #                           mdict=save_dictionary)
        # else:
        #     i = 0
        #     while os.path.exists('simulations/CUBN/injected_current/twice_ref/lfp_injected_current_%s_%s.mat' % (floor(input_current),  i)):
        #         i += 1
                
        #     scipy.io.savemat('simulations/CUBN/injected_current/twice_ref/lfp_injected_current_%s_%s.mat' % (floor(input_current), i),
        #                       mdict=save_dictionary)
        
    else:
        print(colored('Attention! Results of simulation were not saved. SAVE = False', 'yellow'))
    
    

# Run iteratively. Need to uncomment the def line at the start of the file.
# ranges = np.arange(500,500.01,0.001)#(-1500, 1500, 25)#(-300, 1500, 5) 
ranges = np.arange(-1500, 1500, 5)#(-300, 1500, 5) 
for iterations in ranges:
    brunel(corriente = iterations)

# brunel(corriente = 0, PLOT=True, SAVE = False)