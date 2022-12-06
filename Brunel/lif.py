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
from termcolor import colored  # Coloured text in the terminal
import matplotlib.pyplot as plt
prefs.codegen.target = 'numpy'  # use the Python fallback instead of C compilation
devices.device.shape = []       # This and the following line remove an annoying warning when brian2 is imported and loaded into RAM
devices.device.size = []

from lif_model import set_params, get_equations

# def brunel(u=0):
# plt.close('all')
    
#%% Options:
RECURRENT_PYRAMIDAL = True    # Self excitation 
RECURRENT_INHIBITORY = True   # Self inhibition
PARAMS_SOURCE = 'three_pop'       # 'brunel' or 'allen' or 'three_pop'
ACTIVE_INTERNEURONS = True    # Inhibitory population
ACTIVE_GABAb = True           # Second inhibitory (slow) population (Wendling-like model)
INHIBIT_INPUT = False         # Excitatory cortical input to inhibitory population
GAUSSIAN_REFRACTORY = False    # If true, the refractory period of each cell is taken from a gaussian distribution, otherwise it is the same for all
GAUSSIAN_THRESHOLD = False     # If true, the refractory period of each cell is taken from a gaussian distribution, otherwise it is the same for all
SAVE = False                   # Save ground truth data
PLOT = True                   # Plot results (main Figure)
PLOT_EXTRA = False             # Plot extra things.
PSP_FR = 0                    # Presynaptic firing rate for TEST_PSP (TEST_PSP needs to be diff to none for this to take effect)                               
TEST_PSP = 'none'             # Testing the post synaptic potential of given synapses to a specified input firing rate. Options: 'pu', 'pp', 'pi', 'ii', 'ip', 'bp', 'bi', 'pb', 'none'. To prevent neurons spiking, make V_thr large.
GABA_A_MULTIPLIER = 1         # GABA_A Agonist applied for the whole duration
MIDWAY_MULTIPLIER = 1         # GABA_A Agonist applied midsimulation
RUNTYPE = 'normal'            # Simulation: 'normal', 'current_pulse' or 'gaba_agonist'

corriente = 0#50
# Balanced-rate network (?) with input currents: Py = 500.01 pA, In = 398 pA
input_current = corriente  # 437.5 # 500.01       # Injected current to Pyramidal population # Use this to calculate the nonlinearity (Vm -> Spike_rate sigmoid) on the disconnected model
input_current_I = corriente # 350 # 398 # 400.01     # Inhibitory interneurons

input_spike_rate = [2]#[0, 2.5, 5, 7.5]#[0, 1, 3, 5] #[u] #[5] #  [0, 2.5, 5] # spikes/ms/cell (driving input)
input_spike_rate_thalamic = 2 # spikes/ms/cell (spontaneous activity)

#%% parameters  --------------------------------------------------------------
simulation_time = 2 * second
dt_ = 100 * usecond
T = np.linspace(0, simulation_time, round(simulation_time/dt_)) # Time vector for plots (in seconds)
   
# populations
N = 2000 # 135 # 675
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
increment_AMPA_P =  params_py.get('weight')
increment_AMPA_ext_P = params_py.get('external_input_weight')
increment_GABA_P = params_py.get('weight')

# Inhibitory interneurons
increment_AMPA_I = params_in.get('weight')
increment_AMPA_ext_I = params_in.get('external_input_weight')
increment_GABA_I = params_in.get('weight')


# Alpha function's parameter (and double exponential) to fix the units in ds/dt
k = 1 / ms # Dimmensionless?, check Nicola and Campbell 2013


# Injected current
I_injected = -input_current * pA # Input current to Pyramidal population. Sets a baseline spiking rate
I_injected_I = -input_current_I * pA # Input current to Pyramidal population. Sets a baseline spiking rate

# I_injected = np.random.normal(loc=-input_current, scale=10, size=int(l))
# I_injected[:] = I_injected
# I_injected = TimedArray(I_injected * pA, dt=dt_)

#%% modeling  ----------------------------------------------------------------
# model equations
eqs_P = get_equations('pyramidal')
eqs_I = get_equations('inhibitory')
eqs_B = get_equations('gabab')


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
    
# Ornstein-Uhlenbeck process:
# sigma_sq_noise = 0.16 * volt * volt   # Ornstein-Uhlenbeck process (cortico-cortical) units need to be volt^2 for units in the population's equation to be consistent. Makes sense, because sigma should be in volts, hence sigma^2 in volts^2. Check cavalleri 2014
# tau_noise = 16 * ms                   # Ornstein-Uhlenbeck 
# input_cortical = NeuronGroup(num_inputs, 'dv/dt = -v/tau_noise + np.sqrt(2 * sigma_sq_noise/tau_noise) * xi : volt', threshold = 'v >= 0.5 * volt', reset = 'v = 0 * volt', dt=dt_, method='euler')

#%% synaptic equations -------------------------------------------------------
# Pyramidal (self-excitatory)
# eqs_glut_P = '''
# s_AMPA_post = s_AMPA_syn : 1 (summed)
# ds_AMPA_syn / dt = - s_AMPA_syn / tau_s_AMPA_P + k * x : 1 (clock-driven)
# dx / dt = - x / tau_s_AMPA_P : 1 (clock-driven)
# '''
eqs_pre_glut_P = '''
s_AMPA += increment_AMPA_P
'''

# Pyramidal (gabaergic, inhibitory to pyramidal)
# eqs_gaba_P = '''
# s_GABA_post = s_GABA_syn : 1 (summed)
# ds_GABA_syn / dt = - s_GABA_syn / tau_s_GABA_P + k * x : 1 (clock-driven)
# dx / dt = - x / tau_s_GABA_P : 1 (clock-driven)
# '''
eqs_pre_gaba_P = '''
s_GABA += increment_GABA_P
'''

# GABAb to Py
eqs_pre_gabab_P = '''
s_GABAb += increment_GABA_P
'''

# Interneurons (glutamate, pyramidal to inhibitory)
# eqs_glut_I = '''
# s_AMPA_post = s_AMPA_syn : 1 (summed)
# ds_AMPA_syn / dt = - s_AMPA_syn / tau_s_AMPA_I + k * x : 1 (clock-driven)
# dx / dt = - x / tau_s_AMPA_I : 1 (clock-driven)
# '''
eqs_pre_glut_I = '''
s_AMPA += increment_AMPA_I
'''

# Interneurons (recurrent inhibiton)
# eqs_gaba_I = '''
# s_GABA_post = s_GABA_syn : 1 (summed)
# ds_GABA_syn / dt = - s_GABA_syn / tau_s_GABA_I + k * x : 1 (clock-driven)
# dx / dt = - x / tau_s_GABA_I : 1 (clock-driven)
# '''
eqs_pre_gaba_I = '''
s_GABA += increment_GABA_I
'''

# GABAb (inhibition)
eqs_pre_gaba_B = '''
s_GABA += increment_GABA_I
'''
# GABAb (excitation)
eqs_pre_glut_B = '''
s_AMPA += increment_AMPA_I
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
C_I_I.active = RECURRENT_INHIBITORY & (not ACTIVE_GABAb)

# I to P
C_I_P = Synapses(In_Pop, Py_Pop, on_pre=eqs_pre_gaba_P, method='rk4', dt=dt_, delay=delay, name='synapses_pi')
C_I_P.connect(p = p_IP)    
C_I_P.active = ACTIVE_INTERNEURONS

# P to GABAb
C_P_B = Synapses(Py_Pop, B_Pop, on_pre=eqs_pre_glut_B, method='rk4', dt=dt_, delay=delay, name='synapses_bp')
C_P_B.connect(p = p_PI)     
C_P_B.active = ACTIVE_GABAb

# I to GABAb
C_I_B = Synapses(In_Pop, B_Pop, on_pre=eqs_pre_gaba_B, method='rk4', dt=dt_, delay=delay, name='synapses_bi')
C_I_B.connect(p = p_II)
C_I_B.active = ACTIVE_GABAb

# GABAb to P
C_B_P = Synapses(B_Pop, Py_Pop, on_pre=eqs_pre_gabab_P, method='rk4', dt=dt_, delay=delay, name='synapses_pb')
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
    
# Poisson input (Cortico-cortical) input to inhibitory interneurons. Controlled by INHIBIT_INPUT
C_Cor_I = PoissonInput(In_Pop, 's_AMPA_cor', num_inputs, (input_spike_rate[0]*1000/num_inputs) * Hz, increment_AMPA_ext_I)
C_Cor_I.active = INHIBIT_INPUT # Innactive cortico-cortical -> interneuron

# Poisson input (Thalamic, baseline spike rate)
C_Tha_P = PoissonInput(Py_Pop, 's_AMPA_tha', num_inputs, (input_spike_rate_thalamic*1000/num_inputs) * Hz, increment_AMPA_ext_P)
C_Tha_I = PoissonInput(In_Pop, 's_AMPA_tha', num_inputs, (input_spike_rate_thalamic*1000/num_inputs) * Hz, increment_AMPA_ext_I)
C_Tha_B = PoissonInput(B_Pop, 's_AMPA_tha', num_inputs, (input_spike_rate_thalamic*1000/num_inputs) * Hz, increment_AMPA_ext_I)

# Poisson population
# C_Cor_P = Synapses(Pop_Cor, Py_Pop, model=eqs_cor_P, on_pre=eqs_pre_cor_P, method='rk4', dt=dt_, delay=delay, name='synapses_pext')
# C_Cor_P.connect(p = 1)    
# C_Cor_I = Synapses(Pop_Cor, In_Pop, model=eqs_cor_I, on_pre=eqs_pre_cor_I, method='rk4', dt=dt_, delay=delay, name='synapses_iext')
# C_Cor_I.connect(p = 1)    



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
elif (TEST_PSP == 'bi'):
    C_Test_PSP = Synapses(Pop_PSP_Test, B_Pop, on_pre=eqs_pre_gaba_B, method='rk4', dt=dt_, delay=delay, name='synapses_bi_test')
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
sp_B = SpikeMonitor(B_Pop[:]) #[:N_activity_plot])


r_P = PopulationRateMonitor(Py_Pop)#[0:N_activity_plot])
r_I = PopulationRateMonitor(In_Pop)
r_B = PopulationRateMonitor(B_Pop)

st_AMPA_P = StateMonitor(Py_Pop, ('s_AMPA'), record = 0)
st_GABA_P = StateMonitor(Py_Pop, 's_GABA', record = 0)
st_AMPA_cor_P = StateMonitor(Py_Pop, 's_AMPA_cor', record = 0)

st_AMPA_I = StateMonitor(In_Pop, 's_AMPA', record = 0)
st_GABA_I = StateMonitor(In_Pop, 's_GABA', record = 0)
st_AMPA_cor_I = StateMonitor(In_Pop, 's_AMPA_cor', record = 0)

Py_monitor = StateMonitor(Py_Pop, ['I_AMPA_cor', 'I_AMPA_rec', 'I_GABA_rec', 'I_AMPA_spi', 'I_GABAb', 'v',  'v_pi', 'I_tot', 'v_pu'], record = True) # Monitoring the AMPA and GABA currents in the Pyramidal population
In_monitor = StateMonitor(In_Pop, ['v', 'v_ip', 'I_tot'], record = True)
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
lfp = sum((Py_monitor.I_GABAb + Py_monitor.I_GABA_rec),0) - sum((Py_monitor.I_AMPA_cor + Py_monitor.I_AMPA_rec + Py_monitor.I_AMPA_spi),0) # Difference of currents
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
    
 
    # Testing PSP on chosen synapses
    if (TEST_PSP == 'pi'):
        v_xx = mean(Py_monitor.v_pi, 0)
    elif (TEST_PSP == 'pp'):
        v_xx = mean(Py_monitor.v_pp, 0)
    elif (TEST_PSP == 'ip'):
        v_xx = mean(In_monitor.v_ip, 0)
    elif (TEST_PSP == 'ii'):
        v_xx = mean(In_monitor.v_ii, 0)
    elif (TEST_PSP == 'bi'):
        v_xx = mean(B_monitor.v_bi, 0)
    elif (TEST_PSP == 'bp'):
        v_xx = mean(B_monitor.v_bp, 0)
    elif (TEST_PSP == 'pb'):
        v_xx = mean(Py_monitor.v_pb, 0)
    elif (TEST_PSP == 'pu'):
        v_xx = mean(Py_monitor.v_pu, 0)
    else:
        v_xx = T
 
    axs[0].set_xlabel('Time (ms)')
    axs[0].set_ylabel('IPSP (mV)')
    axs[0].plot(T*1000, np.transpose(v_xx)*1000)
    
    axs[1].plot( T*1000, 1000 * (-V_leak + mean(Py_monitor.v,0)),lw=0.5, label='Pyramidal')
    axs[1].plot( T*1000, 1000 * (-V_leak + mean(In_monitor.v,0)),lw=0.5, label='GABA_A')
    axs[1].plot( T*1000, 1000 * (-V_leak + mean(B_monitor.v,0)),lw=0.5, label='GABA_B')
    axs[1].legend()
    
    # axs[1].set_title('Inhibitory population (v_xx)')
    # axs[1].set_xlabel('Time (ms)')
    # axs[1].set_ylabel('PSP (mV)')
    # axs[1].plot(T*1000, np.transpose(v_xx)*1000)

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
    
    # # Third figure. Membrane potential
    # f3, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    
    # axs[0].set_title('Pyramidal Vm (selected cells)')
    # axs[0].set_xlabel('Time (ms)')
    # axs[0].set_ylabel('mV')
    # axs[0].plot(T*1000, np.transpose(Py_monitor.v[0:5])*1e3, lw=0.5, c=c_py)

    # axs[1].set_title('Interneurons Vm (selected cells)')
    # axs[1].set_xlabel('Time (ms)')
    # axs[1].set_ylabel('mV')
    # axs[1].plot(T*1000, np.transpose(In_monitor.v[0:5])*1e3, lw=0.5, c=c_inter)
    
    # f3.tight_layout()
   
    # f4, axs = plt.subplots(1, 1, figsize=(6,6))
    # axs.plot(np.transpose(v_pi) * 1000, np.transpose(v_ip) * 1000)
    # axs.set_xlabel('x1 = V_pi')
    # axs.set_ylabel('x3 = V_ip')
    
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
                    'p_PP': p_PP,
                    'p_II': p_II,
                    'R_py': r_P_rate,
                    'R_in': r_I_rate,
                    'R_b': r_B_rate,
                    'RECURRENT_PYRAMIDAL': RECURRENT_PYRAMIDAL,
                    'RECURRENT_INHIBITORY': RECURRENT_INHIBITORY,
                    'INHIBIT_INPUT': INHIBIT_INPUT,
                    'ACTIVE_INTERNEURONS': ACTIVE_INTERNEURONS,
                    'ACTIVE_GABAb': ACTIVE_GABAb,
                    'GAUSSIAN_REFRACTORY': GAUSSIAN_REFRACTORY,
                    'GAUSSIAN_THRESHOLD': GAUSSIAN_THRESHOLD,
                    'GABA_A_MULTIPLIER': GABA_A_MULTIPLIER,
                    'input_spike_rate': input_spike_rate,
                    'input_spike_rate_thalamic': input_spike_rate_thalamic,
                    'input_current': input_current} 

    
    # Save as lfp_last
    scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/lfp_last.mat',
                     mdict = save_dictionary)
    
    # Save at the end:
    i = 0
    while os.path.exists('C://Users/artemios/Documents/Multiscale_Models_Data/lfp_%s.mat' % i):
        i += 1
    
    save_str = format('lfp_%s.mat' %(i))
    scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/lfp_%s.mat' % (i),
                     mdict = save_dictionary)
    
    print(colored('Results of simulation saved as: ' + save_str, 'green'))

else:
    print(colored('Attention! Results of simulation were not saved. SAVE = False', 'yellow'))



# Run iteratively. Need to uncomment the def line at the start of the file.
# # a = np.arange(500, 501, 0.05)
# # b = np.arange(501, 510, 1)
# # ranges = np.concatenate((a, b))
# for iterations in ranges:
#     brunel(corriente=iterations)

# brunel(input_spike_rate = 3.7, SAVE = False, PLOT = True)