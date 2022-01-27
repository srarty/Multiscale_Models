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
"""

"""
TODO: Calculate the PSP for each population and save it in the output file. Instead of the membrane potentials.
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

# Options:
RECURRENT_PYRAMIDAL = False     # Self excitation 
RECURRENT_INHIBITORY = False    # Self inhibition
INHIBIT_INPUT = False            # Excitatory cortical input to inhibitory population
ACTIVE_INTERNEURONS = True      # Inhibitory population
ACTIVE_SPINY = False            # Spiny Stellate population
SAVE = False                     # Save ground truth data
PLOT = True                    # Plot results

# Balanced-rate network (?) with input currents: Py = 500.01 pA, In = 398 pA
input_current = 500.01       # Injected current to Pyramidal population # Use this to calculate the nonlinearity (Vm -> Spike_rate sigmoid) on the disconnected model
input_current_I = 398 # 400.01     # Inhibitory interneurons
input_current_E = 0     # Excitatory interneurons (Spiny Stellate)         

input_spike_rate = 0 # spikes/ms/cell

#%% parameters  --------------------------------------------------------------
simulation_time = 1 * second
dt_ = 100 * usecond
T = linspace(0, simulation_time, round(simulation_time/dt_)) # Time vector for plots (in seconds)
# T_u = linspace(0, simulation_time, round(simulation_time/u_dt)) # Time vector for u for plots (in seconds)
   
# populations
N = 100 # 135 # 675
N_P = int(N*4)  # pyramidal neurons
N_E = int(N)    # excitatory neurons (spiny stellate) 
N_I = int(N)    # interneurons

# Probability of connection
p_IP = 0.2 * 100/N #* 135/N      # Inhibitory to Pyramidal
p_PI = 0.2 * 100/N #* 135/N      # Pyramidal to Inhibitory
p_PE =  0.2 * 100/N #* 135/N      # Pyramidal to Excitatory
p_EP =  0.2 * 100/N #* 135/N      # Excitatory to Pyramidal
p_PP =  0.2 * 100/N #* 135/N      # recurrent excitation (pyramidal) # Generally less than PI, IP connectivity (Bryson et al., 2021)
p_II =  0.2 * 100/N #* 135/N      # recurrent inhibition
p_input = 0.2 * 100/N #* 135/N 

# voltage
V_leak = -70. * mV      # Resting membrane potential
V_thr = -50 * mV        # Threshold
V_reset = -70 * mV # -59 * mV      # Reset voltage. Equal to V_leak-> To use Burkitt's, 2006 Eq. (12)

# membrane capacitance
C_m_P = 0.5 * nF
C_m_E = 0.2 * nF #0.5 * nF
C_m_I = 0.2 * nF

# membrane leak
g_m_P = 25. * nS
g_m_E = 10. * nS # 25. * nS
g_m_I = 20. * nS

# membrane time constants
tau_m_P = C_m_P / g_m_P
tau_m_E = C_m_E / g_m_E
tau_m_I = C_m_I / g_m_I

# connectivity time constants
# Pyramidal
tau_d_AMPA_P = 2 * ms       # Decay time constant (From Brunel and Wang 2001)
tau_r_AMPA_P = 0.05 * ms    # Rising time constant (< 1 ms), set to 0.05 isntead of 0.4 to match the ratio from GABA
tau_d_GABA_P = 10 * ms      # From Brunel and Wang 2001
tau_r_GABA_P =  0.25 * ms
# Spiny Stellate
tau_d_AMPA_E = 2 * ms # 2 * ms       # Decay time constant
tau_r_AMPA_E = 0.05 * ms # 0.4 * ms     # Rising time constant
#Inhibitory Interneurons
tau_d_AMPA_I = 2 * ms # 1 * ms       # Decay time constant
tau_r_AMPA_I = 0.05 * ms # 0.2 * ms     # Rising time constant
tau_d_GABA_I = 10 * ms
tau_r_GABA_I = 0.25 * ms

# refractory period
tau_rp_P = 2 * ms
tau_rp_E = 2 * ms
tau_rp_I = 1 * ms

# Synaptic delay
delay = None # 1 * ms # 0.5 * ms # 0.5 * ms in Brunel and Wang 2001


# Cortical input
num_inputs = 800                    # Both thalamo-cortical and cortico-cortical 


# Timed array -> input firing rate
l = simulation_time/dt_
segments = 5 # Number of segments of different spike rate
timed_rate = ones(int(l))
timed_rate[:] = input_spike_rate
## Uncomment the following lines for a timed array of rates, i.e. different input rates trhoughout a single simulation
# timed_rate[0 : int(l/segments)] = 5
# timed_rate[int(1*l/segments) : int(2*l/segments)] = 10
# timed_rate[int(2*l/segments) : int(3*l/segments)] = 20
# timed_rate[int(3*l/segments) : int(4*l/segments)] = 30
# timed_rate[int(4*l/segments) : int(5*l/segments)] = 40
# timed_rate[int(5*l/segments) : int(6*l/segments)] = 50
# timed_rate[int(6*l/segments) : int(7*l/segments)] = 100
# timed_rate[int(7*l/segments) : int(8*l/segments)] = 125
# timed_rate[int(8*l/segments) : int(9*l/segments)] = 150
# timed_rate[int(9*l/segments) :] = 200
timed_rate = timed_rate * Hz
timed_rate = TimedArray(timed_rate, dt=dt_)


# Synaptic efficacies
# AMPA (excitatory)
j_AMPA_rec_P = -10.5 * pA #* 100/N # * p_PP * 4000 / N_P
j_AMPA_rec_P_spi = -10.5 * pA # * 100/N
j_AMPA_rec_E = -10.5 * pA #* 100/N
j_AMPA_rec_I = -10.5 * pA # * 100/N # -14 * pA

    
j_AMPA_cor_P = -10.5 * pA # -13.75 * pA
j_AMPA_cor_I = -10.5 * pA # -19 * pA

# GABAergic (inhibitory)
j_GABA_P = 71.4 * pA # * 100/N # 54 * pA # 42.5 * pA # * p_IP * 1000 / N_I #* p_IP * (N/135)
j_GABA_I = 71.4 * pA # * 100/N # 54 * pA # * p_II * 1000 / N_I #* p_II * (N/135)

# Alpha function's parameter
alpha = 1 / ms # 0.62 / ms # Dimmensionless?, check Nicola and Campbell 2013
alpha_pi = 1 * alpha
alpha_ip = 1 * alpha

# Injected current
I_injected = -input_current * pA # Input current to Pyramidal population. Sets a baseline spiking rate
I_injected_I = -input_current_I * pA # Input current to Pyramidal population. Sets a baseline spiking rate
I_injected_E = -input_current_E * pA # Input current to Pyramidal population. Sets a baseline spiking rate

# I_injected = np.random.normal(loc=-input_current, scale=10, size=int(l))
# I_injected[:] = I_injected
# I_injected = TimedArray(I_injected * pA, dt=dt_)

#%% modeling  ----------------------------------------------------------------
# model equations
eqs_P = '''
    dv / dt = (-v + V_leak - (I_tot/g_m_P)) / tau_m_P : volt (unless refractory)
    
    I_tot = I_AMPA_cor + I_AMPA_rec + I_AMPA_spi + I_GABA_rec + I_injected : amp
    
    I_AMPA_cor = j_AMPA_cor_P * s_AMPA_cor : amp
    ds_AMPA_cor / dt = - s_AMPA_cor / (tau_d_AMPA_P + tau_r_AMPA_P) : 1
    
    I_GABA_rec = j_GABA_P * s_GABA : amp
    s_GABA : 1
    
    I_AMPA_rec = j_AMPA_rec_P * s_AMPA : amp
    s_AMPA : 1    
    
    I_AMPA_spi = j_AMPA_rec_P_spi * s_AMPA_spi : amp
    s_AMPA_spi : 1
'''
    # dv_e/dt = (-v + V_leak - ((I_AMPA_spi + I_AMPA_cor)/g_m_P)) / tau_m_P : volt (unless refractory)
    # dv_i/dt = (-v + V_leak - (I_GABA_rec/g_m_P)) / tau_m_P : volt (unless refractory)
    

eqs_E = '''
    dv / dt = (-v + V_leak - (I_tot/g_m_E)) / tau_m_E : volt (unless refractory)
    
    I_tot = I_AMPA_rec + I_injected_E : amp        
    
    I_AMPA_rec = j_AMPA_rec_E * s_AMPA : amp
    s_AMPA : 1
'''

eqs_I = '''
    dv / dt = (-v + V_leak - (I_tot/g_m_I)) / tau_m_I : volt (unless refractory)
    
    I_tot = I_AMPA_cor + I_AMPA_rec + I_GABA_rec + I_injected_I : amp
    
    I_AMPA_cor = j_AMPA_cor_I * s_AMPA_cor : amp
    ds_AMPA_cor / dt = - s_AMPA_cor / (tau_d_AMPA_I + tau_r_AMPA_I) : 1
    
    I_GABA_rec = j_GABA_I * s_GABA : amp
    s_GABA : 1
    
    I_AMPA_rec = j_AMPA_rec_I * s_AMPA : amp
    s_AMPA : 1
'''

# Neuron groups
P_P = NeuronGroup(N_P, eqs_P, threshold='v > V_thr', reset='v = V_reset', refractory=tau_rp_P, method='rk4', dt=dt_) # Pyramidal population
P_P.v = V_leak

P_E = NeuronGroup(N_E, eqs_E, threshold='v > V_thr', reset='v = V_reset', refractory=tau_rp_E, method='rk4', dt=dt_) # Pyramidal population
P_E.v = V_leak

P_I = NeuronGroup(N_I, eqs_I, threshold='v > V_thr', reset='v = V_reset', refractory=tau_rp_I, method='rk4', dt=dt_) # Interneuron population
P_I.v = V_leak

# sigma_sq_noise = 0.16 * volt * volt   # Ornstein-Uhlenbeck process (cortico-cortical) units need to be volt^2 for units in the population's equation to be consistent. Makes sense, because sigma should be in volts, hence sigma^2 in volts^2. Check cavalleri 2014
# tau_noise = 16 * ms                   # Ornstein-Uhlenbeck 
# input_cortical = NeuronGroup(num_inputs, 'dv/dt = -v/tau_noise + sqrt(2 * sigma_sq_noise/tau_noise) * xi : volt', threshold = 'v >= 0.5 * volt', reset = 'v = 0 * volt', dt=dt_, method='euler')

#%% synaptic equations -------------------------------------------------------
# Pyramidal (autoexcitatory)
eqs_glut_P = '''
s_AMPA_post = s_AMPA_syn : 1 (summed)
ds_AMPA_syn / dt = - s_AMPA_syn / tau_d_AMPA_P + alpha * x * (1 - s_AMPA_syn) : 1 (clock-driven)
dx / dt = - x / tau_r_AMPA_P : 1 (clock-driven)
'''    

eqs_pre_glut_P = '''
x += 1
'''

# Pyramidal (AMPA from Spiny Stellate to Pyramidal)
eqs_spi_P = '''
s_AMPA_spi_post = s_AMPA_spi_syn : 1 (summed)
ds_AMPA_spi_syn / dt = - s_AMPA_spi_syn / tau_d_AMPA_P + alpha * x * (1 - s_AMPA_spi_syn) : 1 (clock-driven)
dx / dt = - x / tau_r_AMPA_P : 1 (clock-driven)
'''    

eqs_pre_spi_P = '''
x += 1
'''

# Pyramidal (gabaergic, inhibitory to pyramidal)
eqs_gaba_P = '''
s_GABA_post = s_GABA_syn : 1 (summed)
ds_GABA_syn / dt = - s_GABA_syn / tau_d_GABA_P + alpha_ip * x * (1 - s_GABA_syn) : 1 (clock-driven)
dx / dt = - x / tau_r_GABA_P : 1 (clock-driven)
'''

eqs_pre_gaba_P = '''
x += 1
'''

# Pyramidal (external input)
eqs_pre_cor_P = '''
s_AMPA_cor += 1
'''

# Excitatory (glutamate, pyramidal to spiny)
eqs_glut_E = '''
s_AMPA_post = s_AMPA_syn : 1(summed)
ds_AMPA_syn / dt = - s_AMPA_syn / tau_d_AMPA_E + alpha * x * (1 - s_AMPA_syn) : 1 (clock-driven)
dx / dt = - x / tau_r_AMPA_E : 1 (clock-driven)
'''

eqs_pre_glut_E = '''
x += 1
'''

# Interneurons (glutamate, pyramidal to inhibitory)
eqs_glut_I = '''
s_AMPA_post = s_AMPA_syn : 1 (summed)
ds_AMPA_syn / dt = - s_AMPA_syn / tau_d_AMPA_I + alpha_pi * x * (1 - s_AMPA_syn) : 1 (clock-driven)
dx / dt = - x / tau_r_AMPA_I : 1 (clock-driven)
''' 

eqs_pre_glut_I = '''
x += 1
'''

# Interneurons (autoinhibiton)
eqs_gaba_I = '''
s_GABA_post = s_GABA_syn : 1 (summed)
ds_GABA_syn / dt = - s_GABA_syn / tau_d_GABA_I + alpha * x * (1 - s_GABA_syn) : 1 (clock-driven)
dx / dt = - x / tau_r_GABA_I : 1 (clock-driven)
'''

eqs_pre_gaba_I = '''
x += 1
'''

# Interneurons (external excitatory input)
eqs_pre_cor_I = '''
s_AMPA_cor += 1
''' 

   
# Synapses
# P to P
C_P_P = Synapses(P_P, P_P, model=eqs_glut_P, on_pre=eqs_pre_glut_P, method='rk4', dt=dt_, delay=delay)
C_P_P.connect('i != j', p = p_PP)
C_P_P.active = RECURRENT_PYRAMIDAL    # Testing no recursive connections to match NMM

# P to E
C_P_E = Synapses(P_P, P_E, model=eqs_glut_E, on_pre=eqs_pre_glut_E, method='rk4', dt=dt_, delay=delay)
C_P_E.connect(p = p_PE)
C_P_E.active = ACTIVE_SPINY   
    
# P to I
C_P_I = Synapses(P_P, P_I, model=eqs_glut_I , on_pre=eqs_pre_glut_I, method='rk4', dt=dt_, delay=delay)
C_P_I.connect(p = p_PI)     
C_P_I.active = ACTIVE_INTERNEURONS

# I to I
C_I_I = Synapses(P_I, P_I, model=eqs_gaba_I, on_pre=eqs_pre_gaba_I, method='rk4', dt=dt_, delay=delay)
C_I_I.connect('i != j', p = p_II)
C_I_I.active = RECURRENT_INHIBITORY    # Testing no recursive connections to match NMM

# I to P
C_I_P = Synapses(P_I, P_P, model=eqs_gaba_P, on_pre=eqs_pre_gaba_P, method='rk4', dt=dt_, delay=delay)
C_I_P.connect(p = p_IP)    
C_I_P.active = ACTIVE_INTERNEURONS

# E to P
C_E_P = Synapses(P_E, P_P, model=eqs_spi_P, on_pre=eqs_pre_spi_P, method='rk4', dt=dt_, delay=delay)
C_E_P.connect(p = p_EP)    
C_E_P.active = ACTIVE_SPINY

# external input
C_Cor_P = PoissonInput(P_P, 's_AMPA_cor', num_inputs, (input_spike_rate*1000/num_inputs) * Hz, '1')
C_Cor_I = PoissonInput(P_I, 's_AMPA_cor', num_inputs, (input_spike_rate*1000/num_inputs) * Hz, '1')
C_Cor_I.active = INHIBIT_INPUT # Innactive cortico-cortical -> interneuron

#%% monitors  -----------------------------------------------------------------
N_activity_plot = 30 # How many neurons in the raster plots (too large takes longer to monitor and plot)
sp_P = SpikeMonitor(P_P[:N_activity_plot])
sp_E = SpikeMonitor(P_E[:N_activity_plot])
sp_I = SpikeMonitor(P_I[:N_activity_plot])
# sp_Cor = SpikeMonitor(input_cortical[:N_activity_plot])


r_P = PopulationRateMonitor(P_P) # [0:N_activity_plot])
r_E = PopulationRateMonitor(P_E) # [0:N_activity_plot])
r_I = PopulationRateMonitor(P_I)
# r_Cor = PopulationRateMonitor(input_cortical)

st_AMPA_P = StateMonitor(P_P, ('s_AMPA'), record = 0)
st_GABA_P = StateMonitor(P_P, 's_GABA', record = 0)
st_AMPA_cor_P = StateMonitor(P_P, 's_AMPA_cor', record = 0)

st_AMPA_E = StateMonitor(P_E, ('s_AMPA'), record = 0)

st_AMPA_I = StateMonitor(P_I, 's_AMPA', record = 0)
st_GABA_I = StateMonitor(P_I, 's_GABA', record = 0)
st_AMPA_cor_I = StateMonitor(P_I, 's_AMPA_cor', record = 0)

# Py_monitor = StateMonitor(P_P, ['I_AMPA_cor', 'I_AMPA_rec', 'I_GABA_rec', 'I_AMPA_spi', 'I_tot', 'v', 'v_e', 'v_i'], record = True) # Monitoring the AMPA and GABA currents in the Pyramidal population
Py_monitor = StateMonitor(P_P, ['I_AMPA_cor', 'I_AMPA_rec', 'I_GABA_rec', 'I_AMPA_spi', 'I_tot', 'v'], record = True) # Monitoring the AMPA and GABA currents in the Pyramidal population
In_monitor = StateMonitor(P_I, ['v'], record = True)
Ex_monitor = StateMonitor(P_E, ['v'], record = True)

#%% simulate  -----------------------------------------------------------------
net = Network(collect())
net.run(simulation_time, report='stdout') # Run first segment, if running more segments, run for a fraction of simulation_time

# update prameteres here, then run the next segment:
# j_GABA_P = j_GABA_P * 5
# I_inj_spiny = I_inj_spiny * 2
# net.run(simulation_time/3, report='stdout') # Run second segment

# # j_GABA_P = j_GABA_P / 5
# I_inj_spiny = I_inj_spiny * 2
# net.run(simulation_time/3, report='stdout') # Run second segment


#%% analysis ------------------------------------------------------------------
# Generate LFP
# current based
# lfp = abs(Py_monitor.I_AMPA_cor) + abs(Py_monitor.I_AMPA_rec) + abs(Py_monitor.I_GABA_rec) # Absolute sum of currents
# lfp = (Py_monitor.I_GABA_rec) - (Py_monitor.I_AMPA_cor + Py_monitor.I_AMPA_rec + Py_monitor.I_AMPA_spi + I_injected) # Difference of currents
lfp = (Py_monitor.I_GABA_rec) - (Py_monitor.I_AMPA_cor + Py_monitor.I_AMPA_rec + Py_monitor.I_AMPA_spi) # Difference of currents
lfp_ = sum(lfp,0) / g_m_P # Sum across all Pyramidal neurons and divide by the leak conductance to get volts

# voltage based
mean_v_Py = np.transpose(mean(Py_monitor.v,0) - V_leak) * 1e3
lfp_v = mean_v_Py/volt


# spike rates
window_size = 100.1 * ms # Size of the window for the smooth spike rate # 100.1 instead of 100 to avoid an annoying warning at the end of the simulation

r_P_rate = r_P.smooth_rate(window='flat', width=window_size)
if shape(r_P_rate) != shape(r_P.t):
    r_P_rate = r_P_rate[1:]

r_E_rate = r_E.smooth_rate(window='flat', width=window_size)
if shape(r_E_rate) != shape(r_E.t):
    r_E_rate = r_E_rate[1:]
    
r_I_rate = r_I.smooth_rate(window='flat', width=window_size)
if shape(r_I_rate) != shape(r_I.t):
    r_I_rate = r_I_rate[1:]
    
# r_Cor_rate = r_Cor.smooth_rate(width = window_size)
# if shape(r_Cor_rate) != shape(r_Cor.t):
#     r_Cor_rate = r_Cor_rate[1:]

# Calculate PSP (NMM states)
Py_EPSP = mean((Py_monitor.I_AMPA_cor + Py_monitor.I_AMPA_rec + Py_monitor.I_AMPA_spi)/g_m_P, 0)
Py_IPSP = mean(Py_monitor.I_GABA_rec/g_m_P, 0)

       
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
    
    axs[0].plot(sp_E.t / ms, sp_E.i + 2 * N_activity_plot, '.', markersize=2, label='Spiny', c=c_ex)
    axs[0].plot(sp_P.t / ms, sp_P.i + 1 * N_activity_plot, '.', markersize=2, label='Pyramidal', c=c_py)
    axs[0].plot(sp_I.t / ms, sp_I.i, '.', markersize=2, label='Inhibitory', c=c_inter)
    axs[0].legend(loc=1)
            
    axs[1].set_title('Population rates, moving average')
    axs[1].set_ylabel('Spike rate (Hz)')
    axs[1].spines["top"].set_visible(False)
    axs[1].spines["right"].set_visible(False)   
        
    axs[1].plot(r_P.t / ms, r_P_rate / Hz, label='Pyramidal', c=c_py)
    axs[1].plot(r_E.t / ms, r_E_rate / Hz, label='Excitatory', c=c_ex)
    axs[1].plot(r_I.t / ms, r_I_rate / Hz, label='Interneuron', c=c_inter)
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
    axs[2].plot(T*1000, np.array(st_AMPA_E.s_AMPA).transpose(), label='AMPA (Ex)', c=c_ex)
    axs[2].plot(T*1000, np.array(st_AMPA_I.s_AMPA).transpose(), label='AMPA (In)', c=c_inter)
    axs[2].plot(T*1000, np.array(st_AMPA_cor_P.s_AMPA_cor).transpose(), lw=0.5, c=c_Cor , label='AMPA_cor (Cortical)')
    axs[2].legend(loc=1)
    
    # LFP
    axs[3].set_title('LFP')
    axs[3].set_xlabel('Time (ms)')
    axs[3].set_ylabel('mV')
    axs[3].spines["top"].set_visible(False)
    axs[3].spines["right"].set_visible(False)
    axs[3].plot(T*1000, np.transpose(lfp_v), label='LFP_V')
    axs[3].legend(loc=1)

    f.tight_layout() # Fixes the positions of subplots and labels
        
    # Second figure. LFPs
    f2, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    
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
    if ACTIVE_SPINY:
        f3, axs = plt.subplots(3, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    else:
        f3, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    
    axs[0].set_title('Pyramidal Vm (selected cells)')
    axs[0].set_xlabel('Time (ms)')
    axs[0].set_ylabel('mV')
    axs[0].plot(T*1000, np.transpose(Py_monitor.v[0:5])*1e3, lw=0.5, c=c_py)

    axs[1].set_title('Interneurons Vm (selected cells)')
    axs[1].set_xlabel('Time (ms)')
    axs[1].set_ylabel('mV')
    axs[1].plot(T*1000, np.transpose(In_monitor.v[0:5])*1e3, lw=0.5, c=c_inter)
    
    if ACTIVE_SPINY:
        axs[2].set_title('Spiny Vm (selected cells)')
        axs[2].set_xlabel('Time (ms)')
        axs[2].set_ylabel('mV')
        axs[2].plot(T*1000, np.transpose(Ex_monitor.v[0:5])*1e3, lw=0.5, c=c_ex)
    
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
    i = 0
    while os.path.exists('C://Users/artemios/Documents/Multiscale_Models_Data/lfp_%s.mat' % i):
        i += 1
        
    alpha_1 = j_GABA_P * N_I * p_IP /pA
    alpha_2 = j_AMPA_rec_I * N_P * p_PI /pA
    
    E_ = np.array(list(sp_E.spike_trains().values()))
    I_ = np.array(list(sp_I.spike_trains().values()))
    for i in range(0,shape(E_)[0]):
        E_[i] = E_[i]/second
        
    for i in range(0,shape(I_)[0]):
        I_[i] = I_[i]/second
    
    save_dictionary={'LFP': lfp_,
                    'LFP_V': lfp_v,
                    'lfp_dt' : dt_,
                    'Vm': -(I_injected/g_m_P), # To calculate the nonlinearity, need to simulate single cell disconnected network 
                    'R_py': r_P_rate, # 1/diff(np.array(sp_P.t)).mean(),
                    'R_ex': r_E_rate,
                    'R_in': r_I_rate,
                    'alpha1': alpha_1,
                    'alpha2': alpha_2,
                    'v_rest': V_leak,
                    'isi_P': isi_P,
                    'V_py': Py_monitor.v,
                    'V_in': In_monitor.v,
                    'RECURRENT_PYRAMIDAL': RECURRENT_PYRAMIDAL,
                    'RECURRENT_INHIBITORY': RECURRENT_INHIBITORY,
                    'INHIBIT_INPUT': INHIBIT_INPUT,
                    'ACTIVE_INTERNEURONS': ACTIVE_INTERNEURONS,
                    'ACTIVE_SPINY': ACTIVE_SPINY,
                    'input_spike_rate': input_spike_rate,
                    'input_current': input_current}   
        
    # Save as lfp_last
    scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/lfp_last.mat',
                     mdict = save_dictionary)
    
    # scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/lfp_%s.mat' % i,
    #                  mdict = save_dictionary)
    
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
# a = np.arange(500, 501, 0.05)
# b = np.arange(501, 510, 1)
# c = np.arange(510, 700, 10)
# d = np.arange(0, 1600, 100)
# ranges = np.concatenate((a, b, c, d))

# ranges = np.arange(0.1, 5.1, 0.1)
# for iterations in ranges:
#     brunel(input_spike_rate = iterations)
    
# brunel(input_spike_rate = 3.7, SAVE = False, PLOT = True)