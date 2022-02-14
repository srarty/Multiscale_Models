# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 12:01:41 2022

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


#%% parameters  --------------------------------------------------------------
input_spike_rate = 1 # spikes/ms/cell # Threshold ~= 11 

simulation_time = 1 * second
dt_ = 100 * usecond
T = linspace(0, simulation_time, round(simulation_time/dt_)) # Time vector for plots (in seconds)

# populations
N_P = 1

# voltage
V_leak = -70. * mV      # Resting membrane potential
V_thr = -50 * mV        # Threshold
V_reset = -59 * mV # -59 * mV      # Reset voltage. Equal to V_leak-> To use Burkitt's, 2006 Eq. (12)

# membrane capacitance
C_m = 0.5 * nF

# membrane leak
g_m = 25. * nS

# membrane time constants
tau_m = C_m / g_m

# connectivity time constants
# Pyramidal
tau_d = 2 * ms      # Decay time constant (From Brunel and Wang 2001)
tau_r = 0.1 * ms    # Rising time constant (< 1 ms)
epsilon = 0.05*ms
tau_s = 0.5 * (tau_d + tau_r) + epsilon      # "Lumped" time constant for alpha function. 

# refractory period
tau_rp = 2 * ms

# Synaptic delay
tau_l = 1 * ms # 0.5 * ms # 0.5 * ms in Brunel and Wang 2001
delay_alpha = 0 * ms

# Cortical input
num_inputs = 1                    # Both thalamo-cortical and cortico-cortical 

# Synaptic efficacies
# AMPA (excitatory)
j = -21 * pA 

alpha_weight = 35.25 # To adjust the weight of each increase of 'x' vs 's' for alpha vs single exp, respectively # single_exp = 1 -> alpha = 35.25
single_exp_weight = 1 # 0.091608 # alpha = 1 -> single_exp = 0.091608
double_exp_weight = 10.413 # single_exp = 1 -> double_exp = 10.413 | theory: 1*ms**2/(tau_d * tau_r)
delayed_exp_weight = 20.8465 # single_exp = 1 -> diff_exp = 20.8465 | can be: (tau_m/tau_r), it is huge, tho.
alpha_simple_weight = 1.72 # theory = 1*ms**2/(tau_d + tau_r)**2 | single_exp = 1 -> alpha_simple_weight = 1.795 (when epsilon = 0)

# Alpha function's parameter (and double exponential) to fix the units in ds/dt
k = 1 / ms # Dimmensionless?, check Nicola and Campbell 2013

#%%
# model equations
eqs = '''
    dv / dt = (-v + V_leak - (I_tot/g_m)) / tau_m : volt (unless refractory)
    
    dv1 / dt = (-v1 - (I_AMPA1/g_m)) / tau_m : volt (unless refractory)
    dv2 / dt = (-v2 - (I_AMPA2/g_m)) / tau_m : volt (unless refractory)
    dv3 / dt = (-v3 - (I_AMPA3/g_m)) / tau_m : volt (unless refractory)
    dv4 / dt = (-v4 - (I_AMPA4/g_m)) / tau_m : volt (unless refractory)
    dv5 / dt = (-v5 - (I_AMPA5/g_m)) / tau_m : volt (unless refractory)
        
    I_tot = I_AMPA1 + I_AMPA2 + I_AMPA3 : amp
    
    I_AMPA1 = j * s_AMPA1 : amp
    ds_AMPA1 / dt = -s_AMPA1 / (tau_d + tau_r) : 1
    
    I_AMPA2 = j * s_AMPA2 : amp
    s_AMPA2 : 1
    
    I_AMPA3 = j * s_AMPA3 : amp
    s_AMPA3 : 1
    
    I_AMPA4 = j * s_AMPA4 : amp
    s_AMPA4 : 1
    
    I_AMPA5 = j * s_AMPA5 : amp
    s_AMPA5 : 1
'''

eqs_pre_ampa1 = '''
s_AMPA1 += single_exp_weight
''' 

eqs_ampa2 = '''
s_AMPA2_post = s_AMPA2_syn : 1 (summed)
ds_AMPA2_syn/dt = (x2 - s_AMPA2_syn) / tau_d : 1 (clock-driven)
dx2/dt = -x2/tau_r : 1 (clock-driven)
'''

eqs_pre_ampa2 = '''
x2 += delayed_exp_weight
''' 

eqs_ampa3 ='''
s_AMPA3_post = s_AMPA3_syn : 1 (summed)
ds_AMPA3_syn / dt = - s_AMPA3_syn / tau_d + k * x3 * (1 - s_AMPA3_syn) : 1 (clock-driven)
dx3 / dt = - x3 / tau_r : 1 (clock-driven)
'''

eqs_pre_ampa3 = '''
x3 += alpha_weight
''' 

eqs_ampa4 ='''
s_AMPA4_post = s_AMPA4_syn : 1 (summed)
ds_AMPA4_syn / dt = - s_AMPA4_syn / tau_r + (k * x4) : 1 (clock-driven)
dx4 / dt = -x4/tau_d: 1 (clock-driven)
'''

eqs_pre_ampa4 = '''
x4 += double_exp_weight
''' 

eqs_ampa5 ='''
s_AMPA5_post = s_AMPA5_syn : 1 (summed)
ds_AMPA5_syn / dt = - s_AMPA5_syn / (tau_s) + k * x5 : 1 (clock-driven)
dx5 / dt = - x5 / (tau_s) : 1 (clock-driven)
'''

eqs_pre_ampa5 = '''
x5 += alpha_simple_weight
''' 

Pyramidal = NeuronGroup(N_P, eqs, threshold='v > V_thr', reset='v = V_reset', refractory=tau_rp, method='rk4', dt=dt_, name='PyramidalPop') # Pyramidal population
Pyramidal.v = V_leak

# Input1 = PoissonInput(Pyramidal, 's_AMPA1', num_inputs, input_spike_rate * Hz, single_exp_weight)
# Input2 = PoissonInput(Pyramidal, 'x2', num_inputs, input_spike_rate * Hz, alpha_weight)

Input = PoissonGroup(num_inputs, rates=input_spike_rate * Hz, dt=dt_)

AMPA1_synapses = Synapses(Input, Pyramidal, on_pre=eqs_pre_ampa1, method='rk4', clock=Input.clock) # 
AMPA1_synapses.connect(p = 1)

# AMPA2_synapses = Synapses(Input, Pyramidal, model=eqs_ampa2, on_pre=eqs_pre_ampa2, method='rk4', clock=Input.clock, delay = tau_l)
# AMPA2_synapses.connect(p = 1)

# AMPA3_synapses = Synapses(Input, Pyramidal, model=eqs_ampa3, on_pre=eqs_pre_ampa3, method='rk4', dt=10*usecond) # When model eqs are used here, we need a smaller dt for resolution.
# AMPA3_synapses.connect(p = 1)

# AMPA4_synapses = Synapses(Input, Pyramidal, model=eqs_ampa4, on_pre=eqs_pre_ampa4, method='rk4', dt=10*usecond)
# AMPA4_synapses.connect(p = 1)

AMPA5_synapses = Synapses(Input, Pyramidal, model=eqs_ampa5, on_pre=eqs_pre_ampa5, method='rk4', dt=10*usecond, delay = delay_alpha)
AMPA5_synapses.connect(p = 1)

Py_monitor = StateMonitor(Pyramidal, ['v1', 'v2', 'v3', 'v4', 'v5', 's_AMPA1', 's_AMPA2', 's_AMPA3', 's_AMPA4', 's_AMPA5'], record = True) # Monitoring the AMPA and GABA currents in the Pyramidal population

#%% Run
net = Network(collect())
net.run(simulation_time, report='stdout')

#%% Plot
f, axs = plt.subplots(1, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    
axs.set_title('Pyramidal EPSP | weight: {}, j: {}'.format(alpha_simple_weight, j/pA))
axs.set_xlabel('Time (ms)')
axs.set_ylabel('mV')
axs.plot((np.transpose(Py_monitor.v1) * 1e3), lw=1, label='v1 (single_exp)')
# axs.plot((np.transpose(Py_monitor.v2) * 1e3), lw=1, label='v2 (delayed diff exp)')
# axs.plot((np.transpose(Py_monitor.v3) * 1e3), lw=1, label='v3 (alpha)')
# axs.plot((np.transpose(Py_monitor.v4) * 1e3), lw=1, label='v4 (double exp)')
axs.plot((np.transpose(Py_monitor.v5) * 1e3), lw=1, label='v5 (alpha NMM)')

f.legend()
# f.tight_layout()

f, axs = plt.subplots(1, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    
axs.set_xlabel('Time (ms)')
axs.plot(np.transpose(Py_monitor.s_AMPA1), lw=1, label='s_AMPA1 (single_exp)')
# axs.plot(np.transpose(Py_monitor.s_AMPA2), lw=1, label='s_AMPA2 (delayed diff exp)')
# axs.plot(np.transpose(Py_monitor.s_AMPA3), lw=1, label='s_AMPA3 (alpha)')
# axs.plot(np.transpose(Py_monitor.s_AMPA4), lw=1, label='s_AMPA4 (double exp)')
axs.plot(np.transpose(Py_monitor.s_AMPA5), lw=1, label='s_AMPA5 (alpha NMM)')

f.legend()
# f.tight_layout()

plt.show()