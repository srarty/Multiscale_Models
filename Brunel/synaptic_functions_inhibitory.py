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

#%% options  --------------------------------------------------------------

neuron_type     = 'inhibitory'   # pyramidal, inhibitory or spiny
synaptic_type   = 'AMPA'        # AMPA or GABA
external        = False         # When AMPA, synapsis can be external or recurrent (local)
input_spike_rate = 1            # spikes/ms/cell 
simulation_time = 0.25 * second



#%% parameters  --------------------------------------------------------------
dt_ = 100 * usecond
T = linspace(0, simulation_time, round(simulation_time/dt_)) # Time vector for plots (in seconds)

# populations
N_P = 1

# voltage
V_leak = -70. * mV      # Resting membrane potential
V_thr = -50 * mV        # Threshold
V_reset = -59 * mV # -59 * mV      # Reset voltage. Equal to V_leak-> To use Burkitt's, 2006 Eq. (12)

params = set_params(neuron_type)

# membrane capacitance
C_m = params["C"]

# membrane leak
g_m = params["g_leak"]

# membrane time constants
tau_m = C_m / g_m

# connectivity time constants
# Interneurons
if synaptic_type == 'AMPA':
    tau_d =  params["tau_AMPA_d"]      # Decay time constant (From Brunel and Wang 2001)
    tau_r =  params["tau_AMPA_r"] # Rising time constant (< 1 ms)
else:
    tau_d =  params["tau_AMPA_d"]      # Decay time constant (From Brunel and Wang 2001)
    tau_r =  params["tau_AMPA_r"] # Rising time constant (< 1 ms)
    
tau_s = 0.5 * (tau_d + tau_r) + 0.05*ms      # "Lumped" time constant for alpha function. 
tau_l = 1 * ms
tau_rp =  params["tau_rp"] # refractory period



# Cortical input
num_inputs = 1                    # Both thalamo-cortical and cortico-cortical 
I_input = -300 * pA

# Synaptic efficacies
# AMPA (external on inhibitory interneurons)
if synaptic_type == 'AMPA':
    if external:
        j =  params["j_AMPA_ext"]
    else:
        j =  params["j_AMPA"]        
else:
    j =  params["j_GABA"]
        

single_exp_weight = 1 # Inverse: 1/1.52 = 0.6578947368421053
alpha_simple_weight = 2.85
delayed_exp_weight = 6.20 # Inverse 1.52/6.0995 = 0.24920075416017706

# Alpha function's parameter (and double exponential) to fix the units in ds/dt
k = 1 / ms # Dimmensionless?, check Nicola and Campbell 2013

#%%
# model equations
# Input
input_eqs = '''
dv / dt = (-v + V_leak - (I_input/g_m)) / tau_m : volt (unless refractory)
'''

# Population
if external:
    eqs = '''
        dv / dt = (-v + V_leak - (I_tot/g_m)) / tau_m : volt (unless refractory)
        
        dv1 / dt = (-v1 - (I_AMPA1/g_m)) / tau_m : volt (unless refractory)
        dv5 / dt = (-v5 - (I_AMPA5/g_m)) / tau_m : volt (unless refractory)
            
        I_tot = I_AMPA1 + I_AMPA5: amp
        
        I_AMPA1 = j * s_AMPA1 : amp
        ds_AMPA1 / dt = -s_AMPA1 / (tau_d + tau_r) : 1
        
        I_AMPA5 = j * s_AMPA5 : amp
        s_AMPA5 : 1
    '''
else:
    eqs = '''
        dv / dt = (-v + V_leak - (I_tot/g_m)) / tau_m : volt (unless refractory)
        
        dv1 / dt = (-v1 - (I_AMPA1/g_m)) / tau_m : volt (unless refractory)
        dv5 / dt = (-v5 - (I_AMPA5/g_m)) / tau_m : volt (unless refractory)
            
        I_tot = I_AMPA1 + I_AMPA5: amp
        
        I_AMPA1 = j * s_AMPA1 : amp
        s_AMPA1 : 1
        
        I_AMPA5 = j * s_AMPA5 : amp
        s_AMPA5 : 1
    '''

eqs_pre_ampa1 = '''
s_AMPA1 += single_exp_weight
''' 

eqs_ampa2 = '''
s_AMPA1_post = s_AMPA1_syn : 1 (summed)
ds_AMPA1_syn/dt = (x2 - s_AMPA1_syn) / tau_d : 1 (clock-driven)
dx2/dt = -x2/tau_r : 1 (clock-driven)
'''

eqs_pre_ampa2 = '''
x2 += delayed_exp_weight
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

if num_inputs > 1:
    Input = PoissonGroup(num_inputs, rates=input_spike_rate * Hz, dt=dt_)
else:
    Input = NeuronGroup(num_inputs, input_eqs, threshold='v > V_thr', reset='v = V_reset', refractory=tau_rp, method='rk4', dt=dt_) # Pyramidal population

if external:
    AMPA1_synapses = Synapses(Input, Pyramidal, on_pre=eqs_pre_ampa1, method='rk4', clock=Input.clock) # 
    AMPA1_synapses.connect(p = 1)
else:
    AMPA1_synapses = Synapses(Input, Pyramidal, model=eqs_ampa2, on_pre=eqs_pre_ampa2, method='rk4', clock=Input.clock, delay = tau_l)
    AMPA1_synapses.connect(p = 1)

AMPA5_synapses = Synapses(Input, Pyramidal, model=eqs_ampa5, on_pre=eqs_pre_ampa5, method='rk4', dt=10*usecond)
AMPA5_synapses.connect(p = 1)

Py_monitor = StateMonitor(Pyramidal, ['v1', 'v5', 's_AMPA1', 's_AMPA5'], record = True) # Monitoring the AMPA and GABA currents in the Pyramidal population

#%% Run
net = Network(collect())
net.run(simulation_time, report='stdout')

#%% Plot
f, axs = plt.subplots(1, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    
if external:
    axs.set_title('{} PSP | weight: {}, j: {}'.format(neuron_type, alpha_simple_weight, j/pA))
    axs.set_xlabel('Time (ms)')
    axs.set_ylabel('mV')
    axs.plot((np.transpose(Py_monitor.v1) * 1e3), lw=1, label='v1 (single_exp)')
else:
    axs.set_title('{} PSP | weight: {}, j: {}'.format(neuron_type, delayed_exp_weight, j/pA))
    axs.set_xlabel('Time (ms)')
    axs.set_ylabel('mV')
    axs.plot((np.transpose(Py_monitor.v1) * 1e3), lw=1, label='v1 (delayed diff exp)')
axs.plot((np.transpose(Py_monitor.v5) * 1e3), lw=1, label='v5 (alpha NMM)')

f.legend()
# f.tight_layout()

f, axs = plt.subplots(1, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    
axs.set_xlabel('Time (ms)')
if external:
    axs.plot(np.transpose(Py_monitor.s_AMPA1), lw=1, label='s_AMPA1 (single_exp)')
else:
    axs.plot(np.transpose(Py_monitor.s_AMPA1), lw=1, label='s_AMPA1 (delayed diff exp)')
axs.plot(np.transpose(Py_monitor.s_AMPA5), lw=1, label='s_AMPA5 (alpha NMM)')

f.legend()
# f.tight_layout()

plt.show()