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

# def brunel(corriente=0):
# plt.close('all')
    
# Options:
RECURRENT_PYRAMIDAL = False     # Self excitation 
RECURRENT_INHIBITORY = False     # Self inhibition
INHIBIT_INPUT = True            # Excitatory cortical input to inhibitory population
ACTIVE_INTERNEURONS = False      # Inhibitory population
ACTIVE_SPINY = False            # Spiny Stellate population
PARAMS_SOURCE = 'allen'        # 'brunel' or 'allen'
SAVE = True                     # Save ground truth data
PLOT = False                     # Plot results (main Figure)
PLOT_EXTRA = False              # Plot extra things.

corriente = 100
# Balanced-rate network (?) with input currents: Py = 500.01 pA, In = 398 pA
input_current = corriente # corriente  # 437.5 # 500.01       # Injected current to Pyramidal population # Use this to calculate the nonlinearity (Vm -> Spike_rate sigmoid) on the disconnected model
input_current_I = corriente # 350 # 398 # 400.01     # Inhibitory interneurons
input_current_E = 0     # Excitatory interneurons (Spiny Stellate)         

input_spike_rate = 0 # spikes/ms/cell (driving input)
input_spike_rate_thalamic = 1.5 # spikes/ms/cell (spontaneous activity)

spiny_constant = 30 # temporal variable to  explore Spiny excitability

#%% parameters  --------------------------------------------------------------
simulation_time = 1 * second
dt_ = 100 * usecond
T = linspace(0, simulation_time, round(simulation_time/dt_)) # Time vector for plots (in seconds)
# T_u = linspace(0, simulation_time, round(simulation_time/u_dt)) # Time vector for u for plots (in seconds)
   
# populations
N = 1 # 135 # 675
N_P = 1#int(N*4)  # pyramidal neurons
N_E = 1#int(N)    # excitatory neurons (spiny stellate) 
N_I = 1#int(N)    # interneurons

# set populations parameters
params_py = set_params('pyramidal', PARAMS_SOURCE)
params_in = set_params('inhibitory', PARAMS_SOURCE)

# Probability of connection
p_IP = params_py.get('p_IP') # * np.sqrt(1000/N) #* 500/N #0.2 #* 100/N # Inhibitory to Pyramidal
p_PI = params_py.get('p_PI') # * np.sqrt(1000/N) # * 500/N #0.2 #* 100/N # Pyramidal to Inhibitory
p_PP = params_py.get('p_PP') # * np.sqrt(1000/N) #* 500/N  #0.2 #* 100/N # recurrent excitation (pyramidal) # Generally less than PI, IP connectivity (Bryson et al., 2021)
p_II = params_py.get('p_II') # * np.sqrt(1000/N) # * 500/N  #0.2 #* 100/N # recurrent inhibition
P_uP = 0.2

# voltage
V_leak = -70. * mV      # Resting membrane potential
V_thr = -50 * mV        # Threshold
V_reset = -59 * mV # -59 * mV      # Reset voltage. Equal to V_leak-> To use Burkitt's, 2006 Eq. (12)

# membrane capacitance
C_m_P = params_py.get('C')
C_m_I = params_in.get('C')

# membrane leak
g_m_P = params_py.get('g_leak')
g_m_I = params_in.get('g_leak')

# membrane time constants
tau_m_P = params_py.get('tau_m')
tau_m_I = params_in.get('tau_m')

# connectivity time constants
# Pyramidal
tau_d_AMPA_P = params_py.get('tau_AMPA_d') # Decay time constant (From Brunel and Wang 2001)
tau_r_AMPA_P = params_py.get('tau_AMPA_r')   # Rising time constant (< 1 ms), set to 0.05 isntead of 0.4 to match the ratio from GABA
tau_s_AMPA_P = tau_d_AMPA_P + tau_r_AMPA_P      # "Lumped" time constant for alpha function. 

tau_d_AMPA_P_ext = params_py.get('tau_AMPA_d_ext') # Decay time constant (From Brunel and Wang 2001)
tau_r_AMPA_P_ext = params_py.get('tau_AMPA_r_ext')   # Rising time constant (< 1 ms), set to 0.05 isntead of 0.4 to match the ratio from GABA
tau_s_AMPA_P_ext = tau_d_AMPA_P_ext + tau_r_AMPA_P_ext      # "Lumped" time constant for alpha function. 

tau_d_GABA_P = params_py.get('tau_GABA_d') # From Brunel and Wang 2001
tau_r_GABA_P = params_py.get('tau_GABA_r')
tau_s_GABA_P = tau_d_GABA_P + tau_r_GABA_P      # "Lumped" time constant for alpha function. 
#Inhibitory Interneurons
tau_d_AMPA_I = params_in.get('tau_AMPA_d')  
tau_r_AMPA_I = params_in.get('tau_AMPA_r')
tau_s_AMPA_I = tau_d_AMPA_I + tau_r_AMPA_I      # "Lumped" time constant for alpha function. 

tau_d_AMPA_I_ext = params_in.get('tau_AMPA_d_ext')  
tau_r_AMPA_I_ext = params_in.get('tau_AMPA_r_ext')
tau_s_AMPA_I_ext = tau_d_AMPA_I_ext + tau_r_AMPA_I_ext      # "Lumped" time constant for alpha function. 

tau_d_GABA_I = params_in.get('tau_GABA_d')
tau_r_GABA_I = params_in.get('tau_GABA_r')
tau_s_GABA_I = tau_d_GABA_I + tau_r_GABA_I      # "Lumped" time constant for alpha function. 

# refractory period
tau_rp_P = params_py.get('tau_rp')
tau_rp_I = params_in.get('tau_rp')

# Synaptic delay
delay = 0.2 * ms # 1 * ms # 0.5 * ms # 0.5 * ms in Brunel and Wang 2001

# Cortical input
num_inputs = 800                    # Both thalamo-cortical and cortico-cortical 

# Synaptic efficacies
# AMPA (excitatory)
j_AMPA_rec_P = params_py.get('j_AMPA') * 1000/N # * np.sqrt(1000/N)
j_AMPA_rec_I = params_in.get('j_AMPA') * 1000/N # * np.sqrt(1000/N)
    
j_AMPA_cor_P = params_py.get('j_AMPA_ext')
j_AMPA_cor_I = params_in.get('j_AMPA_ext')

j_AMPA_tha_P = params_py.get('j_AMPA_tha')
j_AMPA_tha_I = params_in.get('j_AMPA_tha')

# GABAergic (inhibitory)
j_GABA_P = params_py.get('j_GABA') * 1000/N # * np.sqrt(1000/N)
j_GABA_I = params_in.get('j_GABA') * 1000/N # * np.sqrt(1000/N) #alpha_ii * -1 * pA * 1000/N # 

# Weight constants. Amplitude of the synaptic input
# Pyramidal 
increment_AMPA_P =  params_py.get('alpha_weight_AMPA') #* np.sqrt(1000/N)
increment_AMPA_ext_P = params_py.get('single_exp') #* np.sqrt(1000/N)
increment_GABA_P = params_py.get('alpha_weight_GABA') #* np.sqrt(1000/N)

# Inhibitory interneurons
increment_AMPA_I = params_in.get('alpha_weight_AMPA') #* np.sqrt(1000/N)
increment_AMPA_ext_I = params_in.get('single_exp') #* np.sqrt(1000/N)
increment_GABA_I = params_in.get('alpha_weight_GABA') #* np.sqrt(1000/N)

# Alpha function's parameter (and double exponential) to fix the units in ds/dt
k = 1 / ms # 0.62 / ms # Dimmensionless?, check Nicola and Campbell 2013


# Injected current
I_injected = -input_current * pA # Input current to Pyramidal population. Sets a baseline spiking rate
I_injected_I = -input_current_I * pA # Input current to Pyramidal population. Sets a baseline spiking rate

# I_injected = np.random.normal(loc=-input_current, scale=10, size=int(l))
# I_injected[:] = I_injected
# I_injected = TimedArray(I_injected * pA, dt=dt_)

## Modeling
# model equations
eqs_P = get_equations('pyramidal')

eqs_I = get_equations('inhibitory')


# Neuron groups
Py_Pop = NeuronGroup(N_P, eqs_P, threshold='v > V_thr', reset='''v = V_reset
                                                                v_pe = V_reset-V_leak
                                                                v_pi = V_reset-V_leak
                                                                ''', refractory=tau_rp_P, method='rk4', dt=dt_, name='PyramidalPop') # Pyramidal population
Py_Pop.v = V_leak


In_Pop = NeuronGroup(N_I, eqs_I, threshold='v > V_thr', reset='''v = V_reset
                                                                v_ip = V_reset-V_leak
                                                                ''', refractory=tau_rp_I, method='rk4', dt=dt_, name='InhibitoryPop') # Interneuron population
In_Pop.v = V_leak

# Pop_Cor = PoissonGroup(num_inputs, rates = (input_spike_rate*1000/num_inputs)*Hz, dt=dt_) # poisson input

# sigma_sq_noise = 0.16 * volt * volt   # Ornstein-Uhlenbeck process (cortico-cortical) units need to be volt^2 for units in the population's equation to be consistent. Makes sense, because sigma should be in volts, hence sigma^2 in volts^2. Check cavalleri 2014
# tau_noise = 16 * ms                   # Ornstein-Uhlenbeck 
# input_cortical = NeuronGroup(num_inputs, 'dv/dt = -v/tau_noise + np.sqrt(2 * sigma_sq_noise/tau_noise) * xi : volt', threshold = 'v >= 0.5 * volt', reset = 'v = 0 * volt', dt=dt_, method='euler')

#%% synaptic equations -------------------------------------------------------
# Pyramidal (autoexcitatory)
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

# Interneurons (glutamate, pyramidal to inhibitory)
# eqs_glut_I = '''
# s_AMPA_post = s_AMPA_syn : 1 (summed)
# ds_AMPA_syn / dt = - s_AMPA_syn / tau_s_AMPA_I + k * x : 1 (clock-driven)
# dx / dt = - x / tau_s_AMPA_I : 1 (clock-driven)
# '''
eqs_pre_glut_I = '''
s_AMPA += increment_AMPA_I
'''

# Interneurons (recurrent inhibition)
# eqs_gaba_I = '''
# s_GABA_post = s_GABA_syn : 1 (summed)
# ds_GABA_syn / dt = - s_GABA_syn / tau_s_GABA_I + k * x : 1 (clock-driven)
# dx / dt = - x / tau_s_GABA_I : 1 (clock-driven)
# '''
eqs_pre_gaba_I = '''
s_GABA += increment_GABA_I
'''

   
# Synapses
# P to P
C_P_P = Synapses(Py_Pop, Py_Pop, on_pre=eqs_pre_glut_P, method='rk4', dt=dt_, delay=delay, name='synapses_pp')
C_P_P.connect('i != j', p = p_PP)
C_P_P.active = RECURRENT_PYRAMIDAL    # Testing no recursive connections to match NMM

# P to I
C_P_I = Synapses(Py_Pop, In_Pop, on_pre=eqs_pre_glut_I, method='rk4', dt=dt_, delay=delay, name='synapses_pi')
C_P_I.connect(p = p_PI)     
C_P_I.active = ACTIVE_INTERNEURONS

# I to I
C_I_I = Synapses(In_Pop, In_Pop, on_pre=eqs_pre_gaba_I, method='rk4', dt=dt_, delay=delay, name='synapses_ii')
C_I_I.connect('i != j', p = p_II)
C_I_I.active = RECURRENT_INHIBITORY    # Testing no recursive connections to match NMM

# I to P
C_I_P = Synapses(In_Pop, Py_Pop, on_pre=eqs_pre_gaba_P, method='rk4', dt=dt_, delay=delay, name='synapses_ip')
C_I_P.connect(p = p_IP)    
C_I_P.active = ACTIVE_INTERNEURONS
    
# External inputs
input1 =  PoissonInput(Py_Pop, 's_AMPA_cor', num_inputs, (input_spike_rate * 1000/num_inputs) * Hz, increment_AMPA_ext_P)

# Poisson input (Cortico-cortical) input to inhibitory interneurons. Controlled by INHIBIT_INPUT
C_Cor_I = PoissonInput(In_Pop, 's_AMPA_cor', num_inputs, (input_spike_rate * 1000/num_inputs) * Hz, increment_AMPA_ext_I)
C_Cor_I.active = INHIBIT_INPUT # Innactive cortico-cortical -> interneuron

# Poisson input (Thalamic, baseline spike rate)
C_Tha_P = PoissonInput(Py_Pop, 's_AMPA_tha', num_inputs, (input_spike_rate_thalamic*1000/num_inputs) * Hz, increment_AMPA_ext_P)
C_Tha_I = PoissonInput(In_Pop, 's_AMPA_tha', num_inputs, (input_spike_rate_thalamic*1000/num_inputs) * Hz, increment_AMPA_ext_I)


# Poisson population
# C_Cor_P = Synapses(Pop_Cor, Py_Pop, model=eqs_cor_P, on_pre=eqs_pre_cor_P, method='rk4', dt=dt_, delay=delay, name='synapses_pext')
# C_Cor_P.connect(p = 1)    
# C_Cor_I = Synapses(Pop_Cor, In_Pop, model=eqs_cor_I, on_pre=eqs_pre_cor_I, method='rk4', dt=dt_, delay=delay, name='synapses_iext')
# C_Cor_I.connect(p = 1)    

#%% monitors  -----------------------------------------------------------------
N_activity_plot = 30 # How many neurons in the raster plots (too large takes longer to monitor and plot)
sp_P = SpikeMonitor(Py_Pop[:]) # N_activity_plot])
sp_I = SpikeMonitor(In_Pop[:]) # N_activity_plot])
# sp_Cor = SpikeMonitor(input_cortical[:N_activity_plot])


r_P = PopulationRateMonitor(Py_Pop) # [0:N_activity_plot])
r_I = PopulationRateMonitor(In_Pop)
# r_Cor = PopulationRateMonitor(input_cortical)

st_AMPA_P = StateMonitor(Py_Pop, ('s_AMPA'), record = 0)
st_GABA_P = StateMonitor(Py_Pop, 's_GABA', record = 0)
st_AMPA_cor_P = StateMonitor(Py_Pop, 's_AMPA_cor', record = 0)

st_AMPA_I = StateMonitor(In_Pop, 's_AMPA', record = 0)
st_GABA_I = StateMonitor(In_Pop, 's_GABA', record = 0)
st_AMPA_cor_I = StateMonitor(In_Pop, 's_AMPA_cor', record = 0)

Py_monitor = StateMonitor(Py_Pop, ['I_AMPA_cor', 'I_AMPA_rec', 'I_GABA_rec', 'I_tot', 'v', 'v_pi'], record = True) # Monitoring the AMPA and GABA currents in the Pyramidal population
In_monitor = StateMonitor(In_Pop, ['v', 'v_ip', 'I_tot'], record = True)

#%% simulate  -----------------------------------------------------------------
net = Network(collect())
net.run(simulation_time, report='stdout') # Run first segment, if running more segments, run for a fraction of simulation_time

# # update prameteres here, then run the next segment:net.run(simulation_time/3, report='stdout') # Run second segment
# C_Cor_P.rate = 1
# C_Cor_I.rate = 1
# net.run(simulation_time/4, report='stdout') # Run second segment

#%% analysis ------------------------------------------------------------------
# spike rates
window_size = 100.1 * ms # Size of the window for the smooth spike rate # 100.1 instead of 100 to avoid an annoying warning at the end of the simulation

r_P_rate = r_P.smooth_rate(window='gaussian', width=window_size)
if shape(r_P_rate) != shape(r_P.t):
    r_P_rate = r_P_rate[5:]

r_I_rate = r_I.smooth_rate(window='gaussian', width=window_size)
if shape(r_I_rate) != shape(r_I.t):
    r_I_rate = r_I_rate[5:]
    
# r_Cor_rate = r_Cor.smooth_rate(width = window_size)
# if shape(r_Cor_rate) != shape(r_Cor.t):
#     r_Cor_rate = r_Cor_rate[1:]

# Calculate mean PSP (NMM states)
v_pi = mean(Py_monitor.v_pi, 0)
v_ip = mean(In_monitor.v_ip, 0)

# Generate LFP
# current based
# lfp = abs(Py_monitor.I_AMPA_cor) + abs(Py_monitor.I_AMPA_rec) + abs(Py_monitor.I_GABA_rec) # Absolute sum of currents
# lfp = (Py_monitor.I_GABA_rec) - (Py_monitor.I_AMPA_cor + Py_monitor.I_AMPA_rec + Py_monitor.I_AMPA_spi + I_injected) # Difference of currents
lfp = (Py_monitor.I_GABA_rec) - (Py_monitor.I_AMPA_cor + Py_monitor.I_AMPA_rec) # Difference of currents
lfp_ = sum(lfp,0) / g_m_P # Sum across all Pyramidal neurons and divide by the leak conductance to get volts

# voltage based
mean_v_Py = np.transpose(mean(Py_monitor.v,0) - V_leak) * 1e3
lfp_v = mean_v_Py/volt
       
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

    P_ = np.array(list(sp_P.spike_trains().values()))
    I_ = np.array(list(sp_I.spike_trains().values()))
    for i in range(0,shape(P_)[0]):
        P_[i] = P_[i]/second
        
    for i in range(0,shape(I_)[0]):
        I_[i] = I_[i]/second

    save_dictionary={'LFP': lfp_,
                    'LFP_V': lfp_v,
                    'lfp_dt' : dt_,
                    'Vm': -(I_injected/g_m_P), # To calculate the nonlinearity, need to simulate single cell disconnected network 
                    'Vm_interneurons': -(I_injected_I/g_m_I), # To calculate the nonlinearity, need to simulate single cell disconnected network 
                    'R_py': r_P_rate, # 1/diff(np.array(sp_P.t)).mean(),
                    'R_in': r_I_rate,
                    'I_in': In_monitor.I_tot,
                    'I_py': Py_monitor.I_tot,
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
    
    i = 0
    while os.path.exists('C://Users/artemios/Documents/Multiscale_Models_Data/nonlinearity/double_exp_v2/lfp_inputCurrent_%s_%s.mat' % (floor(input_current),i)):
        i += 1
    
    scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/nonlinearity/double_exp_v2/lfp_inputCurrent_%s_%s.mat' % (floor(input_current),i),
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
# a = np.arange(500.5, 510, 0.5)
# b = np.arange(-200, 610, 10)
# ranges = np.concatenate((a,b))
# ranges = np.arange(100, 610, 10)
# for iterations in ranges:
    # brunel(corriente = iterations)


# ranges = np.arange(0.1, 5.1, 0.1)
# for iterations in ranges:
#     brunel(input_spike_rate = iterations)
    
# brunel(input_spike_rate = 3.7, SAVE = False, PLOT = True)