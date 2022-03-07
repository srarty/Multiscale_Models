# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 15:05:01 2022

@author: Artemio Soto-Breceda [artemios]

Loads saved data and plots it in the same way lif.py plots the output.
"""
        
import scipy.io
import numpy as np
from brian2 import *
from scipy import signal
import matplotlib.pyplot as plt
prefs.codegen.target = 'cython' #'numpy'  # use the Python fallback instead of C compilation
devices.device.shape = []       # This and the following line remove an annoying warning when brian2 is imported and loaded into RAM
devices.device.size = []

#%% load data -----------------------------------------------------------------
# folder = "C://Users/artemios/Documents/Multiscale_Models_Data/"
folder = "C://Users/artemios/Documents/Multiscale_Models_Data/Spartan/"
file = "lfp_13.mat"

data = scipy.io.loadmat(folder + file)

# Data
r_P_rate =  data.get('R_py')
r_E_rate =  data.get('R_ex')
r_I_rate =  data.get('R_in')
lfp_v =     data.get('LFP_V')
lfp_ =      data.get('LFP')
v_pi =      data.get('v_pi')
v_ip =      data.get('v_ip')
V_py =      data.get('V_py')
V_in =      data.get('V_in')
dt = data.get('lfp_dt').item()

# Flag(s)
ACTIVE_SPINY = data.get('ACTIVE_SPINY')

# Time array
T = linspace(0, size(r_P_rate)*dt, size(r_P_rate))

#%% plotting  -----------------------------------------------------------------
# plt.close('all')

f, axs = plt.subplots(4, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots

# colors
c_inter = 'C6'  # pink
c_py = 'C9'     # light blue
c_ex = 'C0'    # blue
c_Cor = 'C1'    # orange
c_gray = '#e0e0e0' # grey

# raster
# axs[0].set_title('Raster ({} neurons/pop) N = {}, u = {} spikes/ms'.format(N_activity_plot, N, input_spike_rate))
# axs[0].set_ylabel('Neuron')
# axs[0].set_yticks([])
# axs[0].spines["top"].set_visible(False)
# axs[0].spines["right"].set_visible(False)

# axs[0].plot(sp_E.t / ms, sp_E.i + 2 * N_activity_plot, '.', markersize=2, label='Spiny', c=c_ex)
# axs[0].plot(sp_P.t / ms, sp_P.i + 1 * N_I, '.', markersize=2, label='Pyramidal', c=c_py)
# axs[0].plot(sp_I.t / ms, sp_I.i, '.', markersize=2, label='Inhibitory', c=c_inter)
# axs[0].legend(loc=1)
        
axs[1].set_title('Population rates, moving average')
axs[1].set_ylabel('Spike rate (Hz)')
axs[1].spines["top"].set_visible(False)
axs[1].spines["right"].set_visible(False)   
    
axs[1].plot(T * 1000, transpose(r_P_rate) / Hz, label='Pyramidal', c=c_py)
axs[1].plot(T * 1000, transpose(r_E_rate) / Hz, label='Excitatory', c=c_ex)
axs[1].plot(T * 1000, transpose(r_I_rate) / Hz, label='Interneuron', c=c_inter)
# axs[1].plot(r_Cor.t / ms, r_Cor_rate / Hz, label='Cortico-cortical (OU)', c=c_Cor)
# axs[1].plot(T_u / ms, u / Hz)
axs[1].legend(loc=1)

# # synaptic currents
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
# axs[2].plot(T*1000, np.array(st_AMPA_E.s_AMPA).transpose(), label='AMPA (Ex)', c=c_ex)
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
axs[3].plot(T*1000, np.transpose(lfp_), lw=0.5, label='LFP_I')
axs[3].legend(loc=1)

f.tight_layout() # Fixes the positions of subplots and labels

#%% Second figure. PSP
f2, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots

axs[0].set_title('Pyramidal population (v_pi)')
axs[0].set_xlabel('Time (ms)')
axs[0].set_ylabel('IPSP (mV)')
axs[0].plot(T*1000, np.transpose(mean(v_pi,0))*1000)

axs[1].set_title('Inhibitory population (v_ip)')
axs[1].set_xlabel('Time (ms)')
axs[1].set_ylabel('EPSP (mV)')
axs[1].plot(T*1000, np.transpose(mean(v_ip,0))*1000)

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
axs[0].plot(T*1000, np.transpose(V_py[0:5])*1e3, lw=0.5, c=c_py)

axs[1].set_title('Interneurons Vm (selected cells)')
axs[1].set_xlabel('Time (ms)')
axs[1].set_ylabel('mV')
axs[1].plot(T*1000, np.transpose(V_in[0:5])*1e3, lw=0.5, c=c_inter)

if ACTIVE_SPINY:
    axs[2].set_title('Spiny Vm (selected cells)')
    axs[2].set_xlabel('Time (ms)')
    axs[2].set_ylabel('mV')
    axs[2].plot(T*1000, np.transpose(Ex_monitor.v[0:5])*1e3, lw=0.5, c=c_ex)

f3.tight_layout()
plt.show()
    