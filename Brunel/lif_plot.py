# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 15:39:00 2022

@author: Artemio Soto-Breceda [artemios]
"""

from brian2 import *
import lif_model_CUBN as cubn
import lif_model_COBN as cobn
import analyse_spike_trains as ast
import pyspike as spk

#%% Plot 
def plot_results(T, sp_P, sp_I, r_P, r_I, lfp_v, v_p, v_i, I_GABA, I_GABAb, I_AMPA, N, input_spike_rate, I_exc_py, I_inh_py, I_exc_in, I_inh_in, OS='local', folder_path='/data/gpfs/projects/punim0643/artemios/simulations/', save_str='last_results'):
    print('Plotting simulation results ...')
    # spike rates
    window_size = 10*ms#%100.1*ms # Size of the window for the smooth spike rate # 100.1 instead of 100 to avoid an annoying warning at the end of the simulation
    # window_size2 = 0.1*ms
    
    r_P_rate = r_P.smooth_rate(window='gaussian', width=window_size)
    if shape(r_P_rate) != shape(r_P.t):
        r_P_rate = r_P_rate[5:]
    
    r_I_rate = r_I.smooth_rate(window='gaussian', width=window_size)
    if shape(r_I_rate) != shape(r_I.t):
        r_I_rate = r_I_rate[5:]
    
    # r_S_rate = r_S.smooth_rate(window='gaussian', width=window_size)
    # if shape(r_S_rate) != shape(r_S.t):
    #     r_S_rate = r_S_rate[5:]
    
    f, axs = plt.subplots(4, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    
    # colors
    c_inter = 'C6'      # pink
    c_py    = 'C9'      # light blue
    c_b    = 'C0'      # blue
    c_sst   = 'C1'      # orange
    c_gray  = '#e0e0e0' # grey
    
    # raster
    axs[0].set_title('Raster N = {}, u = {} spikes/ms'.format(N, input_spike_rate))
    axs[0].set_ylabel('Neuron')
    axs[0].set_yticks([])
    axs[0].spines["top"].set_visible(False)
    axs[0].spines["right"].set_visible(False)
    
    axs[0].plot(sp_P.t / ms, sp_P.i + N, '.', markersize=2, label='Pyramidal', c=c_py)
    axs[0].plot(sp_I.t / ms, sp_I.i, '.', markersize=2, label='PV', c=c_inter)
    # axs[0].plot(sp_S.t / ms, sp_S.i, '.', markersize=2, label='SST', c=c_sst)
    axs[0].legend(loc=1)
            
    axs[1].set_title('Population rates, moving average')
    axs[1].set_ylabel('Spike rate (Hz)')
    axs[1].spines["top"].set_visible(False)
    axs[1].spines["right"].set_visible(False)   
        
    axs[1].plot(r_P.t / ms, r_P_rate / Hz, label='Pyramidal', c=c_py)
    axs[1].plot(r_I.t / ms, r_I_rate / Hz, label='PV', c=c_inter)
    # axs[1].plot(r_S.t / ms, r_S_rate / Hz, label='SST', c=c_sst)
    axs[1].legend(loc=1)
    
    # synaptic currents
    axs[2].set_title('Synaptic currents')
    axs[2].set_ylabel('Amplitude (A)')
    axs[2].set_xlabel('Time (ms)')
    axs[2].spines["top"].set_visible(False)
    axs[2].spines["right"].set_visible(False)
    # Others
    # axs[2].plot(T*1000, np.array(In_monitor.I_GABA_rec).transpose(), lw=0.5, c=c_gray) # , label='GABA (Inter)'
    # axs[2].plot(T*1000, np.array(Py_monitor.I_AMPA_rec).transpose(), lw=0.5, c=c_gray) # , label='AMPA (Py)'
    # # alphas
    axs[2].plot(T*1000, I_GABA, label='GABAa (Py)', c=c_inter)
    axs[2].plot(T*1000, I_GABAb, label='GABAb (Py)', c=c_py)
    axs[2].plot(T*1000, I_AMPA, label='AMPA (In)', c=c_b)
    # axs[2].plot(T*1000, np.array(Py_monitor.I_AMPA_cor).transpose(), lw=0.5, c=c_Cor , label='AMPA (Cortical)')
    axs[2].legend(loc=1)
    
    # LFP
    axs[3].set_title('LFP')
    axs[3].set_xlabel('Time (ms)')
    axs[3].set_ylabel('mV')
    axs[3].spines["top"].set_visible(False)
    axs[3].spines["right"].set_visible(False)
    axs[3].plot(T*1000, np.transpose(lfp_v), lw=0.5, label='LFP_V')
    # axs[3].plot(T*1000, np.transpose(lfp_), lw=0.5, label='LFP_I')
    axs[3].legend(loc=1)

    f.tight_layout() # Fixes the positions of subplots and labels

    # Second figure. PSP
    f2, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    
    axs[0].set_title('Pyramidal Vm (selected cells)')
    axs[0].set_xlabel('Time (ms)')
    axs[0].set_ylabel('mV')
    axs[0].plot(T*1000, v_p, lw=0.5, c=c_py)

    axs[1].set_title('Interneurons Vm (selected cells)')
    axs[1].set_xlabel('Time (ms)')
    axs[1].set_ylabel('mV')
    axs[1].plot(T*1000, v_i, lw=0.5, c=c_inter)
    
    # Second figure. PSP
    f3, axs = plt.subplots(2, 1, sharex=True, figsize=(10, 6.25)) # New figure with two subplots
    
    axs[0].set_title('Pyramidal synaptic input')
    axs[0].set_xlabel('Time (ms)')
    axs[0].set_ylabel('current (nA)')
    axs[0].plot(T*1000, I_exc_py / nA, lw=1, c=c_py, label='E')
    axs[0].plot(T*1000, I_inh_py / nA, lw=1, c=c_inter, label='I')
    axs[0].plot(T*1000, (I_exc_py + I_inh_py) / nA, lw=2, c=c_gray, label='E+I')
    axs[0].legend(loc=1)
    
    axs[1].set_title('Interneurons synaptic input')
    axs[1].set_xlabel('Time (ms)')
    axs[1].set_ylabel('current (nA)')
    axs[1].plot(T*1000, I_exc_in / nA, lw=1, c=c_py, label='E')
    axs[1].plot(T*1000, I_inh_in / nA, lw=1, c=c_inter, label='I')
    axs[1].plot(T*1000, (I_exc_in + I_inh_in) / nA, lw=2, c=c_gray, label='E+I')
    
    if OS=='local':
        plt.show() 
    else:
        plt.savefig(folder_path + save_str + '.png')
        print('Results saved as:' + save_str + '.png')
        plt.close('all')
    
def plot_spike_stats(sp_P, sp_I, t_start=0, OS='local'):
    #%% Plot ISI distance and synchronization 
    stp = ast.get_spike_trains(sp_P.t/second, sp_P.i, t_start=t_start) # Spike trains of Pyramidal neurons
    sti = ast.get_spike_trains(sp_I.t/second, sp_I.i, t_start=t_start) # Spike trains of inhibitory neurons
    
    # ISI distance
    f1 = plt.figure()
    isidist_py = 0 # ast.isi_distance(stp, population='Py', fig=f1)
    isidist_in = 0 # ast.isi_distance(sti, population='In', fig=f1)
    
    # SPIKE distance
    # f2 = plt.figure()
    spkdist_py = 0#ast.spike_distance(stp, population='Py', fig=f2)
    spkdist_in = 0#ast.spike_distance(sti, population='In', fig=f2)
    
    # Spike syncrhony
    si_py = ast.spike_si(stp, 'Py')
    si_in = ast.spike_si(sti, 'In')
    
    # ISI-CV
    cv_py, cvstd_py = ast.mean_isi_cv(stp, population='Py')
    cv_in, cvstd_in = ast.mean_isi_cv(sti, population='In')
    # cv_py = 0
    # cvstd_py = 0
    # cv_in = 0
    # cvstd_in = 0
    
    return cv_py, cvstd_py, cv_in, cvstd_in, si_py, si_in, spkdist_py, spkdist_in, isidist_py, isidist_in
