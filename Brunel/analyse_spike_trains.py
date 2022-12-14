# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 12:51:37 2022

@author: Artemio Soto-Breceda [artemios]

Uses PySpike to analyse the spiking output from lif_gaba.py or lif.py
"""

import os
import scipy.io
import pyspike as spk
import brian2 as b2
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

def main(sp_P, sp_I):
    
    stp = get_spike_trains(sp_P.t/b2.second, sp_P.i) # Spike trains of Pyramidal neurons
    sti = get_spike_trains(sp_I.t/b2.second, sp_I.i) # Spike trains of inhibitory neurons
    
    # merge spike trains
    merged_p = spk.merge_spike_trains(stp)
    merged_i = spk.merge_spike_trains(sti)
    
    # Define the same interval for both merged spike trains
    if merged_p.t_end > merged_i.t_end:
        merged_p.t_end = merged_i.t_end
    else:
        merged_i.t_end = merged_p.t_end
        
    # Get the distance between Py and In
    d = distance(merged_p, merged_i)
    
    return 0  

  
def get_spike_trains(t, i):
    
    idx, frequency = np.unique(i, return_counts=True)
    t_end = max(t)
    spike_trains = []
    
    for j in range(len(idx)):
        subset = t[i==idx[j]]
        spike_train = spk.SpikeTrain(subset, edges=t_end, is_sorted=True)
        spike_trains.append(spike_train)    
        
    return spike_trains
            


def distance(st1, st2):
    spike_profile = spk.spike_profile(st1, st2)
    x, y = spike_profile.get_plottable_data()
    
    plt.figure()
    plt.plot(x, y, '--k')
    print("SPIKE distance: %.8f" % spike_profile.avrg())
    plt.show()
    
    return spike_profile.avrg()
    
    
def distance_matrix(spike_trains, population=''):
    plt.figure()
    isi_distance = spk.isi_distance_matrix(spike_trains)
    plt.imshow(isi_distance, interpolation='none')
    plt.title("ISI-distance %s" %population)
    
    return isi_distance
    
    
    
    
# Spike synchronization index    
def spike_si(spike_trains, population=''):
    plt.figure()
    spike_sync = spk.spike_sync_matrix(spike_trains)
    plt.imshow(spike_sync, interpolation='none')
    plt.title("SPIKE-Sync %s" %population)
    
    plt.show()
    
    return spike_sync
    



# PSTH
def plot_psth(spike_trains, bins=0.01, fig=""):
    # TODO, this is weird. How does the funciton know when there is a stimulus?
    psth = spk.psth(spike_trains, bins)
    x, y = psth.get_plottable_data()
    
    # Check if fig is an input
    if not fig:
        plt.figure()
    else:
        plt.figure(fig)
    
    plt.plot(x,y)
    plt.show()
    
    return x, y

def isi_cv(st):
    # CV(isi) = std_/mean_isi
    if len(st) > 2:
        isi = np.diff(st.spikes)
        cv = np.std(isi)/np.mean(isi)
    else:
        cv = 0
    
    return cv

def mean_isi_cv(spike_trains):
    cv_vect = np.zeros(len(spike_trains))
    for i in range(0,len(spike_trains)-1):
        cv_vect[i] = isi_cv(spike_trains[i])
        
    return np.mean(cv_vect), np.std(cv_vect)
    
    






    
    
    