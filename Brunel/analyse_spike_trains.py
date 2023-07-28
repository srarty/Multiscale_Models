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

# def main(sp_P, sp_I):
    
#     stp = get_spike_trains(sp_P.t/b2.second, sp_P.i) # Spike trains of Pyramidal neurons
#     sti = get_spike_trains(sp_I.t/b2.second, sp_I.i) # Spike trains of inhibitory neurons
    
#     # merge spike trains
#     merged_p = spk.merge_spike_trains(stp)
#     merged_i = spk.merge_spike_trains(sti)
    
#     # Define the same interval for both merged spike trains
#     if merged_p.t_end > merged_i.t_end:
#         merged_p.t_end = merged_i.t_end
#     else:
#         merged_i.t_end = merged_p.t_end
        
#     # Get the distance between Py and In
#     d = distance(merged_p, merged_i)
    
#     return 0  
  
def get_spike_trains(t, i, t_start=0):
    
    idx, frequency = np.unique(i, return_counts=True)
    t_end = max(t)
    spike_trains = []
    
    for j in range(len(idx)):
        subset = t[i==idx[j]]
        # spike_train = spk.SpikeTrain(subset, edges=t_end, is_sorted=True)
        spike_train = spk.SpikeTrain(subset, edges=(t_start, t_end), is_sorted=True)
        spike_trains.append(spike_train)    
        
    return spike_trains
            

# ISI distance profile
def isi_distance(spike_trains, population='', fig='', OS='local'):
    
    s = spk.isi_profile(spike_trains)
    
    if OS=='local':
        # Check if fig is an input
        if not fig:
            plt.figure()
        else:
            plt.figure(fig)
            
        x, y = s.get_plottable_data()
        plt.plot(x, y, label=population)
        plt.title("ISI distance")
        print("%s ISI distance: %.8f" %(population, s.avrg()))
        plt.show()
    
    return s.avrg()

# ISI distance matrix
def isi_distance_matrix(spike_trains, population='', fig='', OS='local'):
    
    s = spk.isi_distance_matrix(spike_trains)
    
    if OS=='local':
        # Check if fig is an input
        if not fig:
            plt.figure()
        else:
            plt.figure(fig)

        plt.imshow(s, interpolation='none')
        plt.title("ISI-distance")
    
# Spike distance profile 
def spike_distance(spike_trains, population='', fig='', OS='local'):
    # Check if fig is an input
    if not fig:
        plt.figure()
    else:
        plt.figure(fig)
    
    s = spk.spike_profile(spike_trains)
    if OS=='local':
        x,y = s.get_plottable_data()
        plt.plot(x,y, label=population)
        plt.title("SPIKE-distance")
        print("%s SPIKE-distance: %.8f" %(population, s.avrg()))    
        plt.show()
    
    return s.avrg()
    
# Spike syncrhonization
def spike_si(spike_trains, population='', OS='local'):
    # Spike Sync: 
    # Computes the spike synchronization value SYNC of the given spike trains.
    # The spike synchronization value is the computed as the total number of 
    # coincidences divided by the total number of spikes. The multivariate
    # SPIKE-Sync is again defined as the overall ratio of all coincidence 
    # values divided by the total number of spikes.    
    s = spk.spike_sync(spike_trains)
    print("%s SPIKE-synchronization: %.8f" %(population, s))
    
    return s

# PSTH
def plot_psth(spike_trains, bins=0.01, fig="", OS='local'):
    # TODO, this is weird. How does the funciton know when there is a stimulus?
    psth = spk.psth(spike_trains, bins)
    x, y = psth.get_plottable_data()
    
    if OS=='local':
        # Check if fig is an input
        if not fig:
            plt.figure()
        else:
            plt.figure(fig)
        
        plt.plot(x,y)
        plt.show()
    
    return x, y, psth

def isi_cv(st):
    # CV(isi) = std_/mean_isi
    if len(st) > 2:
        isi = np.diff(st.spikes)
        cv = np.std(isi)/np.mean(isi)
    else:
        cv = 0
    
    return cv

def fano(st):
    # CV(isi) = std_/mean_isi
    if len(st) > 2:
        isi = np.diff(st.spikes)
        fano = (np.std(isi)**2)/np.mean(isi)
    else:
        fano = 0
    
    return fano

def mean_isi_cv(spike_trains, ignore_zeros=True, population=''):
    cv_vect = np.zeros(len(spike_trains))
    for i in range(0,len(spike_trains)-1):
        cv_vect[i] = isi_cv(spike_trains[i])
    
    if ignore_zeros:
        cv_vect = np.delete(cv_vect, np.where(cv_vect==0))
        
    print("%s CV: %.8f +- %.8f" %(population, np.mean(cv_vect), np.std(cv_vect)))
            
    return np.mean(cv_vect), np.std(cv_vect)


def mean_fano(spike_trains, ignore_zeros=True, population=''):
    fano_vect = np.zeros(len(spike_trains))
    
    for i in range(0,len(spike_trains)-1):
        fano_vect[i] = fano(spike_trains[i])
    
    if ignore_zeros:
        fano_vect = np.delete(fano_vect, np.where(fano_vect==0))
        
    print("%s Fanofactor: %.8f" %(population, np.mean(fano_vect)))
            
    return np.mean(fano_vect)
    
    






    
    
    