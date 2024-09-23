# -*- coding: utf-8 -*-
"""
Created on Sep 6 15:05:01 2023

@author: Artemio Soto-Breceda [artemios]

Loops through all the saved simulations from _0.mat to _100.mat of each combination of parameters and gets the average values.
"""

from os import walk
import scipy.io
from scipy.fft import fft, fftfreq
import numpy as np
from brian2 import *
prefs.codegen.target = 'numpy'  # use the Python fallback instead of C compilation
devices.device.shape = []       # This and the following line remove an annoying warning when brian2 is imported and loaded into RAM
devices.device.size = []
from scipy import signal
import matplotlib.pyplot as plt

# Read the folder
MODEL = 'cobn_stats'
FIRST = 'e'
SECOND = 'i'

folder = '/home/unimelb.edu.au/artemios/simulations/' + FIRST + '_vs_' + SECOND + '_' + MODEL + '/'
#save_folder = '/home/unimelb.edu.au/artemios/simulations/averages_' + MODEL + '/'
save_folder = '/home/unimelb.edu.au/artemios/Dropbox/University of Melbourne/LinuxRemoteDesktop/averages_' + MODEL + '/'
 
f = []

for (root, dirs, files) in walk(folder):
    f.extend(files)

# Set constants
g = 25 * nsiemens # Membrane conductance of pyramidal cells

# --- get values ---
lfp = {}
R_py_avg = {}
R_in_avg = {}
cv_avg = {}
si_avg = {}
spectra_avg = {}
balance_avg = {}
values = np.arange(0.5, 2.1, 0.1)
for e in values:
    for i in values:
        lst = list(filter(lambda k: '_%s%.2f_%s%.2f'%(FIRST,e,SECOND,i) in k, f))
        y = []
        spectra = []
        cv_in = []
        si_in = []
        R_py = []
        R_in = []
        balance = []
        for current_file in lst:
            data = scipy.io.loadmat(folder + current_file)

            # Get values
            i_pe =  data.get('i_pe')
            i_pi =  data.get('i_pi')
            SAMPLE_RATE =  1/data.get('lfp_dt')
            N = size(i_pe)
            
            try:
                cv_in.extend(data.get('cv_in'))
                si_in.extend(data.get('si_in'))

                R_py.extend(data.get('R_py'))
                R_in.extend(data.get('R_in'))
                
                local_field_potential = amp*(-(i_pe - i_pi)) / g
                y.extend(local_field_potential)

                b = 1e9*(i_pi[0,1000:] + i_pe[0,1000:]).mean()
                balance.extend([b])
                
                # FFT
                local_field_potential = local_field_potential / volt
                local_field_potential = local_field_potential[0,1000:] - local_field_potential[0,1000:].mean()
                yf = fft(local_field_potential)
                xf = fftfreq(N-1000, 1 / SAMPLE_RATE)
                
                #plt.plot(xf, np.abs(yf))
                #plt.show()
                
                spectra.extend([np.abs(yf)**2])
            except Exception as exc:
                print(str(exc))
                #print("Exception in: " + current_file)

        y = np.array(y)
        cv_in = np.array(cv_in)
        si_in = np.array(si_in) # [si_in!=0]
        # si_in = si_in[si_in>0]
        balance = np.array(balance)

        try:
            string = '%s%.2f_%s%.2f'%(FIRST,e,SECOND,i)
            string = string.replace('.','')
            spectra_avg[string] = np.array(spectra).mean(axis=0)
            lfp[string] = y.mean(axis=0)
            R_py_avg[string] = np.array(R_py).mean(axis=0)
            R_in_avg[string] = np.array(R_in).mean(axis=0)
            cv_avg[string] = cv_in.mean()
            si_avg[string] = si_in.mean()
            balance_avg[string] = balance.mean()

        except Exception as exc:
            print('%s%.2f_%s%.2f'%(FIRST,e,SECOND,i))
            print(str(exc))

# Add the frequency domain vector (x axis in the spectra)
spectra_avg['xf'] = xf

print('Analysis done, saving data...')

#scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'fft_x_axis.mat', mdict = xf)
scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'fft_avg.mat', mdict = spectra_avg)
scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'lfp_avg.mat', mdict = lfp)
scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'R_py_avg.mat', mdict = R_py_avg)
scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'R_in_avg.mat', mdict = R_in_avg)
scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'cv_avg.mat', mdict = cv_avg)
scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'si_avg.mat', mdict = si_avg)
scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'balance_avg.mat', mdict = balance_avg)

print('Done!')
