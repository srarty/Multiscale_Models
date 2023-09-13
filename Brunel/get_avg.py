# -*- coding: utf-8 -*-
"""
Created on Sep 6 15:05:01 2023

@author: Artemio Soto-Breceda [artemios]

Loops through all the saved simulations from _0.mat to _100.mat of each combination of parameters and gets the average values.
"""

from os import walk
import scipy.io
import numpy as np
from brian2 import *
prefs.codegen.target = 'numpy'  # use the Python fallback instead of C compilation
devices.device.shape = []       # This and the following line remove an annoying warning when brian2 is imported and loaded into RAM
devices.device.size = []
from scipy import signal
import matplotlib.pyplot as plt

# Read the folder
MODEL = 'cubn'
FIRST = 'ri'
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
values = np.arange(0.5, 2.1, 0.1)
for e in values:
    for i in values:
        lst = list(filter(lambda k: '_%s%.2f_%s%.2f'%(FIRST,e,SECOND,i) in k, f))
        y = []
        cv_in = []
        si_in = []
        R_py = []
        R_in = []
        for current_file in lst:
            data = scipy.io.loadmat(folder + current_file)

            # Get values
            i_pe =  data.get('i_pe')
            i_pi =  data.get('i_pi')
            
            try:
                cv_in.extend(data.get('cv_in'))
                si_in.extend(data.get('si_in'))

                R_py.extend(data.get('R_py'))
                R_in.extend(data.get('R_in'))

                y.extend(amp*(-(i_pe - i_pi)) / g)
            except:
                print(current_file)

        y = np.array(y)
        cv_in = np.array(cv_in)
        si_in = np.array(si_in) # [si_in!=0]
        try:
            string = '%s%.2f_%s%.2f'%(FIRST,e,SECOND,i)
            string = string.replace('.','')
            lfp[string] = y.mean(axis=0)
            R_py_avg[string] = np.array(R_py).mean(axis=0)
            R_in_avg[string] = np.array(R_in).mean(axis=0)
            cv_avg[string] = cv_in.mean()
            si_avg[string] = si_in.mean()
        except Exception as exc:
            print('%s%.2f_%s%.2f'%(FIRST,e,SECOND,i))
            print(str(exc))

scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'lfp_avg.mat', mdict = lfp)
scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'R_py_avg.mat', mdict = R_py_avg)
scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'R_in_avg.mat', mdict = R_in_avg)
scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'cv_avg.mat', mdict = cv_avg)
scipy.io.savemat(save_folder + FIRST +'vs' + SECOND + 'si_avg.mat', mdict = si_avg)
