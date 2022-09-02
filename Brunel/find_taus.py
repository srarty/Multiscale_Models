# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 13:59:51 2022

@author: Artemio Soto-Breceda [artemios]
"""
from synaptic_functions_exploration import *

for alpha in arange(-61,0.1,1.2):
    monitor = synaptic_functions_exploration(alpha_ii = alpha, PLOT = False)
    
    # RECURSIVE Interneurons:
    save = {'ipsp': monitor.v1*1e3}
    scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/taus_fit/inhibitory_ipsp_%s.mat' %round(alpha,3), mdict=save)

for alpha in arange(-335, 0.1, 6.5):
    monitor = synaptic_functions_exploration(alpha_ie = alpha, PLOT = False)
        
    # Pyramidal:
    save = {'ipsp': monitor.v1*1e3}
    scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/taus_fit/pyramidal_ipsp_%s.mat' %round(alpha,3), mdict=save)

for alpha in arange(0, 691, 13.53):
    monitor = synaptic_functions_exploration(alpha_ei = alpha, PLOT = False)
        
    # Interneurons:
    save = {'epsp': monitor.v1*1e3}
    scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/taus_fit/inhibitory_epsp_%s.mat' %round(alpha,3), mdict=save)

for alpha in arange(0, 891, 35.6):
    monitor = synaptic_functions_exploration(alpha_ee = alpha, PLOT = False)
        
    # RECURSIVE Pyramidal:
    save = {'epsp': monitor.v1*1e3}
    scipy.io.savemat('C://Users/artemios/Documents/Multiscale_Models_Data/taus_fit/pyramidal_epsp_%s.mat' %round(alpha,3), mdict=save)

