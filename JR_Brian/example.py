# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 10:20:50 2022

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

