'''
Created on Oct 8, 2014

@author: cusgadmin
'''

import numpy as np
from matplotlib import pyplot as plt

def compute_input_fft(input_sig, dt):
    mean_value = np.mean(input_sig)
    fft = np.fft.rfft(input_sig - np.mean(input_sig))
    modules = np.absolute(fft)
    args = np.angle(fft)
    freqs = np.fft.rfftfreq(len(input_sig), dt)
    return mean_value, modules, args, freqs

def compute_inv_fft(signal_mean, modules, args, freqs, dt, n_points):
    t_values = dt * np.arange(n_points)
    values = np.zeros(n_points)
    for i in range(len(freqs)):
        values += modules[i] * np.cos(2 * np.pi * freqs[i] * t_values + args[i])
    return values / float(len(freqs)) + np.ones(n_points) * signal_mean