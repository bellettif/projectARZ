'''
Created on Oct 10, 2014

@author: cusgadmin
'''

import numpy as np
from matplotlib import pyplot as plt

from fourier_transform import *

f_1 = 2.3
f_2 = 3.7

T = 10.0 / f_1

n_points = 500

dt = T / float(n_points)

t_values = np.arange(n_points) * dt

input_signal = np.cos(2 * np.pi * f_1 * t_values + 0.5) + \
                np.sin(2 * np.pi * f_2 * t_values - 0.2)

mean_value, modules, angles, freqs = compute_input_fft(input_signal, dt)
reconstruct = compute_inv_fft(mean_value, modules, angles, freqs, dt, n_points)
    
plt.plot(t_values, input_signal)
plt.plot(t_values, reconstruct)
plt.legend(('Input', 'Reconstruct'), 'upper right')
plt.show()
