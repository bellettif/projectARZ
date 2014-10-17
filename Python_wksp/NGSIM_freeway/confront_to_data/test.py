'''
Created on Oct 10, 2014

@author: cusgadmin
'''

import numpy as np
from matplotlib import pyplot as plt

f_1 = 2.3
f_2 = 3.7

T = 10.0 / f_1

n_points = 500

dt = T / float(n_points)

t_values = np.arange(n_points) * dt

input_signal = np.cos(2 * np.pi * f_1 * t_values + 0.5) + \
                np.sin(2 * np.pi * f_2 * t_values - 0.2)

transform = np.fft.rfft(input_signal)
freqs = np.fft.rfftfreq(len(input_signal), dt)

modules = np.absolute(transform)
angles = np.angle(transform)

reconstruct = np.zeros(n_points)
for i in range(len(freqs)):
    reconstruct += modules[i] * np.cos(2 * np.pi * freqs[i] * t_values + angles[i])
reconstruct /= len(freqs)
    
plt.plot(t_values, input_signal)
plt.plot(t_values, reconstruct)
plt.legend(('Input', 'Reconstruct'), 'upper right')
plt.show()
