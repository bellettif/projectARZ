'''
Created on Oct 8, 2014

@author: cusgadmin
'''

import numpy as np
from matplotlib import pyplot as plt
import cPickle as pickle

from time_domain_responses import *
from fourier_transform import compute_input_fft

w = 0.03

n_grid = 100

params = pickle.load(open('system_params/%d_%d_params.pi' % (n_grid, n_grid), 'rb'))

lambda_1 = params['lambda_1']
lambda_2 = params['lambda_2']
rho_star = params['rho_star']

mat_dict = pickle.load(open('matrices/mat_dict_%d_%d.pi' % (n_grid, n_grid), 'rb'))

rho_matrix = mat_dict['rho']
q_matrix = mat_dict['q']
v_matrix = mat_dict['v']
dx = mat_dict['dx']
dt = mat_dict['dt']

tau = 1.0
L = 100.0
T = 100.0

x_values = np.linspace(0, L, 100)
t_values = np.linspace(0, T, 100)

xsi_1 = [[fund_cos_1_1(2 * np.pi * w, t, x, lambda_1, lambda_2, tau, L)
          for t in t_values] for x in x_values]
xsi_2 = [[fund_cos_1_2(2 * np.pi * w, t, x, lambda_1, lambda_2, tau, L)
          for t in t_values] for x in x_values]

xsi_1 = np.asanyarray(xsi_1)
xsi_2 = np.asanyarray(xsi_2)

plt.plot(xsi_1[0])
plt.show()
plt.plot(xsi_2[-1])
plt.show()

q = xsi_1 - lambda_2 / lambda_1 * xsi_2
v = (lambda_1 - lambda_2) / (rho_star * lambda_1) * xsi_2

plt.subplot(221)
plt.title('xsi_1')
plt.imshow(xsi_1, interpolation = 'None')
plt.subplot(222)
plt.title('xsi_2')
plt.imshow(xsi_2, interpolation = 'None')
plt.subplot(223)
plt.title('q')
plt.imshow(q, interpolation = 'None')
plt.subplot(224)
plt.title('v')
plt.imshow(v, interpolation = 'None')
plt.show()

xsi_1 = [[fund_cos_2_1(2 * np.pi * w, t, x, lambda_1, lambda_2, tau, L)
          for t in t_values] for x in x_values]
xsi_2 = [[fund_cos_2_2(2 * np.pi * w, t, x, lambda_1, lambda_2, tau, L)
          for t in t_values] for x in x_values]

xsi_1 = np.asanyarray(xsi_1)
xsi_2 = np.asanyarray(xsi_2)

q = xsi_1 - lambda_2 / lambda_1 * xsi_2
v = (lambda_1 - lambda_2) / (rho_star * lambda_1) * xsi_2

plt.plot(xsi_1[0])
plt.show()
plt.plot(xsi_2[-1])
plt.show()

plt.subplot(221)
plt.title('xsi_1')
plt.imshow(xsi_1, interpolation = 'None')
plt.subplot(222)
plt.title('xsi_2')
plt.imshow(xsi_2, interpolation = 'None')
plt.subplot(223)
plt.title('q')
plt.imshow(q, interpolation = 'None')
plt.subplot(224)
plt.title('v')
plt.imshow(v, interpolation = 'None')
plt.show()