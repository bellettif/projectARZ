'''
Created on Oct 8, 2014

@author: cusgadmin
'''

import numpy as np
from matplotlib import pyplot as plt
import cPickle as pickle

from time_domain_responses import *
from fourier_transform import compute_input_fft
from fourier_transform import compute_inv_fft

n_grid = 80

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
n_grid_x = q_matrix.shape[0]
n_grid_t = q_matrix.shape[1]

tau = 200.0

lambda_1 = np.median(v_matrix)

print lambda_1
print lambda_2

print n_grid_x
print n_grid_t

print dt

x_values = np.arange(n_grid_x) * dx
t_values = np.arange(n_grid_t) * dt

T = t_values[-1]
L = x_values[-1]

xsi_1_matrix = rho_star * lambda_2 / (lambda_1 - lambda_2) * v_matrix + q_matrix
xsi_2_matrix = rho_star * lambda_1 / (lambda_1 - lambda_2) * v_matrix

xsi_1_0 = xsi_1_matrix[0]
xsi_2_0 = xsi_2_matrix[0]

xsi_1_L = xsi_1_matrix[-1]
xsi_2_L = xsi_2_matrix[-1]

mean_0, modules_0, args_0, freqs_0 = compute_input_fft(xsi_1_0, dt)
mean_L, modules_L, args_L, freqs_L = compute_input_fft(xsi_2_L, dt)

reconstruct_xsi_1 = compute_inv_fft(mean_0,
                                     modules_0,
                                     args_0,
                                     freqs_0, 
                                     dt,
                                     len(xsi_1_0)) 
reconstruct_xsi_2 = compute_inv_fft(mean_L,
                                     modules_L,
                                     args_L,
                                     freqs_L, 
                                     dt,
                                     len(xsi_2_L))
plt.subplot(211)
plt.plot(xsi_1_0)
plt.subplot(212)
plt.plot(xsi_2_L)
plt.subplot(211)
plt.plot(reconstruct_xsi_1)
plt.subplot(212)
plt.plot(reconstruct_xsi_2)
plt.savefig('Boundary_cond.png')
plt.close()

xsi_1 = np.zeros((n_grid_x, n_grid_t))
xsi_2 = np.zeros((n_grid_x, n_grid_t))

for i in range(len(modules_0))[1:]:
    xsi_1 += modules_0[i] * \
                np.asanyarray([[fund_cos_1_1(2 * np.pi * freqs_0[i], t + args_0[i] / (2 * np.pi * freqs_0[i]),
                                           x, lambda_1, lambda_2, tau, L)
                                  for t in t_values]
                                 for x in x_values])
    xsi_2 += modules_0[i] * \
                np.asanyarray([[fund_cos_1_2(2 * np.pi * freqs_0[i], t + args_0[i] / (2 * np.pi * freqs_0[i]),
                                           x, lambda_1, lambda_2, tau, L)
                                  for t in t_values]
                                 for x in x_values])

for i in range(len(modules_L))[1:]:
    xsi_1 += modules_L[i] * \
                np.asanyarray([[fund_cos_2_1(2 * np.pi * freqs_L[i], t + args_L[i] / (2 * np.pi * freqs_L[i]),
                                           x, lambda_1, lambda_2, tau, L)
                                  for t in t_values]
                                 for x in x_values])
    xsi_2 += modules_L[i] * \
                np.asanyarray([[fund_cos_2_2(2 * np.pi * freqs_L[i], t + args_L[i] / (2 * np.pi * freqs_L[i]),
                                           x, lambda_1, lambda_2, tau, L)
                                  for t in t_values]
                                 for x in x_values])

xsi_1 /= float(len(freqs_0))
xsi_2 /= float(len(freqs_L))

xsi_1 += mean_0 * np.asanyarray([[fund_H_1_1(t, x, lambda_1, lambda_2, tau, L)
                                  for t in t_values]
                                 for x in x_values])
xsi_1 += mean_L * np.asanyarray([[fund_H_2_1(t, x, lambda_1, lambda_2, tau, L)
                                  for t in t_values]
                                 for x in x_values])

xsi_2 += mean_0 * np.asanyarray([[fund_H_1_2(t, x, lambda_1, lambda_2, tau, L)
                                  for t in t_values]
                                 for x in x_values])[-1]

xsi_2 += mean_L * np.asanyarray([[fund_H_2_2(t, x, lambda_1, lambda_2, tau, L)
                                  for t in t_values]
                                 for x in x_values])

q_pred = lambda_1 * xsi_2 - lambda_2 * xsi_1
v_pred = (lambda_1 - lambda_2) / (rho_star * lambda_1) * xsi_2

plt.plot(xsi_1_0)
plt.plot(xsi_1[0])
plt.show()

plt.plot(xsi_2_L)
plt.plot(xsi_2[-1])
plt.show()

plt.subplot(221)
plt.imshow(xsi_1_matrix, interpolation = 'None')
plt.subplot(222)
plt.imshow(xsi_1, interpolation = 'None')
plt.subplot(223)
plt.imshow(xsi_2_matrix, interpolation = 'None')
plt.subplot(224)
plt.imshow(xsi_2, interpolation = 'None')
plt.show()







