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

for n_grid in [80, 100, 120, 140]:
    params = pickle.load(open('../system_params/%d_%d_params.pi' % (n_grid, n_grid), 'rb'))
    #
    lambda_1 = params['lambda_1']
    lambda_2 = params['lambda_2']
    rho_star = params['rho_star']
    v_star = params['v_star']
    q_star = params['q_star']
    #
    print 'N_grid = %d' % n_grid
    print 'Lambda_1 = %.2f' % lambda_1
    print 'Lambda_2 = %.2f' % lambda_2
    print 'rho_star = %.2f' % rho_star
    print 'v_star = %.2f' % v_star
    print 'q_star = %.2f' % q_star
    #
    mat_dict = pickle.load(open('../matrices/mat_dict_%d_%d.pi' % (n_grid, n_grid), 'rb'))
    #
    rho_data = mat_dict['rho'] - rho_star
    q_data = mat_dict['q'] - q_star
    v_data = mat_dict['v'] - v_star
    dx = mat_dict['dx']
    dt = mat_dict['dt']
    n_grid_x = q_data.shape[0]
    n_grid_t = q_data.shape[1]
    #
    x_values = np.arange(n_grid_x) * dx
    t_values = np.arange(n_grid_t) * dt
    #
    T = t_values[-1]
    L = x_values[-1]
    #
    xsi_1_data = rho_star * lambda_2 / (lambda_1 - lambda_2) * v_data + q_data
    xsi_2_data = rho_star * lambda_1 / (lambda_1 - lambda_2) * v_data
    #
    #    Downstream boundary condition
    #
    xsi_1_data_0 = xsi_1_data[0]
    xsi_2_data_0 = xsi_2_data[0]
    #
    #    Upstream boundary condition
    #
    xsi_1_data_L = xsi_1_data[-1]
    xsi_2_data_L = xsi_2_data[-1]
    #
    mean_0, modules_0, args_0, freqs_0 = compute_input_fft(xsi_1_data_0, dt)
    mean_L, modules_L, args_L, freqs_L = compute_input_fft(xsi_2_data_L, dt)
    #
    reconstruct_xsi_1_data = compute_inv_fft(mean_0,
                                         modules_0,
                                         args_0,
                                         freqs_0, 
                                         dt,
                                         len(xsi_1_data_0)) 
    reconstruct_xsi_2_data = compute_inv_fft(mean_L,
                                         modules_L,
                                         args_L,
                                         freqs_L, 
                                         dt,
                                         len(xsi_2_data_L))
    plt.subplot(211)
    plt.plot(xsi_1_data_0)
    plt.plot(reconstruct_xsi_1_data)
    plt.title('Xsi_1(0,t) condition n = %d' % n_grid)
    plt.xlabel('t')
    plt.ylabel('Xsi_1')
    plt.legend(('Xsi_1', 'FFT inv'), 'upper right')
    plt.subplot(212)
    plt.plot(xsi_2_data_L)
    plt.plot(reconstruct_xsi_2_data)
    plt.title('Xsi_2(L,t) condition n = %d' % n_grid)
    plt.xlabel('t')
    plt.ylabel('Xsi_2')
    plt.legend(('Xsi_2', 'FFT inv'), 'upper right')
    plt.savefig('plots/Boundary_conditions_FFT_n_%d.png' % n_grid)
    plt.close()
    #    
    for TAU in [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]:
        print '\ttau = %.2f' % TAU
        #
        xsi_1_sim = np.zeros((n_grid_x, n_grid_t))
        xsi_2_sim = np.zeros((n_grid_x, n_grid_t))
        #
        for i in range(len(modules_0))[1:]:
            xsi_1_sim += modules_0[i] * \
                        np.asanyarray([[fund_cos_1_1(2 * np.pi * freqs_0[i], t + args_0[i] / (2 * np.pi * freqs_0[i]),
                                                   x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
            xsi_2_sim += modules_0[i] * \
                        np.asanyarray([[fund_cos_1_2(2 * np.pi * freqs_0[i], t + args_0[i] / (2 * np.pi * freqs_0[i]),
                                                   x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
        
        for i in range(len(modules_L))[1:]:
            xsi_1_sim += modules_L[i] * \
                        np.asanyarray([[fund_cos_2_1(2 * np.pi * freqs_L[i], t + args_L[i] / (2 * np.pi * freqs_L[i]),
                                                   x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
            xsi_2_sim += modules_L[i] * \
                        np.asanyarray([[fund_cos_2_2(2 * np.pi * freqs_L[i], t + args_L[i] / (2 * np.pi * freqs_L[i]),
                                                   x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
        #
        xsi_1_sim /= float(len(freqs_0))
        xsi_2_sim /= float(len(freqs_L))
        #
        xsi_1_sim += mean_0 * np.asanyarray([[fund_H_1_1(t, x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
        xsi_1_sim += mean_L * np.asanyarray([[fund_H_2_1(t, x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
        
        xsi_2_sim += mean_0 * np.asanyarray([[fund_H_1_2(t, x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])[-1]
        
        xsi_2_sim += mean_L * np.asanyarray([[fund_H_2_2(t, x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
        #
        q_sim = lambda_1 * xsi_1_sim - lambda_2 * xsi_2_sim
        v_sim = (lambda_1 - lambda_2) / (rho_star * lambda_1) * xsi_2_sim
        #
        plt.subplot(211)
        plt.plot(xsi_1_data_0)
        plt.plot(xsi_1_sim[0])
        plt.title('Xsi_1(0,t) condition n = %d' % n_grid)
        plt.xlabel('t')
        plt.ylabel('Xsi_1')
        plt.legend(('Xsi_1_data', 'Xsi_1_sim'), 'upper right')
        plt.subplot(212)
        plt.plot(xsi_2_data_L)
        plt.plot(xsi_2_sim[-1])
        plt.title('Xsi_2(L,t) condition n = %d' % n_grid)
        plt.xlabel('t')
        plt.ylabel('Xsi_2')
        plt.legend(('Xsi_2', 'Xsi_2_sim'), 'upper right')
        plt.savefig('plots/Boundary_conditions_check_n=%d_tau=%.2f.png' % (n_grid, TAU))
        plt.close()
        #
        plt.subplot(221)
        plt.imshow(xsi_1_data[::-1], interpolation = 'None')
        plt.title('Xsi_1 data')
        plt.xlabel('t')
        plt.ylabel('x')
        plt.colorbar()
        plt.subplot(222)
        plt.imshow(xsi_1_sim[::-1], interpolation = 'None')
        plt.title('Xsi_1 sim')
        plt.xlabel('t')
        plt.ylabel('x')
        plt.colorbar()
        plt.subplot(223)
        plt.imshow(xsi_2_data[::-1], interpolation = 'None')
        plt.title('Xsi_2 data')
        plt.xlabel('t')
        plt.ylabel('x')
        plt.colorbar()
        plt.subplot(224)
        plt.imshow(xsi_2_sim[::-1], interpolation = 'None')
        plt.title('Xsi_2 sim')
        plt.xlabel('t')
        plt.ylabel('x')
        plt.colorbar()
        plt.savefig('plots/Xsi_map_n=%d_tau=%.2f.png' % (n_grid, TAU))
        plt.close()
        #
        plt.subplot(221)
        plt.imshow(v_data[::-1], interpolation = 'None')
        plt.title('v data')
        plt.xlabel('t')
        plt.ylabel('x')
        plt.colorbar()
        plt.subplot(222)
        plt.imshow(v_sim[::-1], interpolation = 'None')
        plt.title('v sim')
        plt.xlabel('t')
        plt.ylabel('x')
        plt.colorbar()
        plt.subplot(223)
        plt.imshow(q_data[::-1], interpolation = 'None')
        plt.title('q data')
        plt.xlabel('t')
        plt.ylabel('x')
        plt.colorbar()
        plt.subplot(224)
        plt.imshow(q_sim[::-1], interpolation = 'None')
        plt.title('q sim')
        plt.xlabel('t')
        plt.ylabel('x')
        plt.colorbar()
        plt.savefig('plots/vq_map_n=%d_tau=%.2f.png' % (n_grid, TAU))
        plt.close()
