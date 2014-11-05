'''
Created on Oct 8, 2014

@author: cusgadmin
'''

import numpy as np
from matplotlib import pyplot as plt
import cPickle as pickle
import os

from time_domain_responses import *
from fourier_transform import compute_input_fft
from fourier_transform import compute_inv_fft

PLOT_FOLDER = 'new_plots'
CALIBRATION_FOLDER = 'new_calibration_4'

if PLOT_FOLDER not in os.listdir('./'):
    os.mkdir(PLOT_FOLDER)
    
if CALIBRATION_FOLDER not in os.listdir('./'):
    os.mkdir(CALIBRATION_FOLDER)

N_TAUS = 80
TAU_VALUES = np.linspace(5, 80, N_TAUS)
#TAU_VALUES = [10, 15, 20, 25, 30]
#N_TAUS = len(TAU_VALUES)

PLOT_ALL = False
CALIBRATE_TAU = True

for n_grid in [80]:
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
    xi_1_data = (rho_star * lambda_2 / (lambda_1 - lambda_2)) * v_data + q_data
    xi_2_data = (rho_star * lambda_1 / (lambda_1 - lambda_2)) * v_data
    #
    #    Downstream boundary condition
    #
    xi_1_data_0 = xi_1_data[0]
    xi_2_data_0 = xi_2_data[0]
    #
    #    Upstream boundary condition
    #
    xi_1_data_L = xi_1_data[-1]
    xi_2_data_L = xi_2_data[-1]
    #
    mean_0, modules_0, args_0, freqs_0 = compute_input_fft(xi_1_data_0, dt)
    mean_L, modules_L, args_L, freqs_L = compute_input_fft(xi_2_data_L, dt)
    #
    reconstruct_xi_1_data = compute_inv_fft(mean_0,
                                         modules_0,
                                         args_0,
                                         freqs_0, 
                                         dt,
                                         len(xi_1_data_0)) 
    reconstruct_xi_2_data = compute_inv_fft(mean_L,
                                         modules_L,
                                         args_L,
                                         freqs_L, 
                                         dt,
                                         len(xi_2_data_L))
    if PLOT_ALL:
        plt.subplot(211)
        plt.plot(xi_1_data_0)
        plt.plot(reconstruct_xi_1_data)
        plt.title('xi_1(0,t) condition n = %d' % n_grid)
        plt.xlabel('t')
        plt.ylabel('xi_1')
        plt.legend(('xi_1', 'FFT inv'), 'upper right')
        plt.subplot(212)
        plt.plot(xi_2_data_L)
        plt.plot(reconstruct_xi_2_data)
        plt.title('xi_2(L,t) condition n = %d' % n_grid)
        plt.xlabel('t')
        plt.ylabel('xi_2')
        plt.legend(('xi_2', 'FFT inv'), 'upper right')
        plt.savefig('%s/Boundary_conditions_FFT_n_%d.png' % (PLOT_FOLDER, n_grid))
        plt.close()
    #
    q_mean_abs_errors = np.zeros(N_TAUS)
    v_mean_abs_errors = np.zeros(N_TAUS)
    xi_1_mean_abs_errors = np.zeros(N_TAUS)
    xi_2_mean_abs_errors = np.zeros(N_TAUS)
    for tau_index, TAU in enumerate(TAU_VALUES):
        print '\ttau = %.2f' % TAU
        #
        xi_1_sim = np.zeros((n_grid_x, n_grid_t))
        xi_2_sim = np.zeros((n_grid_x, n_grid_t))
        #
        for i in range(len(modules_0))[1:]:
            xi_1_sim += modules_0[i] * \
                        np.asanyarray([[fund_cos_1_1(2 * np.pi * freqs_0[i], t, 
                                                     args_0[i],# / (2 * np.pi * freqs_0[i]),
                                                     x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
            xi_2_sim += modules_0[i] * \
                        np.asanyarray([[fund_cos_1_2(2 * np.pi * freqs_0[i], t, 
                                                     args_0[i],# / (2 * np.pi * freqs_0[i]),
                                                     x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
        
        for i in range(len(modules_L))[1:]:
            xi_1_sim += modules_L[i] * \
                        np.asanyarray([[fund_cos_2_1(2 * np.pi * freqs_L[i], t,
                                                     args_L[i],# / (2 * np.pi * freqs_L[i]),
                                                     x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
            xi_2_sim += modules_L[i] * \
                        np.asanyarray([[fund_cos_2_2(2 * np.pi * freqs_L[i], t,
                                                     args_L[i],# / (2 * np.pi * freqs_L[i]),
                                                     x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
        #
        xi_1_sim /= float(len(freqs_0))
        xi_2_sim /= float(len(freqs_L))
        #
        xi_1_sim += mean_0 * np.asanyarray([[fund_H_1_1(t, x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
        xi_1_sim += mean_L * np.asanyarray([[fund_H_2_1(t, x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
        
        xi_2_sim += mean_0 * np.asanyarray([[fund_H_1_2(t, x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])[-1]
        
        xi_2_sim += mean_L * np.asanyarray([[fund_H_2_2(t, x, lambda_1, lambda_2, TAU, L)
                                          for t in t_values]
                                         for x in x_values])
        #
        q_sim = xi_1_sim - lambda_2 / lambda_1 * xi_2_sim
        v_sim = (lambda_1 - lambda_2) / (rho_star * lambda_1) * xi_2_sim
        if PLOT_ALL:
            #
            #    Check that boundary conditions do match
            #
            plt.subplot(211)
            plt.plot(xi_1_data_0)
            plt.plot(xi_1_sim[0])
            plt.title('xi_1(0,t) condition n = %d' % n_grid)
            plt.xlabel('t')
            plt.ylabel('xi_1')
            plt.legend(('xi_1_data', 'xi_1_sim'), 'upper right')
            plt.subplot(212)
            plt.plot(xi_2_data_L)
            plt.plot(xi_2_sim[-1])
            plt.title('xi_2(L,t) condition n = %d' % n_grid)
            plt.xlabel('t')
            plt.ylabel('xi_2')
            plt.legend(('xi_2', 'xi_2_sim'), 'upper right')
            plt.savefig('%s/Boundary_conditions_check_n=%d_tau=%.2f.png' % (PLOT_FOLDER, n_grid, TAU))
            plt.close()
            #
            #    Check maps in xi_1, xi_2 domain
            #
            #    Xi_1
            #        Values from data
            min_value = min(np.min(xi_1_data),
                            np.min(xi_1_sim),
                            np.min(xi_1_data - xi_1_sim))
            max_value = max(np.max(xi_1_data),
                            np.max(xi_1_sim),
                            np.max(xi_1_data - xi_1_sim))
            #
            #
            height = 14; width = 10
            #
            #
            plt.subplot(231)
            plt.imshow(xi_1_data[::-1], 
                       #vmin = min_value,
                       #vmax = max_value,
                       interpolation = 'None')
            plt.title('xi_1 data (veh/m)')
            plt.xlabel('t (seconds)')
            plt.ylabel('x (meters)')
            plt.colorbar()
            #        Simulated values
            plt.subplot(232)
            plt.imshow(xi_1_sim[::-1],
                       #vmin = min_value,
                       #vmax = max_value,
                       interpolation = 'None')
            plt.title('xi_1 sim (veh/m)')
            plt.xlabel('t (seconds)')
            plt.ylabel('x (meters)')
            plt.colorbar()
            #        Error
            plt.subplot(233)
            plt.imshow(xi_1_data - xi_1_sim[::-1],
                       #vmin = min_value,
                       #vmax = max_value,
                       interpolation = 'None')
            plt.title('xi_1 data - xi_1 sim (veh/m)')
            plt.xlabel('t (seconds)')
            plt.ylabel('x (meters)')
            plt.colorbar()
            #    Xi_2
            #        Values from data
            min_value = min(np.min(xi_2_data),
                            np.min(xi_2_sim),
                            np.min(xi_2_data - xi_2_sim))
            max_value = max(np.max(xi_2_data),
                            np.max(xi_2_sim),
                            np.max(xi_2_data - xi_2_sim))
            plt.subplot(234)
            plt.imshow(xi_2_data[::-1], 
                       #vmin = min_value,
                       #vmax = max_value,
                       interpolation = 'None')
            plt.title('xi_2 data (veh/m)')
            plt.xlabel('t (seconds)')
            plt.ylabel('x (meters)')
            plt.colorbar()
            #        Simulated values
            plt.subplot(235)
            plt.imshow(xi_2_sim[::-1], 
                       #vmin = min_value,
                       #vmax = max_value,
                       interpolation = 'None')
            plt.title('xi_2 sim (veh/m)')
            plt.xlabel('t (seconds)')
            plt.ylabel('x (meters)')
            plt.colorbar()
            #        Error
            plt.subplot(236)
            plt.imshow(xi_2_data[::-1] - xi_2_sim[::-1],
                       #vmin = min_value,
                       #vmax = max_value,
                       interpolation = 'None')
            plt.title('xi_2 data - xi_2 sim (veh/m)')
            plt.xlabel('t (seconds)')
            plt.ylabel('x (meters)')
            # Done
            plt.colorbar()
            fig = plt.gcf()
            fig.set_size_inches((height, width))
            plt.savefig('%s/xi_map_n=%d_tau=%.2f.png' % (PLOT_FOLDER, n_grid, TAU))
            plt.close()
            #
            #    Plot histogram of error in xi_1, xi_2 domain
            #
            plt.hist(np.ravel(xi_1_data[::-1] - xi_1_sim[::-1]), bins = 100)
            plt.xlabel('xi_1 error (veh/m)')
            plt.title('Histogram of xi_1 error')
            plt.ylabel('xi_1 data - xi_1 sim')
            plt.savefig('%s/xi_1_error_%d_%.2f.png' % (PLOT_FOLDER, n_grid, TAU))
            plt.close()
            #
            plt.hist(np.ravel(xi_2_data[::-1] - xi_2_sim[::-1]), bins = 100)
            plt.xlabel('xi_2 error (veh/m)')
            plt.title('Histogram of xi_2 error')
            plt.ylabel('xi_2 data - xi_2 sim')
            plt.savefig('%s/xi_2_error_%d_%.2f.png' % (PLOT_FOLDER, n_grid, TAU))
            plt.close()
            #
            #    Check maps in v, q domain
            #
            #    v
            #        Values from data
            #        Determining color scale
            min_value = min(np.min(v_data),
                            np.min(v_sim),
                            np.min(v_data - v_sim))
            max_value = max(np.max(v_data),
                            np.max(v_sim),
                            np.max(v_data - v_sim))
            #
            plt.subplot(231)
            plt.imshow(v_data[::-1], 
                       #vmin = min_value,
                       #vmax = max_value,
                       interpolation = 'None')
            plt.title('v data (m/s)')
            plt.xlabel('t (seconds)')
            plt.ylabel('x (meters)')
            plt.colorbar()
            #        Simulated values
            plt.subplot(232)
            plt.imshow(v_sim[::-1],
                       #vmin = min_value,
                       #vmax = max_value, 
                       interpolation = 'None')
            plt.title('v sim (m/s)')
            plt.xlabel('t (seconds)')
            plt.ylabel('x (meters)')
            plt.colorbar()
            #        Error
            plt.subplot(233)
            plt.imshow(v_data[::-1] - v_sim[::-1],
                       #vmin = min_value,
                       #vmax = max_value,
                       interpolation = 'None')
            plt.title('v data - v sim (m/s)')
            plt.xlabel('t (seconds)')
            plt.ylabel('x (meters)')
            plt.colorbar()
            #    q
            #        Values from data
            min_value = min(np.min(q_data),
                            np.min(q_sim),
                            np.min(q_data - q_sim))
            max_value = max(np.max(q_data),
                            np.max(q_sim),
                            np.max(q_data - q_sim))
            plt.subplot(234)
            plt.imshow(q_data[::-1], interpolation = 'None')
            plt.title('q data (veh/s)')
            plt.xlabel('t (seconds)')
            plt.ylabel('x (meters)')
            #        Simulated values
            plt.subplot(235)
            plt.imshow(q_sim[::-1], interpolation = 'None')
            plt.title('q sim (veh/s)')
            plt.xlabel('t (seconds)')
            plt.ylabel('x (meters)')
            #        Error
            plt.subplot(236)
            plt.imshow(q_data[::-1] - q_sim[::-1], interpolation = 'None')
            plt.title('q data - q sim (veh/s)')
            plt.xlabel('t (seconds)')
            plt.ylabel('x (meters)')
            plt.colorbar()
            # Done
            fig = plt.gcf()
            fig.set_size_inches((height,width))
            plt.savefig('%s/vq_map_n=%d_tau=%.2f.png' % (PLOT_FOLDER, n_grid, TAU))
            plt.close()
            #
            #    Plot histogram of error in v, q domain
            #
            plt.hist(np.ravel(v_data[::-1] - v_sim[::-1]), bins = 100)
            plt.xlabel('v error (m/s)')
            plt.title('Histogram of v error')
            plt.ylabel('v data - v sim (m/s)')
            plt.savefig('%s/v_error_%d_%.2f.png' % (PLOT_FOLDER, n_grid, TAU))
            plt.close()
            #
            plt.hist(np.ravel(q_data[::-1] - q_sim[::-1]), bins = 100)
            plt.xlabel('q error (veh/m)')
            plt.title('Histogram of q error')
            plt.ylabel('q data - q sim (veh/m)')
            plt.savefig('%s/q_error_%d_%.2f.png' % (PLOT_FOLDER, n_grid, TAU))
            plt.close()
        #
        #    Recording errors
        #
        xi_1_mean_abs_errors[tau_index] = np.mean(np.abs(xi_1_data - xi_1_sim))
        xi_2_mean_abs_errors[tau_index] = np.mean(np.abs(xi_2_data - xi_2_sim))
        q_mean_abs_errors[tau_index] = np.mean(np.abs(q_data - q_sim))
        v_mean_abs_errors[tau_index] = np.mean(np.abs(v_data - v_sim))
    if CALIBRATE_TAU:
        #
        #    Error for xi_1 and xi_2
        #
        plt.subplot(211)
        plt.title('Mean Absolute Error for xi_1 and xi_2')
        plt.plot(TAU_VALUES, xi_1_mean_abs_errors)
        plt.ylabel('MAE on xi_1 (veh/s)')
        plt.subplot(212)
        plt.plot(TAU_VALUES, xi_2_mean_abs_errors)
        plt.ylabel('MAE on xi_2 (veh/s)')
        plt.xlabel('Tau')
        #
        plt.savefig('%s/xi_1_xi_2_error_%d' % (CALIBRATION_FOLDER, n_grid))
        plt.close()
        #
        #    Error for v and q
        #
        plt.subplot(211)
        plt.title('Mean Absolute Error for q and v')
        plt.plot(TAU_VALUES, q_mean_abs_errors)
        plt.ylabel('MAE on q (veh/s)')
        plt.subplot(212)
        plt.plot(TAU_VALUES, v_mean_abs_errors)
        plt.ylabel('MAE on v (m/s)')
        plt.xlabel('Tau')
        #
        plt.savefig('%s/q_v_error_%d' % (CALIBRATION_FOLDER, n_grid))
        plt.close()
    

    