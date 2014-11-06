'''
Created on Oct 10, 2014

@author: cusgadmin
'''

import os
import numpy as np
from matplotlib import pyplot as plt
import cPickle as pickle
import pandas as pd

PLOT_DIR = 'plots'

if PLOT_DIR not in os.listdir('./'):
    os.mkdir(PLOT_DIR)

for n_grid in [80]:
    buckets = pd.read_pickle('../binned_data/buckets_%d_%d.pi' % (n_grid, n_grid))
    q_values = buckets['q'].values
    rho_values = buckets['rho'].values
    v_values = buckets['v'].values
    #
    #    FITTING
    #
    coeffs = np.linalg.lstsq(np.vstack([rho_values, np.ones(len(rho_values))]).T, q_values)
    x_values = np.linspace(np.min(rho_values),
                           np.max(rho_values),
                           num = 100)
    y_values = coeffs[0][0] * x_values + coeffs[0][1]
    v_star = np.mean(v_values)
    q_star = np.mean(q_values)
    rho_star = q_star / v_star
    lambda_1 = v_star
    lambda_2 = coeffs[0][0]
    #
    #    COMPUTE R^2
    #
    q_pred = coeffs[0][0] * rho_values + coeffs[0][1]
    u = np.sum((q_values - q_pred) ** 2)
    v = np.sum((q_values - np.mean(q_values)) ** 2)
    r_2 = 1 - u/v
    #
    #    PLOT
    #
    plt.plot(x_values, y_values, c = 'red')
    plt.scatter(rho_values, q_values, alpha = 0.2, lw = 0)
    plt.scatter(rho_star, q_star, c = 'red', s = 50)
    plt.title('Fundamental diagram fitting n = %d' % n_grid)
    plt.xlabel('rho (veh/m)')
    plt.ylabel('q_count (veh/s)')
    plt.legend(('Fitted model', 'Data scatter', '(rho_star, q_star)'), 'upper right')
    plt.savefig('plots/Fundamental diagram fitting n = %d' % n_grid)
    plt.close()
    #
    #
    #
    print 'Lambda_1 = %.2f, lambda_2 = %.2f' % (lambda_1, lambda_2)
    print 'v_star = %.2f' % v_star
    print 'q_star = %.2f' % q_star
    print 'rho_star = %.2f' % rho_star
    print 'r^2 = %.2f' % r_2
    parameters = {'v_star' : v_star,
                 'rho_star' : rho_star,
                 'q_star' : q_star,
                 'lambda_1' : lambda_1,
                 'lambda_2' : lambda_2}
    pickle.dump(parameters, open('../system_params/%d_%d_params.pi' % (n_grid, n_grid), 'wb'))