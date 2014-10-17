'''
Created on Oct 10, 2014

@author: cusgadmin
'''

import numpy as np
from matplotlib import pyplot as plt
import cPickle as pickle
import pandas as pd

for n_grid in [80, 100, 120, 140]:
    buckets = pd.read_pickle('../binned_data/buckets_%d_%d.pi' % (n_grid, n_grid))
    q_values = buckets['q'].values
    rho_values = buckets['rho'].values
    v_values = buckets['v'].values
    #
    coeffs = np.linalg.lstsq(np.vstack([rho_values, np.ones(len(rho_values))]).T, q_values)
    x_values = np.linspace(np.min(rho_values),
                           np.max(rho_values),
                           num = 100)
    y_values = coeffs[0][0] * x_values + coeffs[0][1]
    plt.plot(x_values, y_values, c = 'red')
    plt.scatter(rho_values, q_values, alpha = 0.2, lw = 0)
    plt.title('Fundamental diagram fitting n = %d' % n_grid)
    plt.xlabel('rho')
    plt.ylabel('q')
    plt.savefig('plots/Fundamental diagram fitting n = %d' % n_grid)
    plt.close()
    v_star = np.mean(v_values)
    rho_star = np.mean(rho_values)
    q_star = v_star * rho_star
    lambda_1 = v_star
    lambda_2 = coeffs[0][0]
    print 'Lambda_1 = %.2f, lambda_2 = %.2f' % (lambda_1, lambda_2)
    parameters = {'v_star' : v_star,
                 'rho_star' : rho_star,
                 'q_star' : q_star,
                 'lambda_1' : lambda_1,
                 'lambda_2' : lambda_2}
    pickle.dump(parameters, open('../system_params/%d_%d_params.pi' % (n_grid, n_grid), 'wb'))