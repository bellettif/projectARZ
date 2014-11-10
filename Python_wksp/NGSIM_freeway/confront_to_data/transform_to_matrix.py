'''
Created on Oct 8, 2014

@author: cusgadmin
'''

import csv
import pandas as pd
import cPickle as pickle
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx


for n_grid in [80]:
    buckets = pd.read_pickle('../binned_data/buckets_%d_%d.pi' % (n_grid, n_grid))
    buckets = buckets[['x_start', 'x_end', 't_start', 't_end', 'rho', 'v', 'q']]
    buckets = buckets.sort(['x_start', 't_start'])
    x_values = list(set(buckets['x_start'].values))
    t_values = list(set(buckets['t_start'].values))
    x_values.sort()
    t_values.sort()
    mat_dict = {}
    for opt in ['rho', 'v', 'q']:
        values = [[buckets[(buckets['x_start'] <= x_values[i])
                               &
                               (buckets['x_end'] >= x_values[i])
                               &
                               (buckets['t_start'] <= t_values[j])
                               &
                               (buckets['t_end'] >= t_values[j])
                               ][opt].mean() for j in range(len(t_values))]
                      for i in range(len(x_values))]
        mat_dict[opt] = np.asanyarray(values)
    mat_dict['dx'] = buckets['x_end'].values[0] - buckets['x_start'].values[0]
    mat_dict['dt'] = buckets['t_end'].values[0] - buckets['t_start'].values[0]    
    pickle.dump(mat_dict, open('../matrices/mat_dict_%d_%d.pi' % (n_grid, n_grid), 'wb'))
    
