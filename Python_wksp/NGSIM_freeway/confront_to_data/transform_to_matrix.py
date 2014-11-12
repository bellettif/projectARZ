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
from fundamental_diagram_fitting.linear_fitting import OUTPUT_FOLDER
import os

INPUT_FOLDER = '../binned_data_I80'
OUTPUT_FOLDER = 'matrices_I80'

if OUTPUT_FOLDER not in os.listdir('../'):
    os.mkdir('../' + OUTPUT_FOLDER)

for n_grid in [80]:
    buckets = pd.read_pickle('%s/buckets_%d_%d.pi' % (INPUT_FOLDER, n_grid, n_grid))
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
    pickle.dump(mat_dict, open('../%s/mat_dict_%d_%d.pi' % (OUTPUT_FOLDER, n_grid, n_grid), 'wb'))
    
