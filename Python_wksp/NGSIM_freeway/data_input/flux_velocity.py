'''
Created on Sep 25, 2014

@author: francois
'''

import csv
import pandas
import cPickle as pickle
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from data_input.n_datapoints import foot_to_meter

periods = ['0750am-0805am',
           '0805am-0820am',
           '0820am-0835am']

target_dir = '/Users/cusgadmin/projectARZ/US-101/vehicle-trajectory-data/'

foot_to_meter = 0.3048

grid_points = 100

for period in periods:
    print 'Doing period %s' % period
    file_path = '%s%s/trajectories-%s.csv' % (target_dir, period, period)
    data = pandas.read_csv(file_path)
    #
    x_values = data['global_x'].values * foot_to_meter
    t_values = data['time_since_epoch_ms'].values
    v_values = data['veh_v'].values * foot_to_meter
    #
    zipped_data = zip(x_values, t_values, v_values)
    #
    min_x = np.min(x_values)
    max_x = np.max(x_values)
    x_buckets = np.linspace(min_x, max_x, num = grid_points, endpoint = True)
    #
    min_t = np.min(t_values)
    max_t = np.max(t_values)
    t_buckets = np.linspace(min_t, max_t, num = grid_points, endpoint = True)
    #
    min_v = np.min(v_values)
    max_v = np.min(v_values)
    #
    v_buckets = np.zeros((grid_points - 1, grid_points - 1), dtype = np.float)
    q_buckets = np.zeros((grid_points - 1, grid_points - 1), dtype = np.float)
    rho_buckets = np.zeros((grid_points - 1, grid_points - 1), dtype = np.float)
    #
    for i in xrange(grid_points - 1):
        start_x = x_buckets[i]
        stop_x = x_buckets[i+1]
        delta_x = float(stop_x - start_x)
        for j in xrange(grid_points - 1):    
            start_t = t_buckets[j]
            stop_t = t_buckets[j+1]
            delta_t = float(stop_t - start_t)
            sub_selection = [x[-1] for x in 
                                      filter(lambda x : x[0] >= start_x 
                                                        and 
                                                        x[0] < stop_x
                                                        and 
                                                        x[1] >= start_t
                                                        and 
                                                        x[1] < stop_t, 
                                                        zipped_data)]
            n_veh = len(sub_selection)
            v_buckets[i,j] = np.mean(sub_selection) if n_veh > 0 else 0
            q_buckets[i,j] = v_buckets[i,j] * n_veh / (delta_x * delta_t)
            rho_buckets[i,j] = n_veh / delta_x
    #       
    pickle.dump(v_buckets, open('v_buckets_%s.pi' % period, 'wb'))        
    pickle.dump(q_buckets, open('q_buckets_%s.pi' % period, 'wb'))
    pickle.dump(rho_buckets, open('rho_buckets_%s.pi' % period, 'wb'))
    pickle.dump(x_buckets, open('x_buckets_%s.pi' % period, 'wb'))        
    pickle.dump(t_buckets, open('t_buckets_%s.pi' % period, 'wb'))
    #
    plt.imshow(v_buckets, interpolation='None')
    plt.savefig('v_buckets_%s.png' % period, dpi = 300)
    plt.close()
    #
    plt.imshow(q_buckets,interpolation = 'None')
    plt.savefig('q_buckets_%s.png' % period, dpi = 300)
    plt.close()
    print 'Done with period %s' % period
    
print 'All done'