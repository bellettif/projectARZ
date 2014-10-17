'''
Created on Sep 25, 2014

@author: francois
'''

import csv
import pandas as pd
import cPickle as pickle
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

periods = ['0750am-0805am',
           '0805am-0820am',
           '0820am-0835am']

columns = ['vehicule_ID', 'frame_ID', 'tot_frames', 'time_since_epoch_ms',  
           'local_x', 'local_y', 'global_x', 'global_y', 'veh_length', 'veh_width', 'veh_class',
           'veh_v', 'veh_acc', 'lane_id', 'prec_veh', 'follow_veh', 'spacing', 'headway']

#
#    Reducing function
#
def compute_metrics(x):
    result = {'rho': sum(x['count']) * 1.0 / (dx * 10.0 * dt) * 0.2,
              'rho_count' : len(set(x['vehicule_ID'])) / (dx * dt) * 0.2,
              'ids' : set(x['vehicule_ID']),
              'v': x['veh_v'].mean()}
    return pd.Series(result, name='metrics')
#
def q_counter(x):
    if (type(x['ids']) is set) and (type(x['shifted_ids']) is set):
        return len(x['ids'].intersection(x['shifted_ids']))
    else:
        return 0

#
#    Utils for panda arrays
#
def compute_start_point(interval):
    begin, end = interval.split(',')
    begin = float(begin[1:])
    return begin
#
def compute_end_point(interval):
    begin, end = interval.split(',')
    end = float(end[1:-1])
    return end

#
#    Local_x, local_y, global_x, global_y, spacing, headway in feet
#    Veh_v in feet/second
#    Veh_acc in fee/second^2
#    Veh class Vehicle type: 1 - motorcycle, 2 - auto, 3 - truck
#

target_dir = '/Users/cusgadmin/projectARZ/US-101/vehicle-trajectory-data/'

foot_to_meter = 0.3048

#
#    Getting trajectory data
#
all_data = []
for period in periods:
    print 'Loading period %s' % period
    file_path = '%s%s/trajectories-%s.csv' % (target_dir, period, period)
    data = pd.read_csv(file_path)
    print data['vehicule_ID'].min()
    print data['vehicule_ID'].max()
    print '\n'
    all_data.append(data[data['veh_class'] == 2])
all_data = pd.concat(all_data)
all_data[['local_x', 'local_y', 'global_x', 
          'global_y', 'veh_length', 'veh_width',
          'veh_v', 'veh_acc']] *= foot_to_meter
all_data['t'] = (all_data['time_since_epoch_ms'] - all_data['time_since_epoch_ms'].min()) / 1000.0 
all_data['x'] = all_data['local_y']
all_data = all_data[['x', 't', 'veh_v', 'vehicule_ID']]
all_data['count'] = 1

print all_data

for n_grid in [80, 100]:
    n_grid_x = n_grid_t = n_grid
    print 'Doing grid with %d points' % n_grid
    #
    #    Grid slicing
    #
    bins_x, dx = np.linspace(all_data['x'].min(),
                             all_data['x'].max(),
                             n_grid_x, retstep = True)
    bins_t, dt = np.linspace(all_data['t'].min(),
                             all_data['t'].max(),
                             n_grid_t, retstep = True)
    groups = all_data.groupby([pd.cut(all_data['x'], bins_x),
                               pd.cut(all_data['t'], bins_t)])
    buckets = groups.apply(compute_metrics).reset_index()[['x',
                                                           't',
                                                           'rho', 
                                                           'v', 
                                                           'ids']]
    buckets['x_start'] = buckets['x'].apply(compute_start_point)
    buckets['x_end'] = buckets['x'].apply(compute_end_point)
    buckets['t_start'] = buckets['t'].apply(compute_start_point)
    buckets['t_end'] = buckets['t'].apply(compute_end_point)
    buckets['q'] = buckets['v'] * buckets['rho']
    buckets['shifted_ids'] = buckets['ids'].shift()
    buckets['q_count'] = buckets.apply(q_counter, axis = 1) * 0.2
    buckets['median_ID'] = buckets.apply(lambda x : x['vehicule_ID'].median())
    #
    #    Deleting empty buckets
    #
    upper_cut_x = np.percentile(buckets['x_end'].values, 95)
    lower_cut_x = np.percentile(buckets['x_start'].values, 5)
    upper_cut_t = np.percentile(buckets['t_end'].values, 95)
    lower_cut_t = np.percentile(buckets['t_start'].values, 5)
    buckets = buckets[(buckets['x_start'] >= lower_cut_x)
                      &
                      (buckets['x_end'] <= upper_cut_x)
                      &
                      (buckets['t_start'] >= lower_cut_t)
                      &
                      (buckets['t_end'] <= upper_cut_t)]
    buckets = buckets[['t_start', 'x_start', 'rho', 'v', 'ids', 'q', 'q_count']]
    buckets = buckets.sort(['t_start', 'x_start'])
    #
    #    Plotting (t,x) map of v, q and rho
    #
    for opt in ['v', 'q', 'q_count', 'rho', 'rho_count', 'median_ID']:
        sc = plt.scatter(buckets['t_start'].values,
                    buckets['x_start'].values,
                    c = buckets[opt],
                    lw = 0)
        plt.colorbar(sc)
        plt.title('Average %s' % opt)
        plt.xlabel('t (seconds)')
        plt.ylabel('x (meters)')
        plt.savefig('plots/%d_%d_%s_map.png' % (n_grid_t, n_grid_x, opt))
        plt.close()
    #
    #    Sanity check q against q_count
    #   
    plt.scatter(buckets['q'].values, buckets['q_count'].values,
                lw = 0, alpha = 0.2)
    plt.title('q against q_count')
    plt.xlabel('q')
    plt.ylabel('q_count')
    plt.savefig('plots/%d_%d_q_q_count.png' % (n_grid_t, n_grid_x))
    plt.close()
    #
    #    Sanity check rho against rho_count
    #
    for q_opt in ['q_count', 'q']:
        #
        #    Plotting fundamental diagram
        #
        plt.scatter(buckets['v'].values,
                    buckets[q_opt].values,
                    lw = 0, alpha = 0.2)
        plt.title('Fundamental diagram v %s' % q_opt)
        plt.xlabel('v')
        plt.ylabel(q_opt)
        plt.savefig('plots/%d_%d_fundamental_diagram_v_%s.png' % (n_grid_t, n_grid_x, q_opt))
        plt.close()
        #
        plt.scatter(buckets['rho'].values,
                    buckets[q_opt].values,
                    lw = 0, alpha = 0.2)
        plt.title('Fundamental diagram rho %s' % q_opt)
        plt.xlabel('rho')
        plt.ylabel(q_opt)
        plt.savefig('plots/%d_%d_fundamental_diagram_rho_%s.png' % (n_grid_t, n_grid_x, q_opt))
        plt.close()
    #
    #    Pickling data
    #
    pickle.dump(buckets, open('data/buckets_%d_%d.pi' % (n_grid_t, n_grid_x), 'wb'))
    print 'Done'