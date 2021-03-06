'''
Created on Sep 25, 2014

@author: francois
'''

import csv
import pandas as pd
import cPickle as pickle
import numpy as np
import os

from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

#
#    Set up global constants
#
FOOT_TO_METER = 0.3048
SAMPLING_RATE = 10.0 #s^-1
N_LANES = 5
PRODUCE_PLOTS = True
CONTROL_GRID = True

UNIT_DICT = {'v' : 'm/s',
             'q' : 'veh/s',
             'q_count' : 'veh/s',
             'median_ID' : 'NA',
             'rho' : 'veh/m'}
Q_STRING = {'q' : 'q',
            'q_count' : '$q_{count}$'}

#===============================================================================
# 
#    For US-101
#===============================================================================
INPUT_DIR = '/Users/cusgadmin/projectARZ/US-101/vehicle-trajectory-data/'
OUTPUT_DIR = '../binned_data_US101/'
#
TRIM_LOW = 578 * FOOT_TO_METER
TRIM_HIGH = (578 + 698) * FOOT_TO_METER
#
periods = ['0750am-0805am', '0805am-0820am', '0820am-0835am']
#
CONTROL_GRID_FOLDER = 'control_grid_US101'
PLOT_FOLDER = 'plots_US101'

#
#    For I-80
#
#===============================================================================
# INPUT_DIR = '/Users/cusgadmin/projectARZ/I-80/vehicle-trajectory-data/'
# OUTPUT_DIR = '../binned_data_I80/'
# #
# TRIM_LOW = 420 * FOOT_TO_METER
# TRIM_HIGH = np.Inf
# #
# periods = ['0500pm-0515pm', '0515pm-0530pm']
# #
# CONTROL_GRID_FOLDER = 'control_grid_I80'
# PLOT_FOLDER = 'plots_I80'
#===============================================================================

if CONTROL_GRID_FOLDER not in os.listdir('./'):
    os.mkdir(CONTROL_GRID_FOLDER)

if PLOT_FOLDER not in os.listdir('./'):
    os.mkdir(PLOT_FOLDER)


#
#    Information about the table
#
#===============================================================================
# columns = ['vehicule_ID', 'frame_ID', 'tot_frames', 'time_since_epoch_ms',  
#            'local_x', 'local_y', 'global_x', 'global_y', 'veh_length', 'veh_width', 'veh_class',
#            'veh_v', 'veh_acc', 'lane_id', 'prec_veh', 'follow_veh', 'spacing', 'headway']
#
#    Local_x, local_y, global_x, global_y, spacing, headway in feet
#    Veh_v in feet/second
#    Veh_acc in fee/second^2
#    Veh class Vehicle type: 1 - motorcycle, 2 - auto, 3 - truck
#
#===============================================================================


#
#    Reducing functions
#
def compute_metrics(x):
    result = {'rho': sum(x['count']) / (float(SAMPLING_RATE) * float(N_LANES) * dx * dt),
              'ids' : set(x['vehicule_ID']),
              'v': x['veh_v'].mean(),
              'n_traces' : len(x),
              'n_ids' : len(set(x['vehicule_ID']))}
    return pd.Series(result, name='metrics')
#
def q_counter(x):
    if (type(x['ids']) is set) and (type(x['shifted_ids']) is set):
        return len(x['ids'].intersection(x['shifted_ids'])) / (float(N_LANES) * dt)
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
#    Getting trajectory data
#
all_data = []
current_offset = 0
for period in periods:
    print 'Loading period %s' % period
    file_path = '%s%s/trajectories-%s.csv' % (INPUT_DIR, period, period)
    data = pd.read_csv(file_path)
    data['vehicule_ID'] += current_offset
    current_offset = data['vehicule_ID'].max()
    all_data.append(data[data['veh_class'] == 2])
all_data = pd.concat(all_data)
all_data[['local_x', 'local_y', 'global_x', 
          'global_y', 'veh_length', 'veh_width',
          'veh_v', 'veh_acc']] *= FOOT_TO_METER
all_data['t'] = (all_data['time_since_epoch_ms'] - all_data['time_since_epoch_ms'].min()) / 1000.0 
all_data['x'] = all_data['local_y']
all_data = all_data[['x', 't', 'veh_v', 'vehicule_ID']]
all_data = all_data[(all_data['x'] >= TRIM_LOW) &
                    (all_data['x'] <= TRIM_HIGH)]
all_data['count'] = 1

print len(all_data)

for n_grid in [80, 100, 120, 140]:
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
    #
    #    Computing estimators for rho and v
    #
    buckets = groups.apply(compute_metrics).reset_index()[['x',
                                                           't',
                                                           'rho', 
                                                           'v', 
                                                           'ids',
                                                           'n_traces',
                                                           'n_ids']]
    #
    #    Compute x and t
    #
    buckets['x_start'] = buckets['x'].apply(compute_start_point)
    buckets['x_end'] = buckets['x'].apply(compute_end_point)
    buckets['t_start'] = buckets['t'].apply(compute_start_point)
    buckets['t_end'] = buckets['t'].apply(compute_end_point)
    buckets['median_ID'] = buckets.apply(lambda x : np.median(list(x['ids'])), axis = 1)
    #
    #    Compute estimatores for q
    #
    buckets = buckets.sort(['t_start', 'x_start'])
    buckets['q'] = buckets['v'] * buckets['rho']
    buckets['shifted_ids'] = buckets['ids'].shift()
    buckets['q_count'] = buckets.apply(q_counter, axis = 1)
    #
    #   Deleting empty buckets
    #
    upper_cut_x = np.percentile(buckets['x_end'].values, 95)
    lower_cut_x = np.percentile(buckets['x_start'].values, 2.5)
    upper_cut_t = np.percentile(buckets['t_end'].values, 95)
    lower_cut_t = np.percentile(buckets['t_start'].values, 5)
    buckets = buckets[(buckets['x_start'] >= lower_cut_x)
                      &
                      (buckets['x_end'] <= upper_cut_x)
                      &
                      (buckets['t_start'] >= lower_cut_t)
                      &
                      (buckets['t_end'] <= upper_cut_t)]
    buckets = buckets[['t_start', 'x_start', 'x_end', 't_end', 'rho', 'v', 'ids', 'q', 'q_count', 'median_ID',
                       'n_traces', 'n_ids']]
    buckets = buckets.sort(['t_start', 'x_start'])
    if CONTROL_GRID:
        fontsize = 20
        #
        #    Plot histogram of number of traces in buckets
        #
        plt.hist(buckets['n_traces'].values, bins = 100)
        plt.title('Histogram of number of traces in buckets',
                  fontsize = fontsize)
        plt.ylabel('Population',
                   fontsize = fontsize)
        plt.xlabel('Number of traces',
                   fontsize = fontsize)
        plt.grid()
        plt.savefig('%s/traces_%d_%d' % (CONTROL_GRID_FOLDER, n_grid_t, n_grid_x))
        plt.close()
        print '%d, %d, 10 percentile of number of traces: %.2f' % (
                                               n_grid_t,
                                               n_grid_x,
                                               buckets['n_traces'].quantile(0.1))
        #
        #    Plot histogram of number of ids in buckets
        #
        plt.hist(buckets['n_ids'].values, bins = 100)
        plt.title('Histogram of number of ids in buckets',
                  fontsize = fontsize)
        plt.ylabel('Population',
                   fontsize = fontsize)
        plt.xlabel('Number of distinct ids',
                   fontsize = fontsize)
        plt.grid()
        plt.savefig('%s/ids_%d_%d' % (CONTROL_GRID_FOLDER, n_grid_t, n_grid_x))
        plt.close()
        print '%d, %d, 10 percentile of number of ids: %.2f' % (
                                               n_grid_t,
                                               n_grid_x,
                                               buckets['n_ids'].quantile(0.1))
    if PRODUCE_PLOTS:
        fontsize = 20
        #
        #    Plotting (t,x) map of v, q and rho
        #
        for opt in ['v', 'q', 'q_count', 'rho', 'median_ID']:
            sc = plt.scatter(buckets['t_start'].values,
                        buckets['x_start'].values,
                        c = buckets[opt],
                        lw = 0)
            plt.xlim((np.min(buckets['t_start'].values),
                      np.max(buckets['t_start'].values)))
            plt.ylim((np.min(buckets['x_start'].values),
                      np.max(buckets['x_start'].values)))
            plt.colorbar(sc)
            plt.title('Average %s (%s)' % (opt, UNIT_DICT[opt]), fontsize = fontsize)
            plt.xlabel('t (seconds)', fontsize = fontsize)
            plt.ylabel('x (meters)', fontsize = fontsize)
            plt.savefig('%s/%d_%d_%s_map.png' % (PLOT_FOLDER, n_grid_t, n_grid_x, opt))
            plt.close()
        #
        #    Sanity check q against q_count
        #   
        plt.scatter(buckets['q'].values, buckets['q_count'].values,
                    lw = 0, alpha = 0.2)
        plt.title(r'$q_{count}$ against q', fontsize = fontsize)
        plt.xlabel('q (veh/s)', fontsize = 16)
        plt.ylabel(r'$q_{count}$ (veh/s)', fontsize = fontsize)
        plt.xlim((0.0, 0.8))
        plt.ylim((0.0, 0.8))
        plt.plot(np.linspace(0.0, 0.8, 100), 
                 np.linspace(0.0, 0.8, 100),
                 c = 'r')
        plt.legend(('y=x', 'Scatter'), 'lower right', fontsize = fontsize)
        plt.savefig('%s/%d_%d_q_q_count.png' % 
                    (PLOT_FOLDER, n_grid_t, n_grid_x))
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
            plt.title(r'Fundamental diagram v %s' % q_opt, fontsize = fontsize)
            plt.xlabel('v (m/s)', fontsize = fontsize)
            plt.ylabel('%s (%s)' % (q_opt, UNIT_DICT[q_opt]), fontsize = fontsize)
            plt.savefig('%s/%d_%d_fundamental_diagram_v_%s.png' % 
                        (PLOT_FOLDER, n_grid_t, n_grid_x, q_opt))
            plt.close()
            #
            plt.scatter(buckets['rho'].values,
                        buckets[q_opt].values,
                        lw = 0, alpha = 0.2)
            plt.title(r'Fundamental diagram $\rho$ %s' % Q_STRING[q_opt], fontsize = fontsize)
            plt.xlabel(r'$\rho$ (veh/m)', fontsize = fontsize)
            plt.ylabel(r'%s (%s)' % (Q_STRING[q_opt], UNIT_DICT[q_opt]), fontsize = fontsize)
            plt.savefig('%s/%d_%d_fundamental_diagram_rho_%s.png' % 
                        (PLOT_FOLDER, n_grid_t, n_grid_x, q_opt))
            plt.close()
        #
        plt.scatter(buckets['rho'].values,
                    buckets['v'].values,
                    lw = 0, alpha = 0.2)
        plt.title(r'Fundamental diagram v $\rho$', fontsize = fontsize)
        plt.xlabel(r'$\rho$ (veh/m)', fontsize = fontsize)
        plt.ylabel('v (m/s)', fontsize = fontsize)
        plt.savefig('%s/%d_%d_fundamental_diagram_v_rho.png' % 
                    (PLOT_FOLDER, n_grid_t, n_grid_x))
        plt.close()
    #
    #    Pickling data
    #
    pickle.dump(buckets, open('%sbuckets_%d_%d.pi' % (OUTPUT_DIR, n_grid_t, n_grid_x), 'wb'))
    print 'Done'