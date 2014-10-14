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
#    Local_x, local_y, global_x, global_y, spacing, headway in feet
#    Veh_v in feet/second
#    Veh_acc in fee/second^2
#    Veh class Vehicle type: 1 - motorcycle, 2 - auto, 3 - truck
#

n_grid = 200
buckets = pickle.load(open('data/buckets_%d_%d.pi' % (n_grid, n_grid), 'rb'))

min_t = buckets['t_start'].min()
max_t = buckets['t_end'].max()

min_x = buckets['x_start'].min()
max_x = buckets['x_end'].max()

target_dir = '/Users/cusgadmin/projectARZ/US-101/vehicle-trajectory-data/'

foot_to_meter = 0.3048

all_data = []

id_offset = 0

for period in periods:
    print 'Loading period %s' % period
    file_path = '%s%s/trajectories-%s.csv' % (target_dir, period, period)
    data = pd.read_csv(file_path)
    data['vehicule_ID'] += id_offset
    id_offset = data['vehicule_ID'].max()
    all_data.append(data[data['veh_class'] == 2])
     
all_data = pd.concat(all_data)

all_data[['local_x', 'local_y', 'global_x', 
          'global_y', 'veh_length', 'veh_width']] *= foot_to_meter

all_data['t'] = (all_data['time_since_epoch_ms'] - all_data['time_since_epoch_ms'].min()) / 1000.0 
all_data['x'] = (all_data['global_x'] - all_data['global_x'].min())
all_data['id'] = all_data['vehicule_ID']
all_data = all_data[['x', 't', 'veh_v', 'id']]
all_data['count'] = 1
all_data = all_data[(all_data['x'] >= min_x)
                    &
                    (all_data['x'] <= max_x)
                    &
                    (all_data['t'] >= min_t)
                    &
                    (all_data['t'] <= max_t)]

ids = list(set(all_data['id']))

plt.scatter(buckets['rho'].values, buckets['q'].values, 
            marker = 'o', alpha = 0.2, lw = 0.0, c = 'grey')
for id in ids[::200][25:30]:
    traj = all_data[all_data['id'] == id]
    traj_x = traj['x'].values
    traj_t = traj['t'].values
    q_along = [buckets[(buckets['t_start'] <= traj_t[i])
                       &
                       (buckets['t_end'] >= traj_t[i])
                       &
                       (buckets['x_start'] <= traj_x[i])
                       &
                       (buckets['x_end'] >= traj_x[i])]['q'].mean()
               for i in range(len(traj))]
    rho_along = [buckets[(buckets['t_start'] <= traj_t[i])
                       &
                       (buckets['t_end'] >= traj_t[i])
                       &
                       (buckets['x_start'] <= traj_x[i])
                       &
                       (buckets['x_end'] >= traj_x[i])]['rho'].mean()
               for i in range(len(traj))]
    plt.scatter(rho_along[0], q_along[0], marker = 'x', c = 'k')
    plt.scatter(rho_along[-1], q_along[-1], marker = 'o', c = 'k')
    plt.plot(rho_along, q_along)
plt.savefig('plots/test4.png', dpi = 600)
plt.close()