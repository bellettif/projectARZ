'''
Created on Sep 23, 2014

@author: francois
'''

import csv
import pandas
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

periods = ['0750am-0805am',
           '0805am-0820am',
           '0820am-0835am']


target_dir = '/Users/cusgadmin/projectARZ/US-101/vehicle-trajectory-data/'

foot_to_meter = 0.3048

all_data = []

for period in periods:
    file_path = '%s%s/trajectories-%s.csv' % (target_dir, period, period)
    data = pandas.read_csv(file_path)
    all_data.append(data)
    #
    min_t = float(data['time_since_epoch_ms'].min())
    min_v = float(data['veh_v'].min())
    max_v = float(data['veh_v'].max())
    print 'v belongs to (%f, %f)' % (min_v, max_v)
    #
    t_values = (data['time_since_epoch_ms'].values - min_t) / float(1000)
    x_values = data['global_x'].values * foot_to_meter
    v_values = data['veh_v'].values * foot_to_meter
    sc = plt.scatter(t_values, x_values, s = 0.1, alpha = 0.2, 
                     c = v_values , lw = 0, cmap = plt.get_cmap('hot'))
    plt.colorbar(sc)
    plt.title('NGSIM traj. US-101 (color is speed m/s)')
    plt.ylabel('Global x position (meters)')
    plt.xlabel('t (seconds)')
    plt.savefig('US-101_traj_subset_%s.png' % period, dpi = 600)
    plt.close()

all_data = pandas.concat(all_data)

#
min_t = float(all_data['time_since_epoch_ms'].min())
min_v = float(all_data['veh_v'].min())
max_v = float(all_data['veh_v'].max())
print 'v belongs to (%f, %f)' % (min_v, max_v)
#
t_values = (all_data['time_since_epoch_ms'].values - min_t) / float(1000)
x_values = all_data['global_x'].values * foot_to_meter
v_values = all_data['veh_v'].values * foot_to_meter
sc = plt.scatter(t_values, x_values, s = 0.1, alpha = 0.2, 
                 c = v_values , lw = 0, cmap = plt.get_cmap('hot'))
plt.colorbar(sc)
plt.title('NGSIM traj. US-101 (color is speed m/s)')
plt.ylabel('Global x position (meters)')
plt.xlabel('t (seconds)')
fig = plt.gcf()
fig.set_size_inches((16,8))
plt.savefig('US-101_all_traj.png', dpi = 600)
plt.close()
