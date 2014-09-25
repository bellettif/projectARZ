'''
Created on Sep 23, 2014

@author: francois
'''

import csv
import pandas

from matplotlib import pyplot as plt

periods = ['0750am-0805am',
           '0805am-0820am',
           '0820am-0835am']

period = periods[0]

target_dir = '/Users/francois/PATH_data/ARZ Datasets/US-101/vehicle-trajectory-data/'

file_path = '%s%s/trajectories-%s.csv' % (target_dir, period, period)

data = pandas.read_csv(file_path)

min_t = float(data['time_since_epoch_ms'].min())
max_t = float(data['time_since_epoch_ms'].max())

min_x = float(data['global_x'].min())
max_x = float(data['global_x'].max())

print min_t
print max_t

print min_x
print max_x

id_values = list(set(data[data['veh_class'] == 2]['vehicule_ID'].values))

for id_value in id_values:
    sub_selection = data[data['vehicule_ID'] == id_value]
    t_values = sub_selection['time_since_epoch_ms'].values
    t_values = (t_values - min_t) / (max_t - min_t)
    x_values = sub_selection['global_x'].values
    x_values = (x_values - min_x) / (max_x - min_x)
    plt.plot(x_values, t_values, lw = 0.1)
plt.title('NGSIM traj. US-101')
plt.xlabel('Global x position')
plt.ylabel('t')
plt.savefig('US-101_traj_subset_2.png', dpi = 600)
plt.close()




