'''
Created on Sep 25, 2014

@author: francois
'''

import csv
import pandas

from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

periods = ['0750am-0805am',
           '0805am-0820am',
           '0820am-0835am']

target_dir = '/Users/cusgadmin/projectARZ/US-101/vehicle-trajectory-data/'

foot_to_meter = 0.3048

for period in periods:
    file_path = '%s%s/trajectories-%s.csv' % (target_dir, period, period)
    data = pandas.read_csv(file_path)
    min_t = float(data['time_since_epoch_ms'].min())
    min_x = float(data['global_x'].min())
    max_x = float(data['global_x'].max())
    t_values = (data['time_since_epoch_ms'].values - min_t) / float(1000)
    x_values = data['global_x'].values * foot_to_meter
    plt.hist2d(x_values, t_values, bins = 200)
    plt.colorbar()
    plt.title('Number of data points')
    plt.xlabel('x (meters)')
    plt.ylabel('t (seconds)')
    plt.savefig('US-101_occupancy_%s.png' % period, dpi = 600)
    plt.close()