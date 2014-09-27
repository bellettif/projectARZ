'''
Created on Sep 26, 2014

@author: francois
'''

import numpy as np
from scipy import ndimage
import cPickle as pickle

from matplotlib import pyplot as plt

periods = ['0750am-0805am',
           '0805am-0820am',
           '0820am-0835am']

period = periods[0]

n_cut = 8 # For period 0

q_file = '../bucket_data_101/q_buckets_%s.pi' % period
v_file = '../bucket_data_101/v_buckets_%s.pi' % period
rho_file = '../bucket_data_101/rho_buckets_%s.pi' % period
x_file = '../bucket_data_101/x_buckets_%s.pi' % period
t_file = '../bucket_data_101/t_buckets_%s.pi' % period

q_buckets = pickle.load(open(q_file, 'rb'))
v_buckets = pickle.load(open(v_file, 'rb'))
x_buckets = pickle.load(open(x_file, 'rb'))
t_buckets = pickle.load(open(t_file, 'rb'))
rho_buckets = pickle.load(open(rho_file, 'rb'))

q_buckets = np.rot90(q_buckets[n_cut:-n_cut,n_cut:-n_cut])
v_buckets = np.rot90(v_buckets[n_cut:-n_cut,n_cut:-n_cut])
rho_buckets = np.rot90(rho_buckets[n_cut:-n_cut,n_cut:-n_cut])

rho_buckets = ndimage.filters.gaussian_filter(rho_buckets, 1.5)

plt.matshow(rho_buckets)
plt.show()

q_flat = np.ravel(q_buckets)
v_flat = np.ravel(v_buckets)
rho_flat = np.ravel(rho_buckets)
 
plt.scatter(rho_flat, q_flat, alpha = 0.2, lw = 0)
plt.show()


