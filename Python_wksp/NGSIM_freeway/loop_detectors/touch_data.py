'''
Created on Sep 26, 2014

@author: francois
'''

import numpy as np
from matplotlib import pyplot as plt
import pandas

columns = ['V_avg', 'Occ_avg', 'Sp_avg']

file_path = '/Users/cusgadmin/projectARZ/US-101/detector-data/detector-data.csv'

data = pandas.read_csv(file_path)

print data['V_avg'] - (data['V1'] + data['V2'] + data['V3'] + data['V4'] + data['V5'])
print data['Occ_avg'] - 0.2 * (data['Occ1'] + data['Occ2'] + data['Occ3'] + data['Occ4'] + data['Occ5'])
print data['Sp_avg'] - 0.2 * (data['Sp1'] + data['Sp2'] + data['Sp3'] + data['Sp4'] + data['Sp5'])

v_values = data['V_avg'].values
occ_values = data['Occ_avg'].values
sp_values = data['Sp_avg'].values * 1609.34 / 3600.0

plt.scatter(occ_values, v_values, alpha = 0.2, lw= 0)
plt.show()


