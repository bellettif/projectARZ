'''
Created on Oct 13, 2014

@author: cusgadmin
'''

import numpy as np
from matplotlib import pyplot as plt

from time_domain_responses import fund_cos_1_2, fund_H_1_2

t_values = np.linspace(0, 3000, 80)

lambda_1 = 10.09
lambda_2 = -2.19

w = 1.0
L = 463.0

tau = 20.0

plt.plot([fund_H_1_2(t, L, lambda_1, lambda_2, tau, L) for t in t_values])
plt.show()
