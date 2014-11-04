'''
Created on Oct 8, 2014

@author: cusgadmin
'''

import numpy as np

def kappa_cos(alpha, w, t, phi):
    coeff = alpha * np.cos(w * phi) + w * np.sin(w * phi)
    return - (coeff * np.exp(-alpha * t)
              - alpha * np.cos(w * (t + phi))
              - w * np.sin(w * (t + phi))) \
              / (alpha ** 2 + w ** 2)

def fund_cos_1_1(w, t, phi, x, lambda_1, lambda_2, tau, L):
    red_t_1 = t - x / lambda_1
    xsi_1 = np.exp( - x / (lambda_1 * tau)) * \
            np.cos(w * (red_t_1 + phi)) * \
            (red_t_1 >= 0)
    return xsi_1

def fund_cos_1_2(w, t, phi, x, lambda_1, lambda_2, tau, L):
    alpha = - lambda_2 / (tau * (lambda_1 - lambda_2))
    red_t_1 = t - x / lambda_1
    red_t_2 = t - (x - L * (lambda_1 - lambda_2) / lambda_1) / lambda_2
    xsi_2 = lambda_1 * alpha / lambda_2 * \
            (np.exp( - x / (lambda_1 * tau)) * 
             kappa_cos(alpha, w, red_t_1, phi) * (red_t_1 >= 0)
             -
             np.exp(- L / (lambda_1 * tau)) *
             kappa_cos(alpha, w, red_t_2, phi) * (red_t_2 >= 0))
    return xsi_2
    
def fund_cos_2_1(w, t, phi, x, lambda_1, lambda_2, tau, L):
    return 0

def fund_cos_2_2(w, t, phi, x, lambda_1, lambda_2, tau, L):
    red_t = t - (x - L) / lambda_2
    xsi_2 = np.cos(w * (red_t + phi)) * (red_t >= 0)
    return xsi_2

def fund_H_1_1(t, x, lambda_1, lambda_2, tau, L):
    return np.exp( - x / (lambda_1 * tau)) * (t - x / lambda_1 >= 0)
            
def fund_H_1_2(t, x, lambda_1, lambda_2, tau, L):
    alpha = - lambda_2 / (tau * (lambda_1 - lambda_2))
    red_t_1 = t - x / lambda_1
    red_t_2 = t - (x - L * (lambda_1 - lambda_2) / lambda_1) / lambda_2
    return (lambda_1 / lambda_2) * (
            np.exp( - x / (lambda_1 * tau))
            * (1.0 - np.exp( - alpha * red_t_1))
            * (red_t_1 >= 0)
            -
            np.exp( - L / (lambda_1 * tau))
            * (1.0 - np.exp( - alpha * red_t_2))
            * (red_t_2 >= 0)
            )

def fund_H_2_1(t, x, lambda_1, lambda_2, tau, L):
    return 0

def fund_H_2_2(t, x, lambda_1, lambda_2, tau, L):
    return (t - (x - L) / lambda_2 >= 0)


            
