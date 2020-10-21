#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 13:29:27 2020

@author: clement
"""

import random
import math
from scipy.stats import norm

# this script contains 3 different functions
# montecarlodigitalcall is a Monte Carlo simulation of the value of a digital call option
# digitalcalloption is the analytic value of the digital call option
# digitalcalloptionGreeks returns a table of the Greeks of the digital call option

def montecarlodigitalcall(S, T, vol, r, K, n):
    
# simulates the value of a digital call option via Monte Carlo methods

    total = 0

    for i in range(n):
    
    # run n simulations

        Y_T = math.log(S) + (r - (1/2)*(vol**2))*T + (vol * (T ** (1/2)))*random.gauss(0,1)
        S_T = math.exp(Y_T)
    
        # Y_T = log of the simulated stock price, S_T = simulated stock price
    
        if S_T > K:
            total += 1
        else:
            total += 0

        # if stock price above strike, add value of option to running total, otherwise option does not exercise

    mc_price = total/n * math.exp(-r*T)
    
    return(mc_price)

    # monte carlo solution is equal to average of the simualted option values, suitably time discounted

def digitalcalloption(S, T, vol, r, K):

# computes the analytical price of a digital call option
    
    d2 = (math.log(S/K) + (r - (1/2)*(vol**2))*T) / (vol * (T ** (1/2)))
    
    analytical_price = norm.cdf(d2) * math.exp(-r*T)
    
    return(analytical_price)

def digitalcalloptiongreeks(S, T, vol, r, K, step):

# computes analytical values of the Greeks and checks them via finite difference method
    
    d2 = (math.log(S/K) + (r - (1/2)*(vol**2))*T) / (vol * (T ** (1/2)))
    
    delta = math.exp(-r*T) * (math.exp(-(d2**2)/2)/((2*math.pi)**(1/2)))/(S * vol * (T ** (1/2)))
    gamma = (1 + d2/(vol * (T ** (1/2)))) * (-1/S) * delta
    vega = math.exp(-r*T) * (math.exp(-(d2**2)/2)/((2*math.pi)**(1/2))) * (-d2/vol - (T ** (1/2)))
    rho = (math.exp(-r*T) * (math.exp(-(d2**2)/2)/((2*math.pi)**(1/2))) * ((T ** (1/2))/vol)) - (T * math.exp(-r*T) * norm.cdf(d2))
    theta = -r * norm.cdf(d2) * math.exp(-r*T) + math.exp(-r*T) * (math.exp(-(d2**2)/2)/((2*math.pi)**(1/2))) * (-d2/(2*T) + (r - (1/2)*(vol**2)) / (vol * (T ** (1/2))))

    delta_diff = (digitalcalloption(S + step, T, vol, r, K) - digitalcalloption(S, T, vol, r, K)) * (1/step)
    gamma_diff = (digitalcalloption(S + step, T, vol, r, K) - 2*digitalcalloption(S, T, vol, r, K) + digitalcalloption(S - step, T, vol, r, K)) * (1/(step ** 2))
    vega_diff = (digitalcalloption(S, T, vol + step, r, K) - digitalcalloption(S, T, vol, r, K)) * (1/step)
    rho_diff = (digitalcalloption(S, T, vol, r + step, K) - digitalcalloption(S, T, vol, r, K)) * (1/step)
    theta_diff = (digitalcalloption(S, T + step, vol, r, K) - digitalcalloption(S, T, vol, r, K)) * (1/step)

    print('Greek\t\tAnalytic\t\tFinite Difference')
    print('Delta\t\t', round(delta,4), '\t\t', round(delta_diff,4))
    print('Gamma\t\t', round(gamma,4), '\t\t', round(gamma_diff,4))
    print('Vega\t\t\t', round(vega,4), '\t', round(vega_diff,4))
    print('Rho\t\t\t', round(rho,4), '\t', round(rho_diff,4))
    print('Theta\t\t', round(theta,4), '\t\t', round(theta_diff,4))
    
S = 100
r = 0.05
T = 1
vol = 0.1
K = 105
n = 1000000
step = 0.0001

# input variables of the B-S equation

print('The Monte Carlo simulation suggests the price of this digital call option is', round(montecarlodigitalcall(S, T, vol, r, K, n),4))
print('The analytic price of the digital call option is', round(digitalcalloption(S, T, vol, r, K),4))
digitalcalloptiongreeks(S, T, vol, r, K, step)