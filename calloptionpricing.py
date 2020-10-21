#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 20:47:43 2020

@author: clement
"""
import random
import math
from scipy.stats import norm

# this script contains 3 different functions
# montecarlocall is a Monte Carlo simulation of the value of a call option
# calloption is the analytic value of the call option
# calloptionGreeks returns a table of the Greeks of the call option

def montecarlocall(S, T, vol, r, K, n):
    
# simulates the value of a call option via Monte Carlo methods

    total = 0

    for i in range(n):
    
    # run n simulations

        Y_T = math.log(S) + (r - (1/2)*(vol**2))*T + (vol * (T ** (1/2)))*random.gauss(0,1)
        S_T = math.exp(Y_T)
    
        # Y_T = log of the simulated stock price, S_T = simulated stock price
    
        if S_T > K:
            total += (S_T - K)
        else:
            total += 0

        # if stock price above strike, add value of option to running total, otherwise option does not exercise

    mc_price = total/n * math.exp(-r*T)
    
    return(mc_price)

    # monte carlo solution is equal to average of the simualted option values, suitably time discounted

def calloption(S, T, vol, r, K):

# computes the analytical price of a call option
    
    d1 = (math.log(S/K) + (r + (1/2)*(vol**2))*T) / (vol * (T ** (1/2)))
    d2 = (math.log(S/K) + (r - (1/2)*(vol**2))*T) / (vol * (T ** (1/2)))
    
    analytical_price = S * norm.cdf(d1) - K * math.exp(-r*T) * norm.cdf(d2)
    
    return(analytical_price)

def calloptiongreeks(S, T, vol, r, K, step):

# computes analytical values of the Greeks and checks them via finite difference method
    
    d1 = (math.log(S/K) + (r + (1/2)*(vol**2))*T) / (vol * (T ** (1/2)))
    d2 = (math.log(S/K) + (r - (1/2)*(vol**2))*T) / (vol * (T ** (1/2)))
    
    delta = norm.cdf(d1)
    gamma = (math.exp(-(d1**2)/2)/((2*math.pi)**(1/2)))/(S * vol * (T ** (1/2)))
    vega = S * (T **(1/2)) * (math.exp(-(d1**2)/2)/((2*math.pi)**(1/2)))
    rho = K * T * math.exp(-r*T) * norm.cdf(d2)
    theta = (S * (math.exp(-(d1**2)/2)/((2*math.pi)**(1/2))) * vol) / (2 * (T ** (1/2))) + (r * K * math.exp(-r*T) * norm.cdf(d2))

    delta_diff = (calloption(S + step, T, vol, r, K) - calloption(S, T, vol, r, K)) * (1/step)
    gamma_diff = (calloption(S + step, T, vol, r, K) - 2*calloption(S, T, vol, r, K) + calloption(S - step, T, vol, r, K)) * (1/(step ** 2))
    vega_diff = (calloption(S, T, vol + step, r, K) - calloption(S, T, vol, r, K)) * (1/step)
    rho_diff = (calloption(S, T, vol, r + step, K) - calloption(S, T, vol, r, K)) * (1/step)
    theta_diff = (calloption(S, T + step, vol, r, K) - calloption(S, T, vol, r, K)) * (1/step)

    print('Greek\t\tAnalytic\t\tFinite Difference')
    print('Delta\t\t', round(delta,4), '\t\t', round(delta_diff,4))
    print('Gamma\t\t', round(gamma,4), '\t\t', round(gamma_diff,4))
    print('Vega\t\t\t', round(vega,4), '\t', round(vega_diff,4))
    print('Rho\t\t\t', round(rho,4), '\t', round(rho_diff,4))
    print('Theta\t\t', round(theta,4), '\t\t', round(theta_diff,4))
    
S = 105
r = 0.05
T = 1
vol = 0.1
K = 105
n = 1000000
step = 0.0001

# input variables of the B-S equation

print('The Monte Carlo simulation suggests the price of this call option is', round(montecarlocall(S, T, vol, r, K, n),4))
print('The analytic price of the call option is', round(calloption(S, T, vol, r, K),4))
calloptiongreeks(S, T, vol, r, K, step)