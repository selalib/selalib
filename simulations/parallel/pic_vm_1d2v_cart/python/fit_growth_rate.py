#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 18:23:58 2019

@author: Katharina Kormann
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema

# Parameter dependent function
def energy_par(x,c,d):
    return np.exp(2*c*x)*d#((1+e*np.cos(x*a+b)))*np.exp(2*c*x)*d

# TODO: read in the data

# set the range wher the growth rate shall be fitted
ran = range(1500,4000)
xdata = data[ran,0]
ydata = data[ran,3]


popt, pcov = curve_fit(energy_par, xdata,ydata, p0=[0.02784,1e-8],absolute_sigma=True)#,[1.7*0.5,0.0,0.02784, 1E-8,1.0])

#plt.plot(xdata, ydata, xdata, energy_par(xdata, *popt))#,xdata, energy_par(xdata, 1.7*0.5,0.0,0.02784, 1E-4,1.0))
plt.semilogy(xdata, ydata, xdata, energy_par(xdata, *popt))

ind=argrelextrema(ydata, np.greater)
plt.semilogy(xdata,ydata,xdata[ind],ydata[ind],'*')
plt.show()
print(popt)

xxdata = xdata[ind];
yydata = ydata[ind];

popt2, pcov2 = curve_fit(energy_par, xxdata,yydata, p0=[0.02784,1e-8],absolute_sigma=True)
plt.semilogy(xdata, ydata, xdata, energy_par(xdata, *popt2))
plt.show()
print(popt)