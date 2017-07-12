# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 15:01:00 2015

@author: GF245157
"""

import matplotlib.pyplot as mpp
from numpy import linspace, meshgrid, cos, pi, sin, arange

r0 = 5
alpha = 3/5
gamma = 1-alpha**2
rmax = r0/10
psi_1 = (alpha * r0 * rmax)**2
f = lambda r,z: (r**2 - gamma * r0**2) * z**2 + alpha**2/4 * (r**2 - r0**2)**2

r = linspace(r0-2*rmax, r0+2*rmax,1000)
z = linspace(-1, 1,1000)
R, Z = meshgrid(r,z)
Psi = f(R,Z)

rc = rmax + 0.0145 #0.01388
Delta = 0.059 #0.05556
E = 0.0105
T = 0.035
omega = linspace(0, 2*pi)
R_t = r0 - Delta + (E - rc) * cos(omega) + T * cos(2*omega)
Z_t = (rc + E) * sin(omega) + T * sin(2*omega)

mpp.figure(figsize=(20,20))
mpp.hold(True)
CS = mpp.contour(R, Z, Psi, [(alpha * r0 * rmax * i/10)**2 for i in range(1, 10 +1)])
mpp.clabel(CS, inline=1, fontsize=10)
mpp.axis([4.4, 5.6, -0.6, 0.6])
mpp.plot(R_t, Z_t)
mpp.hold(False)
print(rc-E)
print(Delta-T)
rt = rc + 0.001246754
print(rt**3/(8*r0**2)+rt*Delta/(2*r0)-E**2/(2*rt)-T**2/rt)
# 0.515746754