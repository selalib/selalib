# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 09:10:21 2015

@author: GF245157
"""

from Cequilibrum import *
from scipy.constants import mu_0, pi
import numpy as np
import matplotlib.pyplot as mpp

A = 1.
alpha = 3./5
gamma = 1 - alpha**2
R0 = 5
a = 0.515746754
N = 100
h = a/N
B0 = 1.
Nomega = 50
Nr = 10

r = h*np.arange((N+1))
p = 10**6 - (A*alpha*R0*r)**2/(2*mu_0*(1+alpha**2))
q = (1+alpha**2)/(A*alpha**2*R0**2)*np.sqrt(B0**2 + (1-alpha**2)*alpha**2/((1+alpha**2)**2)*(A*R0*r)**2)
Diff_Delta_C, Delta_C, E_C, Diff_E_C, T_C, Diff_T_C, P_C, f_C, g_C, grad_r_C, gr_cdot_gomega_C, grad_omega_C = Culham_equilibrum(q, p, R0, B0, 0.0105, 0.035, h, Nomega, Nr)

r0 = 5.
rmax = r0/10
psi_1 = (alpha * r0 * rmax)**2
f = lambda r,z: (r**2 - gamma * r0**2) * z**2 + alpha**2/4 * (r**2 - r0**2)**2

r = np.linspace(r0-2*rmax, r0+2*rmax,1000)
z = np.linspace(-1, 1,1000)
R, Z = np.meshgrid(r,z)
Psi = f(R,Z)

mpp.figure(figsize=(10,10))
mpp.hold(True)
CS = mpp.contour(R, Z, Psi, [(alpha * r0 * rmax * i/Nr)**2 for i in range(1, Nr +1)])
mpp.clabel(CS, inline=1, fontsize=10)


r = h*np.arange((N+1))
omega = np.linspace(0, 2*pi)
#mpp.plot(R0 + (E_C[-1] - r[-1] + P_C[-1]) * np.cos(omega) - Delta_C[-1] + T_C[-1] * np.cos(2*omega), (r[-1] + E_C[-1] - P_C[-1]) * np.sin(omega) + T_C[-1] * np.sin(2 * omega), 'o-')
Ndr = int(N/Nr)
ri = range(Ndr, N+1, Ndr)

for rj in ri:
    mpp.plot(R0 + (E_C[rj] - r[rj] + P_C[rj]) * np.cos(omega) - Delta_C[rj] + T_C[rj] * np.cos(2*omega), (r[rj] + E_C[rj] - P_C[rj]) * np.sin(omega) + T_C[rj] * np.sin(2*omega), 'o-')
#        print('R = ', (E - r + P)[rj] * np.cos(omega) - Delta[rj] + T[rj] * np.cos(2*omega))
#        print('Z = ', (r + E - P)[rj] * np.sin(omega) + T[rj] * np.sin(2*omega))
    #end for
mpp.hold(False)
mpp.axis("equal")
mpp.axis([r0-0.6, r0+0.6, -0.6, 0.6])
mpp.show()
    
