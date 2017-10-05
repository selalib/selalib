# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 11:36:20 2015

@author: GF245157
"""

import scipy.constants as sc
import numpy as np
import matplotlib.pyplot as mpp
import sympy as sp
from warnings import warn

def Culham_computation(q, p, R0, B0, Ea, Ta, h, r):
    """
        Call the computation of Delta, E and T
    """
    f, g = f_g_computation(q, p, R0, B0, h, r)
    Diff_Delta = Diff_Delta_computation(f, p, R0, B0, h, r)
    Delta = primitive(Diff_Delta, 0, h)
    E, Diff_E, C_E = E_computation(q, Ea, h, r)
    T, Diff_T, C_T = T_computation(q, Ta, h, r)
    
    return Diff_Delta, Delta, E, Diff_E, T, Diff_T, f, g, C_E, C_T
#end def

def F(r0, y0):
    """
        Define the function F such that (f**2)' + 2*g' = F(r, f(r), p', B0)
        knowing that f(r)**2/r --> 0 when r --> 0
    """
    if r0 != 0:
        return -2 * y0**2 / r0
    elif y0 == 0:
        return 0
    else:
        raise TypeError('Divided by 0')
    #end if
#end def

def f_g_step(f0, g0, q1, p0, p1, R0, B0, h, r0):
    """
        Computation of the next value of f and g
    """
    r1 = r0 + h
    interm = f0**2 + 2*g0 + 2 * sc.mu_0 * p0 / B0**2 + h*F(r0, f0)
    f_i, g_i = sp.symbols('f_i,g_i')
    Results = sp.solve([2 * sc.mu_0 * p1 / B0**2 + f_i**2 + 2*g_i - interm, f_i - g_i * r1 / (R0 * q1)], [f_i, g_i], dict=True)
#    print(Results)
    sR = len(Results)
    if sR == 0:
        raise TypeError('No solution !! Either step too high or either there is a problem')
    elif sR == 1:
        sol = 0
    elif sR == 2:
        sol_f1 = Results[0][f_i]
        sol_f2 = Results[1][f_i]
        if abs(sol_f1 - f0) > abs(sol_f2 - f0):
            sol = 1
            if abs(sol_f1 - f0) < 10 * abs(sol_f2 - f0):
                warn('Solutions are close. Step may be not low enough')
            #end if
        else:
            sol = 0
            if abs(sol_f2 - f0) < 10 * abs(sol_f1 - f0):
                warn('Solutions are close. h maybe not low enough')
            #end if
        #end if
    else:
        print(sR)
        raise TypeError('Too many solutions for a quadratic equation')
    #end if
    
    return Results[sol][f_i], Results[sol][g_i]
#end def
        

def f_g_computation(q, p, R0, B0, h, r):
    """
        Computation of f and g
    """
    f = np.zeros_like(q)
    g = np.zeros_like(q)
    g[0] = 1
    
    for i in range(1, np.size(f)):
        sol_f, sol_g = f_g_step(f[i-1], g[i-1], q[i], p[i-1], p[i], R0, B0, h, r[i-1])
        f[i] = sol_f
        g[i] = sol_g
    #end for
    
    return f, g
#end def

def primitive(f, F0, h):
    """
        Compute the primitive F of f such that F(0) = F0
    """
    F = np.zeros_like(f)
    F[0] = F0
    for i in range(1, np.size(f)):
        F[i] = F[i-1] + h/2*(f[i]+f[i-1])
    #end for
    
    return F
#end def

def Diff_Delta_computation(f, p, R0, B0, h, r):
    """
        Compute Delta'
    """
    Betap = primitive(r*f**2, 0, h)
    li = 4 * sc.mu_0 / B0**2 * primitive(r * p, 0, h)
    
    Diff_Delta = Betap + li - 2 * sc.mu_0 / B0**2 * p * r**2
    Diff_Delta[0] = 0
    Diff_Delta[1::] = Diff_Delta[1::]/(R0 * r[1::] * f[1::]**2)
    
    return Diff_Delta
#end def

def E_computation(q, Ea, h, r):
    """
        Compute E
    """
    delta_E = primitive(np.log(q/q[0]) * r**3, 0, h)
    C_E = Ea / (r[-1] + 2 * delta_E[-1] / (r[-1]**3))
    E = r
    E[1::] = E[1::] + 2 * delta_E[1::] / (r[1::]**3)
    E = C_E * E
    Diff_E = 1 + 2 * np.log(q/q[0])
    Diff_E[1::] = Diff_E[1::] - 6 * delta_E[1::] / (r[1::]**4)
    Diff_E = C_E * Diff_E
    
    return E, Diff_E, C_E
#end def

def T_computation(q, Ta, h, r):
    """
        Compute T
    """
    delta_T = primitive(np.log(q/q[0]) * r**5, 0, h)
    C_T = Ta / (r[-1]**2 + 4 * delta_T[-1] / (r[-1]**4))
    T = r**2
    T[1::] = T[1::] + 4 * delta_T[1::] / (r[1::]**4)
    T = C_T * T
    Diff_T = 2 * r + 4 * r * np.log(q/q[0])
    Diff_T[1::] = Diff_T[1::] - 16 * delta_T[1::] / (r[1::]**5)
    Diff_T = C_T * Diff_T
    
    return T, Diff_T, C_T
#end def

def Culham_equilibrum(q, p, R0, B0, Ea, Ta, h, Nomega, Nr):
    """
        Build the function Delta, E and T, compute P and plot the magnetic fields
    """
    r = h*np.arange(np.size(q))
    Diff_Delta, Delta, E, Diff_E, T, Diff_T, f, g, C_E, C_T = Culham_computation(q, p, R0, B0, Ea, Ta, h, r)
    P = r**3 / (8 * R0**2) + r * Delta / (2*R0)
    P[1::] = P[1::] - E[1::]**2 / (2*r[1::]) - T[1::]**2 / r[1::]
    print('C_E assez petit ? C_E = ', C_E, ' ; C_E/epsilon_max = ', C_E*R0/r[-1])
    
    omega = np.linspace(0, 2*sc.pi, num=Nomega+1)
    
    grad_r = grad_r_square_computation(Diff_Delta, Delta, E, Diff_E, T, Diff_T, r, omega, R0)
    gr_cdot_gom = gr_cd_gom_computation(Diff_Delta, Delta, E, Diff_E, T, Diff_T, r, omega, R0)
    grad_omega = grad_omega_square_computation(Diff_Delta,  E, Diff_E, T, Diff_T, r, omega, R0)
    
    Ndr = int(np.size(q)/Nr)
    ri = range(Ndr-1, np.size(q), Ndr)
    
    mpp.figure()
    mpp.hold(True)
    for rj in ri:
        mpp.plot((- E + r - P)[rj] * np.cos(omega) - Delta[rj] + T[rj] * np.cos(2*omega), (r + E - P)[rj] * np.sin(omega) - T[rj] * np.sin(2 * omega), 
'o-')
#        print('R = ', (E - r + P)[rj] * np.cos(omega) - Delta[rj] + T[rj] * np.cos(2*omega))
#        print('Z = ', (r + E - P)[rj] * np.sin(omega) + T[rj] * np.sin(omega))
    #end for
    mpp.hold(False)
    mpp.axis("equal")
    mpp.show()
    
    return Diff_Delta, Delta, E, Diff_E, T, Diff_T, P, f, g, grad_r, gr_cdot_gom, grad_omega
#end def

def grad_r_square_computation(Diff_Delta, Delta, E, Diff_E, T, Diff_T, r, omega, R0):
    """
        Compute |grad r|**2 for every position, except r=0.
        The result is a 2D-array A: the computation of (grad r)**2 (i*h, j*domega)
        is stored in A[j, i-1]
    """
    grad_r_square = 1 - 2 * np.kron(np.cos(omega)[np.newaxis].T, Diff_Delta[1::]) + 2 * np.kron(np.cos(2 * omega)[np.newaxis].T, Diff_E[1::]) + 2 * np.kron(np.cos(3 * omega)[np.newaxis].T, Diff_T[1::]) + np.kron(np.ones_like(omega)[np.newaxis].T, ((r/R0)**2 + 2*Delta/R0 + r*Diff_Delta/R0 + 3/2*Diff_Delta**2 + 3/2*Diff_E**2 + 3/2*Diff_T**2)[1::] + 3/2*(E[1::]/r[1::])**2 - E[1::]*Diff_E[1::]/r[1::] + 4*(T[1::]/r[1::])**2 - 2*T[1::]*Diff_T[1::]/r[1::])
    
    return grad_r_square

def DD_Delta_computation(Diff_Delta, p, f, r, R0, B0, h):
    """
        Compute the second derivative of Delta
    """
    DD_Delta = 1/R0 + 0*Diff_Delta
    DD_Delta[1::][:-1:] = DD_Delta[1::][:-1:] - (1/r[1::][:-1:] + (f[2::] - f[0:-2:])/f[1::][:-1:]) * Diff_Delta[1::][:-1:] - 2*r[1::][:-1:]/R0*sc.mu_0/(B0**2*f[1::][:-1:]**2)*(p[2::] - p[0:-2:])
    DD_Delta[-1] = DD_Delta[-1] - (1/r[-1] + (f[-1] - f[-2])/(2*f[-1])) * Diff_Delta[-1] - r[-1]/R0*sc.mu_0/(B0**2*f[-1]**2)*(p[-1] - p[-2])
    
    return DD_Delta
    
def DD_E_computaion(Diff_E, E, r, f):
    """
        Compute the second derivative of E
    """
    DD_E = 0*Diff_E
    DD_E[1::][:-1:] = - (1/r[1::][:-1:] + (f[2::] - f[0:-2:])/f[1::][:-1:]) * Diff_E[1::][:-1:] + 3*E[1::][:-1:]/r[1::][:-1:]**2
    DD_E[-1] =  - (1/r[-1] + (f[-1] - f[-2])/f[-1]) * Diff_E[-1] + 3*E[-1]/r[-1]**2
    
    return DD_E

def DD_T_computation(Diff_T, T, r, f):
    """
        Compute the second derivative of T
    """
    DD_T = 0*Diff_T
    DD_T[1::][:-1:] = - (1/r[1::][:-1:] + (f[2::] - f[0:-2:])/f[1::][:-1:]) * Diff_T[1::][:-1:] + 8*T[1::][:-1:]/r[1::][:-1:]**2
    DD_T[-1] =  - (1/r[-1] + (f[-1] - f[-2])/f[-1]) * Diff_T[-1] - 3*T[-1]/r[-1]**2
    
    return DD_T

def gr_cd_gom_computation(Diff_Delta, Delta, E, Diff_E, T, Diff_T, r, omega, R0):
    """
        Compute r * (grad r . grad omega) for every position, except r=0.
        The result is a 2D-array A: the computation of r * (grad r . grad omega) (i*h, j*domega)
        is stored in A[j, i-1]
    """
    gr_cd_gom = - np.kron(np.sin(omega)[np.newaxis].T, Diff_Delta[1::]) - np.kron(np.sin(2*omega)[np.newaxis].T, (Diff_E[1::]+E[1::]/r[1::])) + np.kron(np.sin(3*omega)[np.newaxis].T, (Diff_T[1::]+2*T[1::]/r[1::]))
    
    return gr_cd_gom

def grad_omega_square_computation(Diff_Delta,  E, Diff_E, T, Diff_T, r, omega, R0):
    """
        Compute r**2 * (grad omega)**2 for every position, except r=0.
        The result is a 2D-array A: the computation of r**2 * (grad omega)**2 (i*h, j*domega)
        is stored in A[j, i-1]
    """
    grad_omega_square = 1 - np.kron(np.cos(2*omega)[np.newaxis].T, 2*E[1::]/r[1::]) + np.kron(np.cos(3*omega)[np.newaxis].T, 4*T[1::]/r[1::]) + np.kron(np.ones_like(omega)[np.newaxis].T, (Diff_Delta**2/2 - r*Diff_Delta/R0 + r**2/(4*R0**2) + Diff_T**2/2)[1::] + (E[1::]/r[1::] + Diff_E[1::])**2/2 + 2*T[1::]*Diff_T[1::]/r[1::] + 4*T[1::]**2/r[1::]**2)
    
    return grad_omega_square

#a = 1
#R0 = 5
#N = 100
#B0 = 1
#q0 = 0.8
#qa = 0.7
#alpha = 1
#p0 = 10**5
#Nomega = 50
#Nr = 10
#pa = p0/10
#h = a/N
#r = h*np.arange((N+1))
#q = q0 + (qa - q0)*(r/a)**2
#p = (p0 - pa)*(1-(r/a)**2)**alpha + pa
#Delta, E, T, P = Culham_equilibrum(q, p, R0, B0, 0.25, 0.1, h, Nomega, Nr)
#print(Delta)
