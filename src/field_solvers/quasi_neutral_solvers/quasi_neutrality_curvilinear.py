# Compute rho=-div(A*grad(phi))+c*phi in curvilinear coordinates
# How to run the script:
# >> ipython3 qtconsole
# >> run quasi_neutrality_curvilinear.py circle

import sys
import sympy as sp
from sympy import *

sp.init_session()
sp.init_printing()

mapping = sys.argv[1]

s = sp.symbols('s')
t = sp.symbols('t')

if ( mapping == 'circle' ):
   x = s*sp.cos(t)
   y = s*sp.sin(t)
elif ( mapping == 'target' ):
   xp = sp.symbols('x_p')
   yp = sp.symbols('y_p')
   k = sp.symbols('kappa')
   d = sp.symbols('Delta')
   x = xp + s*sp.cos(t) - k*s*sp.cos(t) - d*s**2
   y = yp + s*sp.sin(t) + k*s*sp.sin(t)
elif ( mapping == 'czarny' ):
   yp = sp.symbols('y_p')
   inv_asp_ratio = sp.symbols('epsilon')
   e  = sp.symbols('e')
   xi = sp.symbols('xi')
   x = ( 1 - sp.sqrt( 1 + inv_asp_ratio * ( inv_asp_ratio + 2*s*sp.cos(t)) ) ) / inv_asp_ratio
   y = yp + e * xi * s * sin(t) / ( 1 + inv_asp_ratio * x )

g_ss = x.diff( s, 1 )**2 + y.diff( s, 1 )**2
g_st = x.diff( s, 1 ) * x.diff( t, 1 ) + y.diff( s, 1 ) * y.diff( t, 1 )
g_ts = x.diff( t, 1 ) * x.diff( s, 1 ) + y.diff( t, 1 ) * y.diff( s, 1 )
g_tt = x.diff( t, 1 )**2 + y.diff( t, 1 )**2

g_ss = sp.simplify( g_ss )
g_st = sp.simplify( g_st )
g_ts = sp.simplify( g_ts )
g_tt = sp.simplify( g_tt )

det_g = g_ss * g_tt - g_st * g_ts
det_J = sp.sqrt( det_g )

det_g = sp.simplify( det_g )
det_J = sp.simplify( det_J )

g__ss =   g_tt / det_g
g__st = - g_st / det_g
g__ts = - g_ts / det_g
g__tt =   g_ss / det_g

phi = s**2 * ( 1 - s**2 ) * sp.cos(t)
A = sp.exp( - sp.tanh( ( s - 0.5 ) / 0.1 ) )
c = 1.0 / sp.exp( - sp.tanh( ( s - 0.5 ) / 0.2 ) )

dphi_ds = phi.diff( s, 1 )
dphi_dt = phi.diff( t, 1 )

qn_ss = det_J * g__ss * dphi_ds
qn_st = det_J * g__st * dphi_dt
qn_ts = det_J * g__ts * dphi_ds
qn_tt = det_J * g__tt * dphi_dt

qn_ss = A * qn_ss
qn_st = A * qn_st
qn_ts = A * qn_ts
qn_tt = A * qn_tt

qn_phi = qn_ss.diff( s, 1 ) + qn_st.diff( s, 1 ) + qn_ts.diff( t, 1 ) + qn_tt.diff( t, 1 )

qn_phi = sp.simplify( qn_phi / det_J )
#qn_phi = sp.collect ( qn_phi, [ s, k, d ] )

rho = - qn_phi + c * phi
