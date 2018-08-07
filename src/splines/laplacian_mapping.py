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

if ( mapping == 'circle' ):
   phi = ( 1 - s**2 ) * sp.cos( 2*sp.pi*x ) * sp.sin( 2*sp.pi*y )
else:
   phi = s**2 * ( 1 - s**2 ) * sp.cos(t)

dphi_ds = phi.diff( s, 1 )
dphi_dt = phi.diff( t, 1 )

l_ss = det_J * g__ss * dphi_ds
l_st = det_J * g__st * dphi_dt
l_ts = det_J * g__ts * dphi_ds
l_tt = det_J * g__tt * dphi_dt

lap_phi = l_ss.diff( s, 1 ) + l_st.diff( s, 1 ) + l_ts.diff( t, 1 ) + l_tt.diff( t, 1 )

lap_phi = sp.simplify( lap_phi / det_J )
#lap_phi = sp.collect ( lap_phi, [ s, k, d ] )

rho = - lap_phi
