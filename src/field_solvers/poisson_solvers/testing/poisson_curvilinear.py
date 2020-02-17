# Compute rho=-div(grad(phi)) in curvilinear coordinates
# How to run the script:
# >> ipython3 qtconsole
# >> run poisson_curvilinear.py circle

import sys
import h5py
import numpy as np
import sympy as sp
from sympy import *

sp.init_session()
sp.init_printing()

mapping = sys.argv[1]

s = sp.symbols( 's', real=True, positive=True )
t = sp.symbols( 'theta', real=True, positive=True )

if ( mapping == 'circle' ):
   x = s*sp.cos(t)
   y = s*sp.sin(t)
elif ( mapping == 'target' ):
   k = sp.symbols( 'kappa', real=True, positive=True )
   d = sp.symbols( 'Delta', real=True, positive=True )
   xp = 0.0
   yp = 0.0
   x = xp + s*sp.cos(t) - k*s*sp.cos(t) - d*s**2
   y = yp + s*sp.sin(t) + k*s*sp.sin(t)
elif ( mapping == 'czarny' ):
    print( 'Not implemented' )
    sys.exit()
#   yp = 0.0
#   epsil = sp.symbols( 'epsilon', real=True, positive=True )
#   ellip = sp.symbols( 'e', real=True, positive=True )
#   xi    = sp.symbols( 'xi', real=True, positive=True )
#   x = ( 1 - sp.sqrt( 1 + epsil * ( epsil + 2*s*sp.cos(t)) ) ) / epsil
#   y = yp + ellip * xi * s * sin(t) / ( 1 + epsil * x )

g_ss = x.diff( s )**2 + y.diff( s )**2
g_st = x.diff( s ) * x.diff( t ) + y.diff( s ) * y.diff( t )
g_ts = x.diff( t ) * x.diff( s ) + y.diff( t ) * y.diff( s )
g_tt = x.diff( t )**2 + y.diff( t )**2

g = sp.Matrix( [ [ g_ss, g_st ], [ g_ts, g_tt ] ] )
g = sp.simplify( g )

det_g = g.det()
det_g = sp.simplify( det_g )

G = g**(-1) * det_g
G = sp.simplify( G )

det_J = sp.sqrt( det_g )
det_J = sp.simplify( det_J )

phi = ( 1 - s**2 ) * sp.cos( 2*sp.pi*x ) * sp.sin( 2*sp.pi*y )

dphi = sp.Matrix( [ phi.diff( s ), phi.diff( t ) ] )

DET_g = det_g / s**2
DET_J = det_J / s
DPHI_DT = sp.cancel( dphi[1] / s )

l_ss = det_J * G[0,0] / det_g * dphi[0]
l_st = det_J * G[0,1] / det_g * dphi[1]
l_ts = det_J * G[1,0] / det_g * dphi[0]
l_tt = DET_J * G[1,1] / DET_g * DPHI_DT

lap_temp = l_ss.diff( s ) + l_st.diff( s ) + l_ts.diff( t ) + l_tt.diff( t )
lap_glob = lap_temp / det_J

if ( sp.simplify( lap_temp.subs( s, 0 ) ) == 0 ):
    # do not divide by det_J because numerator is 0
    lap_pole = lap_temp.subs( s, 0 ) 
    lap_pole = sp.simplify( lap_pole )
else:
    print( 'Error in evaluation of Laplacian at s=0' )
    sys.exit()

rho_glob = - lap_glob
rho_pole = - lap_pole

if ( mapping == 'circle' ):
    RHO_glob = lambdify( [ s, t ], rho_glob )
    RHO_pole = lambdify( [ s, t ], rho_pole )
elif ( mapping == 'target' ):
    kv = 0.3
    dv = 0.2
    RHO_glob = lambdify( [ s, t ], ( rho_glob.subs( k, kv ) ).subs( d, dv ) )
    RHO_pole = lambdify( [ s, t ], ( rho_pole.subs( k, kv ) ).subs( d, dv ) )
elif ( mapping == 'czarny' ):
    epsilv = 0.3
    ellipv = 1.4
    RHO_glob = lambdify( [ s, t ], ( rho_glob.subs( epsil, epsilv ) ).subs( ellip, ellipv ) )
    RHO_pole = lambdify( [ s, t ], ( rho_pole.subs( epsil, epsilv ) ).subs( ellip, ellipv ) )

# Open file
file_name = 'test_poisson_2d_fem_sps_stencil_grid.h5'
h5 = h5py.File( file_name, mode='r' )

# Load data
eta1_grid = h5['eta1_grid'].value
eta2_grid = h5['eta2_grid'].value

# Close file
h5.close()

n1 = eta1_grid.size
n2 = eta2_grid.size

eta1_meshgrid, eta2_meshgrid = np.meshgrid( eta1_grid, eta2_grid )

rho_grid = np.zeros( ( n2, n1 ) )

rho_grid[:,0 ] = RHO_pole( eta1_meshgrid[:,0 ], eta2_meshgrid[:,0 ] )
rho_grid[:,1:] = RHO_glob( eta1_meshgrid[:,1:], eta2_meshgrid[:,1:] )

# Open file
file_name = 'test_poisson_2d_fem_sps_stencil_rho.h5'
h5 = h5py.File( file_name, mode='w' )

# Write data
h5.create_dataset( 'rho_grid', data=rho_grid )

# Close file
h5.close()
