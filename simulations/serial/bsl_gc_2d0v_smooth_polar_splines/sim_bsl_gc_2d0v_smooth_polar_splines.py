import sys
import h5py
import numpy as np

np.set_printoptions( precision=15 )

# Read folder name from standard input
fold_name = sys.argv[1]

# Open HDF5 file
file_name = fold_name+'sim_bsl_gc_2d0v_smooth_polar_splines.h5'
h5 = h5py.File( file_name, mode='r' )

# Load x1 and x2 mesh
x1 = h5['x1'].value
x2 = h5['x2'].value

# Load attributes
n1 = h5.attrs['n1']
n2 = h5.attrs['n2']
nx = h5.attrs['nx1']
ny = h5.attrs['nx2']
nc = h5.attrs['nc']
dt = h5.attrs['time_step']
Ni = h5.attrs['iterations']
df = h5.attrs['diag_freq']

Nd = (Ni-1) // df + 1
tt = np.ndarray( Nd, dtype=int )
for i in range(0,Ni,df):
    j = i // df
    tt[j] = i

# Load data
rho = dict()
phi = dict()
Ex  = dict()
Ey  = dict()
Ex_cart = dict()
Ey_cart = dict()
for t in tt:
    rho    [t] = h5['rho_'    +str(t)].value
    phi    [t] = h5['phi_'    +str(t)].value
    Ex     [t] = h5['Ex_'     +str(t)].value
    Ey     [t] = h5['Ey_'     +str(t)].value
    Ex_cart[t] = h5['Ex_cart_'+str(t)].value
    Ey_cart[t] = h5['Ey_cart_'+str(t)].value

if ( nc > 0 ):
   point_charges = dict()
   for t in range(Ni):
       point_charges[t] = h5['point_charges_'+str(t)].value

# Load equilibrium density
rho_eq = h5['rho_eq'].value

# Load Cartesian grids
x1_cart = h5['x1_cart'].value
x2_cart = h5['x2_cart'].value

# Close HDF5 file
h5.close()

# magnitude of electric field
Em = dict()
for t in tt:
    Em[t] = np.sqrt( Ex[t]**2 + Ey[t]**2 )

# to plot fine grid
nr = n1 // 8

#-------------------------------------------------------------------------------
# Plot time evolution of RHO
#-------------------------------------------------------------------------------

# find minimum and maximum over all times
def minmax( Q ):
    min_Q = Q[0].min()
    max_Q = Q[0].max()
    for t in tt:
        if ( Q[t].min() < min_Q ):
           min_Q = Q[t].min()
        if ( Q[t].max() > max_Q ):
           max_Q = Q[t].max()
    return [ min_Q, max_Q ]

min_phi, max_phi = minmax( phi )
min_Ex , max_Ex  = minmax( Ex )
min_Ey , max_Ey  = minmax( Ey )
min_Em , max_Em  = minmax( Em )

min_rho, max_rho = minmax( rho )
# if max_rho = 0.0, contour levels need to be increasing
if ( max_rho == min_rho ):
   max_rho = 1.0

min_rho_eq = rho_eq.min()
max_rho_eq = rho_eq.max()
# if max_rho_eq = 0.0, contour levels need to be increasing
if ( max_rho_eq == min_rho_eq ):
   max_rho_eq = 1.0

min_delta_rho = ( rho[0] - rho_eq ).min()
max_delta_rho = ( rho[0] - rho_eq ).max()
for t in tt:
    delta_rho = rho[t] - rho_eq
    if ( delta_rho.min() < min_delta_rho ):
       min_delta_rho = delta_rho.min()
    if ( delta_rho.max() > max_delta_rho ):
       max_delta_rho = delta_rho.max()
# if max_delta_rho = 0.0, contour levels need to be increasing
if ( max_delta_rho == min_delta_rho ):
   max_delta_rho = 1.0

#-------------------------------------------------------------------------------
# Scalar diagnostics
#-------------------------------------------------------------------------------

# Analytical growth rate
l = 9
gamma = 0.179630959411
coeff = 3.5e-08

# Load data from file
# y1 = mass, y2 = energy, y3 = L2-norm of phi
file_name = fold_name+'scalar_diagnostics.dat'
xt = np.loadtxt( file_name, usecols=[0] )
y1 = np.loadtxt( file_name, usecols=[1] )
y2 = np.loadtxt( file_name, usecols=[2] )
y3 = np.loadtxt( file_name, usecols=[3] )

ya = coeff * np.exp( gamma * xt )

## Limit t in some range
#mask = ( xt >= 0. )
#xt = xt[mask]
#y1 = y1[mask]
#ya = ya[mask]

#-------------------------------------------------------------------------------
# Vorticity
#-------------------------------------------------------------------------------

Ax = dict()
Ay = dict()
for t in tt:
    Ax[t] = - Ey_cart[t].transpose()
    Ay[t] =   Ex_cart[t].transpose()

dx = 2.0 / nx
dy = 2.0 / ny

# shifted grids
x_grid = np.ndarray( nx-1 )
y_grid = np.ndarray( ny-1 )

for jx in range(nx-1):
    x_grid[jx] = ( x1_cart[jx] + x1_cart[jx+1] ) / 2.
for jy in range(ny-1):
    y_grid[jy] = ( x2_cart[jy] + x2_cart[jy+1] ) / 2.

x_meshgrid, y_meshgrid = np.meshgrid( x_grid, y_grid, indexing='ij' )

dAx_dx = dict()
dAx_dy = dict()
dAy_dx = dict()
dAy_dy = dict()
for t in tt:
    dAx_dx[t] = np.ndarray( [ nx-1, ny-1 ] )
    dAx_dy[t] = np.ndarray( [ nx-1, ny-1 ] )
    dAy_dx[t] = np.ndarray( [ nx-1, ny-1 ] )
    dAy_dy[t] = np.ndarray( [ nx-1, ny-1 ] )

for t in tt:
    # dAx/dx
    dAx_dx[t][:,:] = ( Ax[t][1:,:-1] - Ax[t][:-1,:-1] \
                     + Ax[t][1:,1: ] - Ax[t][:-1,1: ] ) / ( 2.*dx )
    # dAx/dy
    dAx_dy[t][:,:] = ( Ax[t][:-1,1:] - Ax[t][:-1,:-1] \
                     + Ax[t][1: ,1:] - Ax[t][1: ,:-1] ) / ( 2.*dy )
    # dAy/dx
    dAy_dx[t][:,:] = ( Ay[t][1:,:-1] - Ay[t][:-1,:-1] \
                     + Ay[t][1:,1: ] - Ay[t][:-1,1: ] ) / ( 2.*dx )
    # dAy/dy
    dAy_dy[t][:,:] = ( Ay[t][:-1,1:] - Ay[t][:-1,:-1] \
                     + Ay[t][1: ,1:] - Ay[t][1: ,:-1] ) / ( 2.*dy )

# expansion dAx/dx + dAy/dy
expansion = dict()
for t in tt:
    expansion[t] = dAx_dx[t] + dAy_dy[t]

# vorticity dAy/dx - dAx/dy
vorticity = dict()
for t in tt:
    vorticity[t] = dAy_dx[t] - dAx_dy[t]

# shear rate (eigenvalues of traceless symmetric part of Jacobian of advection field)
shearrate = dict()
for t in tt:
    shearrate[t] = 0.5 * np.sqrt( ( dAx_dx[t] - dAy_dy[t] )**2 + ( dAx_dy[t] + dAy_dx[t] )**2 )

min_expansion, max_expansion = minmax( expansion )
min_vorticity, max_vorticity = minmax( vorticity )
min_shearrate, max_shearrate = minmax( shearrate )

#-------------------------------------------------------------------------------
# stream lines in inertial and rotating frame (works for circle only)
#-------------------------------------------------------------------------------

t = 0

x1_mesh, x2_mesh = np.meshgrid( x1_cart, x2_cart )

# advection field in rotating frame
omega = 0.3332
Ax_rot = - Ey_cart[t] + omega*x2_mesh
Ay_rot =   Ex_cart[t] - omega*x1_mesh

# set zero outside domain
for ix in range(nx):
    for iy in range(ny):
        if ( np.sqrt(x1_cart[ix]**2+x2_cart[iy]**2) > 1.0 ):
           Ex_cart[t][ix,iy] = 0.0
           Ey_cart[t][ix,iy] = 0.0
           Ax_rot[ix,iy] = 0.0
           Ay_rot[ix,iy] = 0.0
