import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

# Read folder name from standard input
fold_name = sys.argv[1]

# Open HDF5 file
file_name = fold_name+'sim_bsl_gc_2d0v_smooth_polar_splines.h5'
h5 = h5py.File( file_name, mode='r' )

# Load x1 and x2 mesh
x1 = h5['x1'].value
x2 = h5['x2'].value

# Load Jacobian
Jd = h5['jacobian'].value

# Load control points along x1 and x2
c_x1 = h5['c_x1'].value
c_x2 = h5['c_x2'].value

# Load other attributes

# degrees of freedom
n1 = h5.attrs['n1']
n2 = h5.attrs['n2']
# spline degrees
p1 = h5.attrs['p1']
p2 = h5.attrs['p2']
# interpolation points
ntau1 = h5.attrs['ntau1']
ntau2 = h5.attrs['ntau2']
# finite elements
Nk1 = h5.attrs['Nk1']
Nk2 = h5.attrs['Nk2']
# quadrature points
Nq1 = h5.attrs['Nq1']
Nq2 = h5.attrs['Nq2']
# Cartesian grid points
nx = h5.attrs['nx1']
ny = h5.attrs['nx2']
# time stepping
dt = h5.attrs['time_step']
Ni = h5.attrs['iterations']
df = h5.attrs['diag_freq']
# point charges
nc = h5.attrs['nc']

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
   for t in tt:
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

# color map for plots
cmap = 'Greys'

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

min_rho_eq, max_rho_eq = minmax( rho_eq )
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

# plot equilibrium density
fg = plt.figure(figsize=[9.0,9.0])
ax = fg.add_subplot(111)
cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )
clevels = np.linspace( min_rho_eq, max_rho_eq, 101 )
im = ax.contourf( x1, x2, rho_eq, clevels, cmap=cmap )
fg.colorbar( im, cax=cax )
ax.set_title( r'Equilibrium density $\rho_{eq}$' )
ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=0.5 )
ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
ax.set_xlabel( r'$x$' )
ax.set_ylabel( r'$y$', rotation=0 )
ax.set_aspect( 'equal' )
fg.tight_layout()
fg.show()

## save figures
#t  = 0
#fg = plt.figure()
#ax = fg.add_subplot(111)
#cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )
#
#clevels = np.linspace( min_rho, max_rho, 101 )
## contour plot
#im = ax.contourf( x1, x2, rho[t], clevels, cmap=cmap )
#fg.colorbar( im, cax=cax )
## grid
#ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
#ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=0.5 )
#ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
## style
#fs = 10
#ax.set_xlabel( r'$x$', fontsize=fs )
#ax.set_ylabel( r'$y$', fontsize=fs, rotation=0 )
#ax.set_title( r'Density $\rho$ at $t = %g$' %(t*dt) )
#ax.set_aspect( 'equal' )
#fg.tight_layout()
#fg.show()
#fg_name = './rho_t=%g.pdf' %(i*dt)
#fg.savefig( fg_name, dpi=300 )

# Animated plot
def plot_time_evolution( arg ):

    fg = plt.figure(figsize=[9.0,9.0])
    ax = fg.add_subplot(111)
    cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )

    for t in tt:
        ax.clear()
        cax.clear()
        if ( arg == 'rho' ):
           clevels = np.linspace( min_rho, max_rho, 101 )
           im = ax.contourf( x1, x2, rho[t], clevels, cmap=cmap )
           if ( nc > 0 ):
              for ic in range(nc):
                  xc = point_charges[t][ic,0]
                  yc = point_charges[t][ic,1]
                  ax.plot( xc, yc, marker='o', color='k', markersize=5. )
           ax.set_title( r'Density $\rho$ at $t = %g$' %(t*dt) )
        if ( arg == 'delta_rho' ):
           clevels = np.linspace( min_delta_rho, max_delta_rho, 101 )
           im = ax.contourf( x1, x2, rho[t]-rho_eq, clevels, cmap=cmap )
           ax.set_title( r'Density $\rho$ at $t = %g$' %(t*dt) )
        if ( arg == 'phi' ):
           clevels = np.linspace( min_phi, max_phi, 101 )
           im = ax.contourf( x1, x2, phi[t], clevels, cmap=cmap )
           ax.set_title( r'Potential $\phi$ at $t = %g$' %(t*dt) )
        if ( arg == 'Ex' ):
           clevels = np.linspace( min_Ex, max_Ex, 101 )
           im = ax.contourf( x1, x2, Ex[t], clevels, cmap=cmap )
           ax.set_title( r'Electric field component $E_x$ at $t = %g$' %(t*dt) )
        if ( arg == 'Ey' ):
           clevels = np.linspace( min_Ey, max_Ey, 101 )
           im = ax.contourf( x1, x2, Ey[t], clevels, cmap=cmap )
           ax.set_title( r'Electric field component $E_y$ at $t = %g$' %(t*dt) )
        if ( arg == 'Em' ):
           clevels = np.linspace( min_Em, max_Em, 101 )
           im = ax.contourf( x1, x2, Em[t], clevels, cmap=cmap )
           ax.set_title( r'Electric field magnitude $|E|$ at $t = %g$' %(t*dt) )
        ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
        ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=0.5 )
        ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
        ax.set_xlabel( r'$x$' )
        ax.set_ylabel( r'$y$', rotation=0 )
        ax.set_aspect( 'equal' )
        fg.colorbar( im, cax=cax )
        fg.canvas.draw()
        plt.pause(0.1)

#-------------------------------------------------------------------------------
# Scalar diagnostics
#-------------------------------------------------------------------------------

# Analytical growth rate
l = 9
omega = 0.179630959411
coeff = 3.5e-08

# Load data from file
# y1 = mass, y2 = energy, y3 = L2-norm of phi
file_name = fold_name+'scalar_diagnostics.dat'
xt = np.loadtxt( file_name, usecols=[0] )
y1 = np.loadtxt( file_name, usecols=[1] )
y2 = np.loadtxt( file_name, usecols=[2] )
y3 = np.loadtxt( file_name, usecols=[3] )

ya = coeff * np.exp( omega * xt )

## Limit t in some range
#mask = ( xt >= 0. )
#xt = xt[mask]
#y1 = y1[mask]
#ya = ya[mask]

# Plot time evolution of L2-norm of phi
def plot_scalar_data( arg ):

    fg = plt.figure( figsize=[9.5,6.0] )
    ax = fg.add_subplot(111)
    fs = 16
    if ( arg == 'mass' ):
      ax.plot( xt, y1, '-r' , lw=1.5, label='Mass conservation' )
      ax.set_ylabel( r'$\int dx\,dy\,\rho$', fontsize=fs )
    if ( arg == 'energy' ):
      ax.plot( xt, y2, '-r' , lw=1.5, label='Energy conservation' )
      ax.set_ylabel( r'$\int dx\,dy\,|\nabla\phi|^2$', fontsize=fs )
    if ( arg == 'l2_norm' ):
      ax.plot( xt, y3, '-r' , lw=1.5, label=r'Numerical instability (mode $m=%d$)' %(l) )
      ax.plot( xt, ya, '--k', lw=1.5, label=r'Analytical growth rate: $\gamma\approx %g$' %(omega) )
      ax.set_ylabel( r'$||\phi-\phi_0||_{L^2}$', fontsize=fs )
      ax.set_yscale( 'log' )
      ymin = ax.get_ylim()[0]
      ymax = 2.*max(y3)
      ax.set_ylim( [ymin,ymax] )
    ax.set_xlabel( r'$t$', fontsize=fs )
    ax.legend( loc='best', fontsize=fs )
    ax.tick_params( labelsize=fs )
    ax.grid()
    fg.tight_layout()
    fg.show()

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

def plot_advection_field( arg ):

    fg = plt.figure(figsize=[9.0,9.0])
    ax = fg.add_subplot(111)
    cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )

    for t in tt:
        ax.clear()
        cax.clear()
        if ( arg == 'expansion' ):
           clevels = np.linspace( min_expansion, max_expansion, 101 )
           im = ax.contourf( x_meshgrid, y_meshgrid, expansion[t], clevels, cmap=cmap )
           ax.set_title( r'Expansion rate of advection field at $t = %g$' %(t*dt) )
        if ( arg == 'vorticity' ):
           clevels = np.linspace( min_vorticity, max_vorticity, 101 )
           im = ax.contourf( x_meshgrid, y_meshgrid, vorticity[t], clevels, cmap=cmap )
           ax.set_title( r'Vorticity of advection field at $t = %g$' %(t*dt) )
        if ( arg == 'shearrate' ):
           clevels = np.linspace( min_shearrate, max_shearrate, 101 )
           im = ax.contourf( x_meshgrid, y_meshgrid, shearrate[t], clevels, cmap=cmap )
           ax.set_title( r'Shear rate of advection field at $t = %g$' %(t*dt) )
        nx = 16
        ax.plot( x_meshgrid[:,::nx], y_meshgrid[:,::nx], color='lightgrey', lw=0.5 )
        ax.plot( x_meshgrid.transpose()[:,::nx], y_meshgrid.transpose()[:,::nx], color='lightgrey', lw=0.5 )
        ax.plot( x1[:,::nr], x2[:,::nr], color='k', lw=0.5 )
        ax.plot( x1[:,n1-1], x2[:,n1-1], color='k', lw=0.5 )
        ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='k', lw=0.5 )
        ax.set_xlim( [ -1., 1. ] )
        ax.set_ylim( [ -1., 1. ] )
        ax.set_xlabel( r'$x$' )
        ax.set_ylabel( r'$y$', rotation=0 )
        ax.set_aspect( 'equal' )
        fg.colorbar( im, cax=cax )
        fg.canvas.draw()
        plt.pause(0.1)

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

fg = plt.figure(figsize=[9.0,9.0])
ax = fg.add_subplot(111)
ax.streamplot( x1_mesh, x2_mesh, -Ey_cart[t], Ex_cart[t], density=6. )
if ( nc > 0 ):
   for ic in range(nc):
       xc = point_charges[t][ic,0]
       yc = point_charges[t][ic,1]
       ax.plot( xc, yc, marker='o', color='k', markersize=5. )
ax.set_title( r'Stream lines of advection field at $t = %g$ (inertial frame)' %(t*dt) )
ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=0.5 )
ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
ax.set_xlabel( r'$x$' )
ax.set_ylabel( r'$y$', rotation=0 )
ax.set_aspect( 'equal' )
fg.show()

fg = plt.figure(figsize=[9.0,9.0])
ax = fg.add_subplot(111)
ax.streamplot( x1_mesh, x2_mesh, Ax_rot, Ay_rot, density=6. )
if ( nc > 0 ):
   for ic in range(nc):
       xc = point_charges[t][ic,0]
       yc = point_charges[t][ic,1]
       ax.plot( xc, yc, marker='o', color='k', markersize=5. )
ax.set_title( r'Stream lines of advection field at $t = %g$ (rotating frame)' %(t*dt) )
ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=0.5 )
ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
ax.set_xlabel( r'$x$' )
ax.set_ylabel( r'$y$', rotation=0 )
ax.set_aspect( 'equal' )
fg.show()

def return_index( minval, maxval, array ):

    imin = np.argmin( abs( array - minval ) )
    imax = np.argmin( abs( array - maxval ) )

    index = slice( imin, imax+1 )

    return index

index_x1 = return_index(  0.3, 0.5, x1_cart )
index_x2 = return_index( -0.1, 0.1, x2_cart )

index = [ index_x2, index_x1 ]

fg = plt.figure(figsize=[9.0,9.0])
ax = fg.add_subplot(111)
ax.streamplot( x1_mesh[index], x2_mesh[index], Ax_rot[index], Ay_rot[index], density=6. )
if ( nc > 0 ):
   for ic in range(nc):
       xc = point_charges[t][ic,0]
       yc = point_charges[t][ic,1]
       ax.plot( xc, yc, marker='o', color='k', markersize=5. )
ax.grid()
ax.set_title( r'Stream lines of advection field at $t = %g$ (rotating frame)' %(t*dt) )
ax.set_xlabel( r'$x$' )
ax.set_ylabel( r'$y$', rotation=0 )
fg.show()
