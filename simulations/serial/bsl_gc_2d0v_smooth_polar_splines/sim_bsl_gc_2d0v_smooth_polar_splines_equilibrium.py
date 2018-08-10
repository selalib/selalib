import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

# Read folder name from standard input
fold_name = sys.argv[1]

# Open HDF5 file
file_name = fold_name+'sim_bsl_gc_2d0v_smooth_polar_splines_equilibrium.h5'
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
Ni = h5.attrs['iterations']+1

tt = np.ndarray( Ni, dtype=int )
for i in range(0,Ni):
    tt[i] = i

# Load data
rho = dict()
phi = dict()
for i in tt:
    t = str(i)
    rho[t] = h5['rho_'+t].value
    phi[t] = h5['phi_'+t].value

# Close HDF5 file
h5.close()

#-------------------------------------------------------------------------------
# Plot time evolution of RHO
#-------------------------------------------------------------------------------

# find minimum and maximum
def minmax( Q ):
    min_Q = 1.0
    max_Q = 0.0
    for i in tt:
        t = str(i)
        if ( Q[t].min() < min_Q ):
           min_Q = Q[t].min()
        if ( Q[t].max() > max_Q ):
           max_Q = Q[t].max()
    return [ min_Q, max_Q ]

min_rho, max_rho = minmax( rho )
min_phi, max_phi = minmax( phi )

# Animated plot
def plot_iterations( arg ):

    fg = plt.figure(figsize=[9.0,9.0])
    ax = fg.add_subplot(111)
    cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )

    for i in tt:
        t = str(i)
        ax.clear()
        cax.clear()
        if ( arg == 'rho' ):
           clevels = np.linspace( min_rho, max_rho, 100 )
           im = ax.contourf( x1, x2, rho[t], clevels, cmap='jet' )
           ax.set_title( r'Density $\rho$ at iteration $%g$' %i )
        if ( arg == 'phi' ):
           clevels = np.linspace( min_phi, max_phi, 100 )
           im = ax.contourf( x1, x2, phi[t], clevels, cmap='jet' )
           ax.set_title( r'Potential $\phi$ at iteration $%g$' %i )
        nr = 16
        ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
        ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=0.5 )
        ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
        ax.set_xlabel( r'$x$' )#, fontsize=fs )
        ax.set_ylabel( r'$y$', rotation=0 )#, fontsize=fs )
        ax.set_aspect( 'equal' )
        fg.colorbar( im, cax=cax )
        fg.canvas.draw()
        plt.pause(1.0)
