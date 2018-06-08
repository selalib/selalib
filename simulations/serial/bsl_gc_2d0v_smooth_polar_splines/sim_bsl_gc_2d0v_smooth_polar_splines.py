import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

# Open HDF5 file
file_name = 'sim_bsl_gc_2d0v_smooth_polar_splines.h5'
h5 = h5py.File( file_name, mode='r' )

# Load x1 and x2 mesh
x1 = h5['x1'].value
x2 = h5['x2'].value

# Load Jacobian
JJ = h5['jacobian'].value

# Load control points along x1 and x2
c_x1 = h5['c_x1'].value
c_x2 = h5['c_x2'].value

# Load attributes
dt = h5.attrs['time_step']
Ni = h5.attrs['iterations']

# Load density rho
rho = dict()
for i in range(Ni):
    s = str(i)
    rho[s] = h5['rho_'+s].value

# Close HDF5 file
h5.close()

#-------------------------------------------------------------------------------
# PLOTS
#-------------------------------------------------------------------------------

# find minimum of rho
min_rho = 1.0
for i in range(Ni):
    if ( rho[str(i)].min() < min_rho ):
       min_rho = rho[str(i)].min()

# find maximum of rho
max_rho = 0.0
for i in range(Ni):
    if ( rho[str(i)].max() > max_rho ):
       max_rho = rho[str(i)].max()

# Animated plot
def plot_time_evolution( arg ):

    fg = plt.figure(figsize=[9.0,9.0])
    ax = fg.add_subplot(111)
    cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )

    for i in range(Ni):
        ax.clear()
        cax.clear()
        clevels = np.linspace( min_rho, max_rho, 100 )
        # contour plot
        im = ax.contourf( x1, x2, rho[str(i)], clevels, cmap='jet' )
        # grid
        nr = 16
        ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
        ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
        # style
        ax.set_xlabel( r'$x$' )#, fontsize=fs )
        ax.set_ylabel( r'$y$', rotation=0 )#, fontsize=fs )
        ax.set_title( r'Density $\rho$ at $t = %g$' %(i*dt) )
        ax.set_aspect( 'equal' )
        fg.colorbar( im, cax=cax )
        fg.canvas.draw()
        plt.pause(0.1)

plot_time_evolution( 'rho' )
