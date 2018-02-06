import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

# Lists of keys
k1_list = ['circle','target','czarny']
k2_list = ['analytic','discrete','intermed','interpol']

# Dictionary of HDF5 files
h5 = dict()
for k1 in k1_list:
    for k2 in k2_list:
        file_name = 'mapping_' + k2 + '_' + k1 + '.h5'
        h5[(k1,k2)] = h5py.File( file_name, mode='r' )

# Load x1 mesh from HDF5 files
x1 = dict()
for k1 in k1_list:
    # cycle over all keys except 'interpol'
    for k2 in k2_list[:-1]:
        x1[(k1,k2)] = h5[(k1,k2)]['x1'].value

# Load x2 mesh from HDF5 files
x2 = dict()
for k1 in k1_list:
    # cycle over all keys except 'interpol'
    for k2 in k2_list[:-1]:
        x2[(k1,k2)] = h5[(k1,k2)]['x2'].value

# Load Jacobian from HDF5 files
JJ = dict()
for k1 in k1_list:
    # cycle over all keys except 'interpol'
    for k2 in k2_list[:-1]:
        JJ[(k1,k2)] = h5[(k1,k2)]['jacobian'].value

# Load control points along x1 from HDF5 files
c_x1 = dict()
for k1 in k1_list:
    c_x1[k1] = h5[(k1,'discrete')]['c_x1'].value

# Load control points along x2 from HDF5 files
c_x2 = dict()
for k1 in k1_list:
    c_x2[k1] = h5[(k1,'discrete')]['c_x2'].value

# Load interpolation error from HDF5 files
interp_funct = dict()
interp_error = dict()
for k1 in k1_list:
    interp_funct[k1] = h5[(k1,'interpol')]['interp_funct'].value
    interp_error[k1] = h5[(k1,'interpol')]['interp_error'].value

# Close HDF5 files
for k1 in k1_list:
    for k2 in k2_list:
        h5[(k1,k2)].close()

#-------------------------------------------------------------------------------
# OUTPUT
#-------------------------------------------------------------------------------

for k1 in k1_list:
    print()
    print( ' %s mapping' %k1 )
    print( ' ==============' )
    print()
    print(' Maximum absolute errors between analytical and discrete mapping:')
    print()
    print( ' x:', np.amax( np.abs( x1[(k1,'analytic')] - x1[(k1,'discrete')] ) ) )
    print( ' y:', np.amax( np.abs( x2[(k1,'analytic')] - x2[(k1,'discrete')] ) ) )
    print( ' J:', np.amax( np.abs( JJ[(k1,'analytic')] - JJ[(k1,'discrete')] ) ) )

#-------------------------------------------------------------------------------
# PLOTS
#-------------------------------------------------------------------------------

# Mesh from analytical, discrete and intermediate mappings plus control points
fg = plt.figure()
i_plot = 1
for k1 in k1_list:
    ax = fg.add_subplot(1,3,i_plot)
    for k2 in k2_list:
        # analytical mapping
        if (k2 == 'analytic'):
           ax.plot( x1[(k1,k2)], x2[(k1,k2)], color='b' )
           ax.plot( x1[(k1,k2)].transpose(), x2[(k1,k2)].transpose(), color='b' )
        # discrete mapping
        elif (k2 == 'discrete'):
           ax.plot( x1[(k1,k2)], x2[(k1,k2)], '.', color='r' )
           ax.plot( x1[(k1,k2)].transpose(), x2[(k1,k2)].transpose(), '.', color='r' )
        # intermediate mapping
        elif (k2 == 'intermed'):
           ax.plot( x1[(k1,k2)], x2[(k1,k2)], color='k', lw=0.5)
           ax.plot( x1[(k1,k2)].transpose(), x2[(k1,k2)].transpose(), color='k', lw=0.5 )
    # control points
    #ax.plot( c_x1[k1].ravel(), c_x2[k1].ravel(), 'x', color='k' )
    # plot style
    ax.set_aspect( 'equal' )
    ax.set_title( '%s mapping: mesh and control points' %k1 )
    i_plot = i_plot + 1
fg.show()

# Contour plot of interpolation function and error
for k1 in k1_list:
    fg = plt.figure()
    # function profile
    ax = fg.add_subplot(1,2,1)
    cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )
    clevels = np.linspace( interp_funct[k1].min(), interp_funct[k1].max(), 50 )
    im = ax.contourf( x1[(k1,'discrete')], x2[(k1,'discrete')], interp_funct[k1], clevels )
    ax.set_aspect( 'equal' )
    ax.set_title( '%s mapping: function profile' %k1 )
    fg.colorbar( im, cax=cax )
    # interpolation error
    ax = fg.add_subplot(1,2,2)
    cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )
    clevels = np.linspace( interp_error[k1].min(), interp_error[k1].max(), 50 )
    im = ax.contourf( x1[(k1,'discrete')], x2[(k1,'discrete')], interp_error[k1], clevels )
    ax.set_aspect( 'equal' )
    ax.set_title( '%s mapping: interpolation error' %k1 )
    fg.colorbar( im, cax=cax )
    fg.show()

#-------------------------------------------------------------------------------
# Advection tests
#-------------------------------------------------------------------------------

h5 = h5py.File( 'mapping_test_advection.h5', mode='r' )

# meshes
x1 = h5['x1'].value
x2 = h5['x2'].value

# number of iterations (including initial time)
Ni = h5.attrs['iterations']+1

f = dict()
f_ex = dict()
for i in range(Ni):
    s = str(i)
    f[s] = h5['f_'+s].value
    f_ex[s] = h5['f_ex_'+s].value

h5.close()

# find minimum of f
min_f = 1.0
for i in range(Ni):
    if ( f[str(i)].min() < min_f ):
       min_f = f[str(i)].min()

# find maximum of f
max_f = 0.0
for i in range(Ni):
    if ( f[str(i)].max() > max_f ):
       max_f = f[str(i)].max()

# find minimum of f exact
min_f_ex = 1.0
for i in range(Ni):
    if ( f_ex[str(i)].min() < min_f_ex ):
       min_f_ex = f_ex[str(i)].min()

# find maximum of f exact
max_f_ex = 0.0
for i in range(Ni):
    if ( f_ex[str(i)].max() > max_f_ex ):
       max_f_ex = f_ex[str(i)].max()

# compute error
e = dict()
for i in range(Ni):
    s = str(i)
    e[s] = f[s] - f_ex[s]

# find minimum of error
min_err = 1.0
for i in range(Ni):
    if ( e[str(i)].min() < min_err ):
       min_err = e[str(i)].min()

# find maximum of error
max_err = 0.0
for i in range(Ni):
    if ( e[str(i)].max() > max_err ):
       max_err = e[str(i)].max()

# Animated plot
def plot_time_evolution( arg ):

    fg = plt.figure(figsize=[9.0,9.0])
    ax = fg.add_subplot(1,1,1)
    cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )

    for i in range(Ni):
        ax.clear()
        cax.clear()
        if ( arg == 'f' ):
           clevels = np.linspace( min_f, max_f, 50 )
           im = ax.contourf( x1, x2, f[str(i)], clevels )
           ax.set_title( 'Distribution function at t = %i' %i )
        elif ( arg == 'f_ex' ):
           clevels = np.linspace( min_f_ex, max_f_ex, 50 )
           im = ax.contourf( x1, x2, f_ex[str(i)], clevels )
           ax.set_title( 'Exact solution at t = %i' %i )
        elif ( arg == 'error' ):
           clevels = np.linspace( min_err, max_err, 50 )
           im = ax.contourf( x1, x2, e[str(i)], clevels )
           ax.set_title( 'Error on solution at t = %i' %i )
        ax.set_aspect( 'equal' )
        fg.colorbar( im, cax=cax )
        fg.canvas.draw()
        plt.pause(1.0e-05)

print()
print( ' Testing 2D advection on Czarny mapping' )
print( ' ======================================' )
print()

plot_time_evolution( 'f' )
plot_time_evolution( 'error' )
