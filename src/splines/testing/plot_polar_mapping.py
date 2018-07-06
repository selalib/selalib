import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib import cm, rc

# Lists of keys
k1_list = ['circle','target','czarny']
k2_list = ['analytic','discrete','intermed','interpol']

# Dictionary of HDF5 files
h5 = dict()
for k1 in k1_list:
    for k2 in k2_list:
        file_name = 'mapping_' + k2 + '_' + k1 + '.h5'
        h5[(k1,k2)] = h5py.File( file_name, mode='r' )

#h5c = dict()
#h5c['discrete'] = h5py.File( 'mapping_discrete_connor.h5', mode='r' )
#h5c['intermed'] = h5py.File( 'mapping_intermed_connor.h5', mode='r' )

# Load x1 mesh from HDF5 files
x1 = dict()
for k1 in k1_list:
    # cycle over all keys except 'interpol'
    for k2 in k2_list[:-1]:
        x1[(k1,k2)] = h5[(k1,k2)]['x1'].value

#x1c = dict()
#x1c['discrete'] = h5c['discrete']['x1'].value
#x1c['intermed'] = h5c['intermed']['x1'].value

# Load x2 mesh from HDF5 files
x2 = dict()
for k1 in k1_list:
    # cycle over all keys except 'interpol'
    for k2 in k2_list[:-1]:
        x2[(k1,k2)] = h5[(k1,k2)]['x2'].value

#x2c = dict()
#x2c['discrete'] = h5c['discrete']['x2'].value
#x2c['intermed'] = h5c['intermed']['x2'].value

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

#c_x1c = dict()
#c_x1c = h5c['discrete']['c_x1'].value

# Load control points along x2 from HDF5 files
c_x2 = dict()
for k1 in k1_list:
    c_x2[k1] = h5[(k1,'discrete')]['c_x2'].value

#c_x2c = dict()
#c_x2c = h5c['discrete']['c_x2'].value

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

#h5c['discrete'].close()
#h5c['intermed'].close()

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

# text style
#plt.rc( 'text', usetex=True )
#plt.rc( 'font', family='serif' )

## Mesh from analytical, discrete and intermediate mappings plus control points
#fg = plt.figure()
#i_plot = 1
#for k1 in k1_list:
#    ax = fg.add_subplot(1,3,i_plot)
#    for k2 in k2_list:
#        # analytical mapping
#        if (k2 == 'analytic'):
#           ax.plot( x1[(k1,k2)], x2[(k1,k2)], color='b' )
#           ax.plot( x1[(k1,k2)].transpose(), x2[(k1,k2)].transpose(), color='b' )
#        # discrete mapping
#        elif (k2 == 'discrete'):
#           ax.plot( x1[(k1,k2)], x2[(k1,k2)], '.', color='r' )
#           ax.plot( x1[(k1,k2)].transpose(), x2[(k1,k2)].transpose(), '.', color='r' )
#        # intermediate mapping
#        elif (k2 == 'intermed'):
#           ax.plot( x1[(k1,k2)], x2[(k1,k2)], color='k', lw=0.5)
#           ax.plot( x1[(k1,k2)].transpose(), x2[(k1,k2)].transpose(), color='k', lw=0.5 )
#    # control points
#    #ax.plot( c_x1[k1].ravel(), c_x2[k1].ravel(), 'x', color='k' )
#    # plot style
#    ax.set_aspect( 'equal' )
#    ax.set_title( '%s mapping: mesh and control points' %k1 )
#    i_plot = i_plot + 1
#fg.show()

## Temporary plot (Connor)
#fs = 10
#fg = plt.figure()
#ax = fg.add_subplot(1,1,1)
## intermediate mapping
#k2 = 'intermed'
#ax.plot( x1c[k2], x2c[k2], color='lightgrey', lw=0.5)
#ax.plot( x1c[k2].transpose(), x2c[k2].transpose(), color='lightgrey', lw=0.5 )
## discrete mapping
#k2 = 'discrete'
##nr = 16
##ax.plot( x1c[k2][:,::nr], x2c[k2][:,::nr], color='b', lw=0.5 )
##ax.plot( x1c[k2][:,127], x2c[k2][:,127], color='b', lw=0.5 )
##ax.plot( x1c[k2].transpose()[:,::nr], x2c[k2].transpose()[:,::nr], color='b', lw=0.5 )
#ax.plot( x1c[k2], x2c[k2], color='b', lw=0.5 )
#ax.plot( x1c[k2].transpose(), x2c[k2].transpose(), color='b', lw=0.5 )
## plot style
#ax.set_xlabel( r'$x$', fontsize=fs )
#ax.set_ylabel( r'$y$', fontsize=fs, rotation=0 )
#ax.set_aspect( 'equal' )
#fg.tight_layout()
#fg.show()
#fg.savefig( '/home/ezoni/Desktop/mapping.pdf', dpi=300 )

## Contour plot of interpolation function and error
#for k1 in k1_list:
#    fg = plt.figure()
#    # function profile
#    ax = fg.add_subplot(1,2,1)
#    cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )
#    clevels = np.linspace( interp_funct[k1].min(), interp_funct[k1].max(), 50 )
#    im = ax.contourf( x1[(k1,'discrete')], x2[(k1,'discrete')], interp_funct[k1], clevels )
#    ax.set_aspect( 'equal' )
#    ax.set_title( '%s mapping: function profile' %k1 )
#    fg.colorbar( im, cax=cax )
#    # interpolation error
#    ax = fg.add_subplot(1,2,2)
#    cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )
#    clevels = np.linspace( interp_error[k1].min(), interp_error[k1].max(), 50 )
#    im = ax.contourf( x1[(k1,'discrete')], x2[(k1,'discrete')], interp_error[k1], clevels )
#    ax.set_aspect( 'equal' )
#    ax.set_title( '%s mapping: interpolation error' %k1 )
#    fg.colorbar( im, cax=cax )
#    fg.show()

## save figures
#fs = 10
##for k1 in k1_list:
#for k1 in ['czarny']:
#    fg = plt.figure()
#    ax = fg.add_subplot(1,1,1)
#    # intermediate mapping
#    k2 = 'intermed'
#    ax.plot( x1[(k1,k2)], x2[(k1,k2)], color='k', lw=0.3)
#    ax.plot( x1[(k1,k2)].transpose(), x2[(k1,k2)].transpose(), color='k', lw=0.3 )
#    # discrete mapping
#    k2 = 'discrete'
#    nr = 16
#    ax.plot( x1[(k1,k2)][:,::nr], x2[(k1,k2)][:,::nr], color='b', lw=0.5 )
#    ax.plot( x1[(k1,k2)][:,127], x2[(k1,k2)][:,127], color='b', lw=0.5 )
#    ax.plot( x1[(k1,k2)].transpose()[:,::nr], x2[(k1,k2)].transpose()[:,::nr], color='b', lw=0.5 )
#    # plot style
#    ax.set_xlim( -1., 1.3 )
#    ax.annotate( r'$X$', xy=(0.1,-1.55) )
#    ax.annotate( r'$Y$', xy=(1.1,0.05) )
#    ax.set_xlabel( r'$x$', fontsize=fs )
#    ax.set_ylabel( r'$y$', fontsize=fs, rotation=0 )
##    ax.set_title( 'Physical domain', fontsize=fs )
#    ax.set_aspect( 'equal' )
#    fg.tight_layout()
#    fg.show()
#    fg.savefig( './'+k1+'_pseudo.pdf', dpi=300 )

#-------------------------------------------------------------------------------
# ADVECTION TESTS
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
    ax = fg.add_subplot(111)
    cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )

    for i in range(Ni):
        ax.clear()
        cax.clear()
        if ( arg == 'f' ):
           clevels = np.linspace( min_f, max_f, 100 )
           # contour plot
           im = ax.contourf( x1, x2, f[str(i)], clevels, cmap='jet' )
           ax.set_title( 'Distribution function at t = %i' %i )
        elif ( arg == 'f_ex' ):
           clevels = np.linspace( min_f_ex, max_f_ex, 100 )
           im = ax.contourf( x1, x2, f_ex[str(i)], clevels, cmap='jet' )
           ax.set_title( 'Exact solution at t = %i' %i )
        elif ( arg == 'error' ):
           clevels = np.linspace( min_err, max_err, 100 )
           im = ax.contourf( x1, x2, e[str(i)], clevels, cmap='jet' )
           ax.set_title( 'Error on solution at t = %i' %i )
        # plot grids
        nr = 16
        ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
        ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
        # style
        ax.set_aspect( 'equal' )
        fg.colorbar( im, cax=cax )
        fg.canvas.draw()
        plt.pause(1.0e-03)

print()
print( ' Testing 2D advection on Czarny mapping' )
print( ' ======================================' )
print()

## save figures
#xc = 0.25
#yc = 0.0
#
#s = 0.4
#t = np.linspace(0.0,2.0*np.pi,100)
#
#xm = s*np.cos(t)+xc
#ym = s*np.sin(t)+yc
#
#fs = 10
#fg = plt.figure()
#ax = fg.add_subplot(111)
#clevels = np.linspace( min_f, max_f, 100 )
#i  = 17
## plot grids
#nr = 16
#ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
#ax.plot( x1[:,127 ], x2[:,127 ], color='lightgrey', lw=0.5 )
#ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
## plot center of rotation
#ax.plot( xc, yc, '.', color='w' )
#ax.annotate( r'$(x_c,y_c)$', xy=[xc,yc-0.15], color='w',
#             horizontalalignment='center', verticalalignment='bottom', fontsize=fs )
## plot analytical trajectory
#ax.plot( xm, ym, '--', color='w' )
## contour plot
#im = ax.contourf( x1, x2, f[str(i)], clevels, cmap='jet' )
## style
#ax.set_xlabel( r'$x$', fontsize=fs )
#ax.set_ylabel( r'$y$', fontsize=fs, rotation=0 )
#ax.set_title( r'$f$ $(x,y)$' )
#ax.set_aspect( 'equal' )
#fg.colorbar( im )
#fg.tight_layout()
#fg.show()
#fg.savefig( './advection.pdf', dpi=300 )
