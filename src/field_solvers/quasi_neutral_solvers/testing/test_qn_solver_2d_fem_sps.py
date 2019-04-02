import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

#-------------------------------------------------------------------------------
# Load data generated from Fortran
#-------------------------------------------------------------------------------

# Open HDF5 file
file_name = 'qn_solver_2d_fem_sps.h5'
h5 = h5py.File( file_name, mode='r' )

# Load parameters
n1 = h5.attrs['n1']
n2 = h5.attrs['n2']
p1 = h5.attrs['p1']
p2 = h5.attrs['p2']

# Load data
x1 = h5['x1'].value
x2 = h5['x2'].value
#A   = h5['A' ].value
#M   = h5['M' ].value
#Ap  = h5['Ap'].value
#Mp  = h5['Mp'].value
##L   = h5['L' ].value
##b   = h5['b' ].value
x   = h5['x' ].value
phi = h5['phi'].value
phi_spl = h5['phi_spl'].value
err = h5['error'].value

# Close HDF5 file
h5.close()

#-------------------------------------------------------------------------------
# TEST 1: domain area, full mass matrix
#-------------------------------------------------------------------------------

#ones = np.ones( M.shape[0] )
#area = M.dot( ones ).sum()
##error = np.pi - area
#
#print()
#print( ' **************************************' )
#print( ' TEST #1: domain area, full mass matrix' )
#print( ' **************************************' )
#print()
#print( ' area  = {}'.format( area  ) )
##print( ' error = {}'.format( error ) )

#-------------------------------------------------------------------------------
# TEST 1: domain area, projected mass matrix
#-------------------------------------------------------------------------------

#ones = np.ones( Mp.shape[0] )
#area = Mp.dot( ones ).sum()
##error = np.pi - area
#
#print()
#print( ' *******************************************' )
#print( ' TEST #1: domain area, projected mass matrix' )
#print( ' *******************************************' )
#print()
#print( ' area  = {}'.format( area  ) )
##print( ' error = {}'.format( error ) )

#-------------------------------------------------------------------------------
# TEST 2: kernel of full stiffness matrix
#-------------------------------------------------------------------------------

#ones = np.ones( A.shape[0] )
#r = A.dot( ones )
#t = M.dot( r )
#
#err_l2_norm = np.sqrt( np.dot( r, t ) )
#err_maxnorm = abs( r ).max()
#
#print()
#print( ' ****************************************' )
#print( ' TEST #2: kernel of full stiffness matrix' )
#print( ' ****************************************' )
#print()
#print( ' L^2-norm of result = {}'.format( err_l2_norm ) )
#print( ' max-norm of result = {}'.format( err_maxnorm ) )

#-------------------------------------------------------------------------------
# TEST 2: kernel of projected stiffness matrix
#-------------------------------------------------------------------------------

#ones = np.ones( Ap.shape[0] )
#r = Ap.dot( ones )
#t = Mp.dot( r )
#
#err_l2_norm = np.sqrt( np.dot( r, t ) )
#err_maxnorm = abs( r ).max()
#
#print()
#print( ' *********************************************' )
#print( ' TEST #2: kernel of projected stiffness matrix' )
#print( ' *********************************************' )
#print()
#print( ' L^2-norm of result = {}'.format( err_l2_norm ) )
#print( ' max-norm of result = {}'.format( err_maxnorm ) )

#-------------------------------------------------------------------------------
# PLOTS
#-------------------------------------------------------------------------------

#def run_plot_M():
#    fig,ax = plt.subplots(1,1)
#    ax.set_title( "Tensor-product mass matrix M" )
#    im = ax.matshow( M, norm=colors.LogNorm(), cmap='hot_r' )
#    cb = fig.colorbar( im, ax=ax )
#    fig.show()
#    return locals()
#
#def run_plot_Mp():
#    fig,ax = plt.subplots(1,1)
#    ax.set_title( "C^1 mass matrix M' (projection of M)" )
#    im = ax.matshow( Mp, norm=colors.LogNorm(), cmap='hot_r' )
#    cb = fig.colorbar( im, ax=ax )
#    fig.show()
#    return locals()
#
#def run_plot_A():
#    fig,ax = plt.subplots(1,1)
#    ax.set_title( "Tensor-product stiffness matrix A" )
#    norm = colors.SymLogNorm( linthresh=1.0e-8, vmin=-A.max(), vmax=A.max() )
#    im = ax.matshow( A, norm=norm, cmap='seismic' )
#    cb = fig.colorbar( im, ax=ax )
#    fig.show()
#    return locals()
#
#def run_plot_Ap():
#    fig,ax = plt.subplots(1,1)
#    ax.set_title( "C^1 stiffness matrix A' (projection of A)" )
#    norm_Ap = colors.SymLogNorm( linthresh=1.0e-6, vmin=-Ap.max(), vmax=Ap.max() )
#    im = ax.matshow( Ap, norm=norm_Ap, cmap='seismic' )
#    cb = fig.colorbar( im, ax=ax )
#    fig.show()
#    return locals()
#
#def run_plot_Mp_diff():
#    fig,ax = plt.subplots(1,1)
#    ax.set_title( "C^1 stiffness matrix M' (difference)" )
#    norm = colors.SymLogNorm( linthresh=1.0e-6, vmin=-Mp.max(), vmax=Mp.max() )
#    im = ax.matshow( Mp-Mp_py, norm=norm, cmap='seismic' )
#    cb = fig.colorbar( im, ax=ax )
#    fig.show()
#    return locals()
#
#def run_plot_Ap_diff():
#    fig,ax = plt.subplots(1,1)
#    ax.set_title( "C^1 stiffness matrix A' (difference)" )
#    norm = colors.SymLogNorm( linthresh=1.0e-6, vmin=-Ap.max(), vmax=Ap.max() )
#    im = ax.matshow( Ap-Ap_py, norm=norm, cmap='seismic' )
#    cb = fig.colorbar( im, ax=ax )
#    fig.show()
#    return locals()
#
##def run_plot_L_diff():
##    fig,ax = plt.subplots(1,1)
##    ax.set_title( "Change of basis matrix L (difference)" )
##    norm = colors.SymLogNorm( linthresh=1.0e-6, vmin=-L.max(), vmax=L.max() )
##    im = ax.matshow( L.T-L_py, norm=norm, cmap='seismic' )
##    cb = fig.colorbar( im, ax=ax )
##    fig.show()
##    return locals()

#-------------------------------------------------------------------------------

#plot_M  = run_plot_M ()
#plot_A  = run_plot_A ()
#plot_Mp = run_plot_Mp()
#plot_Ap = run_plot_Ap()
#plot_Mp_diff = run_plot_Mp_diff()
#plot_Ap_diff = run_plot_Ap_diff()
##plot_L_diff  = run_plot_L_diff()

phi_plt = np.ndarray([n2+1,n1])
phi_plt[:n2,:] = phi_spl[:,:]
phi_plt[ n2,:] = phi_spl[0,:]

fs = 10
fg = plt.figure()
ax = fg.add_subplot(111)
# contour plot
clevels = np.linspace( phi_plt.min(), phi_plt.max(), 101)
im = ax.contourf( x1, x2, phi_plt, clevels, cmap='jet' )
for c in im.collections:
    c.set_edgecolor('face')
fg.colorbar( im, ax=ax )
# grid
nr = n1//8
ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=0.5 )
ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
# style
ax.set_xlabel( r'$x$', fontsize=fs )
ax.set_ylabel( r'$y$', fontsize=fs, rotation=0 )
ax.set_title( r'Numerical solution: $\phi(x,y)$' )
ax.set_aspect( 'equal' )
fg.tight_layout()
fg.show()
#fg_name = 'qn_solver.pdf'
#fg.savefig( fg_name, dpi=300 )

err_plt = np.ndarray([n2+1,n1])
err_plt[:n2,:] = err[:,:]
err_plt[ n2,:] = err[0,:]

fs = 10
fg = plt.figure()
ax = fg.add_subplot(111)
# contour plot
clevels = np.linspace( err_plt.min(), err_plt.max(), 101)
im = ax.contourf( x1, x2, err_plt, clevels, cmap='seismic' )
for c in im.collections:
    c.set_edgecolor('face')
cb = fg.colorbar( im, ax=ax )
cb.formatter.set_powerlimits((0,0))
cb.formatter.useMathText = True
cb.update_ticks()
# grid
nr = n1//8
ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=0.5 )
ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
# style
ax.set_xlabel( r'$x$', fontsize=fs )
ax.set_ylabel( r'$y$', fontsize=fs, rotation=0 )
ax.set_title( r'Numerical error: $\phi-\phi_{ex}$' )
ax.set_aspect( 'equal' )
fg.tight_layout()
fg.show()
#fg_name = 'qn_solver_error.pdf'
#fg.savefig( fg_name, dpi=300 )
