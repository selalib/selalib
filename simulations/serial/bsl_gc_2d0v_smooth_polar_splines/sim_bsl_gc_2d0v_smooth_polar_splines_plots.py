# run interactively ("-i") after running sim_bsl_gc_2d0v_smooth_polar_splines.py

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

# TeX rendering
import matplotlib
from matplotlib import rc
rc( 'font', **{'family':'sans-serif','sans-serif':['Helvetica']} )
rc( 'text', usetex=True )
matplotlib.rcParams[ 'text.latex.preamble' ] = [ r"\usepackage{amsmath}" ]

# Plot style
fs_small = 10
fs_large = 12
nr = n1//8
lw = 0.5

# color map for plots
cmap = 'seismic'

# angular velocity of background
omega = 0.3332

# time conversion factor (only for vortex in cell simulation)
tconv = 0.417*0.4

# auxiliary functions
#-------------------------------------------------------------------------------
def rotate( x, y, theta ):
    x_new = x*np.cos( theta ) - y*np.sin( theta )
    y_new = x*np.sin( theta ) + y*np.cos( theta )
    return x_new, y_new

def return_index( array, minval, maxval ):

    imin = np.argmin( abs( array - minval ) )
    imax = np.argmin( abs( array - maxval ) )

    index = slice( imin, imax+1 )

    return index
#-------------------------------------------------------------------------------

# point charge trajectory
def plot_point_charge_trajectory( frame, tmax, save=False ):
    fg = plt.figure()
    ax = fg.add_subplot(111)
    if ( frame == 'inertial' ):
        xc = np.fromiter( ( point_charges[t][0,0] for t in range(tmax+1) ), dtype=float )
        yc = np.fromiter( ( point_charges[t][0,1] for t in range(tmax+1) ), dtype=float )
        ax.plot( xc, yc, '-', color='k', markersize=1. )
        ax.set_title( 'Point charge trajectory (inertial frame)', fontsize=fs_large )
    if ( frame == 'rotating' ):
        xc = np.fromiter( ( point_charges[t][0,0] for t in range(tmax+1) ), dtype=float )
        yc = np.fromiter( ( point_charges[t][0,1] for t in range(tmax+1) ), dtype=float )
        th = np.fromiter( ( -omega*(t*dt) for t in range(tmax+1) ), dtype=float )
        xr, yr = rotate( xc, yc, th )
        ax.plot( xr, yr, '-', color='k', markersize=1. )
        ax.set_title( 'Point charge trajectory (rotating frame)', fontsize=fs_large )
    ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=lw )
    ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=lw )
    ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=lw )
    ax.set_xlabel( r'$x$', fontsize=fs_large )
    ax.set_ylabel( r'$y$', fontsize=fs_large, rotation=0 )
    ax.tick_params( axis='both', labelsize=fs_small )
    ax.set_aspect( 'equal' )
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )

# time evolution of radial position (rotating frame)
def plot_point_charge_radial_position( tmax, save=False ):
    fg = plt.figure()
    ax = fg.add_subplot(111)
    xc = np.fromiter( ( point_charges[t][0,0] for t in range(tmax+1) ), dtype=float )
    yc = np.fromiter( ( point_charges[t][0,1] for t in range(tmax+1) ), dtype=float )
    th = np.fromiter( ( -omega*(t*dt) for t in range(tmax+1) ), dtype=float )
    xr, yr = rotate( xc, yc, th )
    rr = np.sqrt( xr**2 + yr**2 )
    tr = np.arctan2( yr, xr ) % ( 2.0 * np.pi )
    ax.plot( th*(-tconv/omega), rr, '-', color='k', markersize=1. )
    ax.grid()
    ax.set_xlabel( r'$t^\prime$', fontsize=fs_large )
    ax.set_ylabel( r'$s$', fontsize=fs_large, rotation=0, labelpad=10 )
    ax.tick_params( axis='both', labelsize=fs_small )
    fg.tight_layout()
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )

# time evolution of angular position (rotating frame)
def plot_point_charge_angular_position( tmax, save=False ):
    fg = plt.figure()
    ax = fg.add_subplot(111)
    xc = np.fromiter( ( point_charges[t][0,0] for t in range(tmax+1) ), dtype=float )
    yc = np.fromiter( ( point_charges[t][0,1] for t in range(tmax+1) ), dtype=float )
    th = np.fromiter( ( -omega*(t*dt) for t in range(tmax+1) ), dtype=float )
    xr, yr = rotate( xc, yc, th )
    rr = np.sqrt( xr**2 + yr**2 )
    tr = np.arctan2( yr, xr ) % ( 2.0 * np.pi )
    ax.plot( th*(-tconv/omega), tr, '-', color='k', markersize=1. )
    ax.grid()
    ax.set_xlabel( r'$T$', fontsize=fs_large )
    ax.set_ylabel( r'$\theta$', fontsize=fs_large, rotation=0 )
    ax.tick_params( axis='both', labelsize=fs_small )
    fg.tight_layout()
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )

# time shot
def plot_time_shot( arg, t=0, save=False ):
    fg = plt.figure()
    ax = fg.add_subplot(111)
    if ( arg == 'rho_eq' ):
        clmin = min_rho_eq
        clmax = max_rho_eq
        clevels = np.linspace( clmin, clmax, 101 )
        im = ax.contourf( x1, x2, rho_eq, clevels, cmap=cmap )
        for c in im.collections:
            c.set_edgecolor('face')
        if ( nc > 0 ):
           for ic in range(nc):
               xc = point_charges[t][ic,0]
               yc = point_charges[t][ic,1]
               ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'$\rho_0(x,y)$', fontsize=fs_large )
    if ( arg == 'rho' ):
        clmin = min_rho
        clmax = max_rho
        clevels = np.linspace( clmin, clmax, 101 )
        im = ax.contourf( x1, x2, rho[t], clevels, cmap=cmap )
        for c in im.collections:
            c.set_edgecolor('face')
        if ( nc > 0 ):
           for ic in range(nc):
               xc = point_charges[t][ic,0]
               yc = point_charges[t][ic,1]
               ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'$\rho(t,x,y)$ at $t = %g$' %(t*dt), fontsize=fs_large )
    if ( arg == 'rho_rot' ):
        clmin = min_rho
        clmax = max_rho
        clevels = np.linspace( clmin, clmax, 101 )
        x1r, x2r = rotate( x1, x2, -omega*(t*dt) )
        im = ax.contourf( x1r, x2r, rho[t], clevels, cmap=cmap )
        for c in im.collections:
            c.set_edgecolor('face')
        if ( nc > 0 ):
           for ic in range(nc):
               xc = point_charges[t][ic,0]
               yc = point_charges[t][ic,1]
               xr, yr = rotate( xc, yc, -omega*(t*dt) )
               ax.plot( xr, yr, marker='o', color='k', markersize=5. )
        ax.set_title( r'$\rho(t^\prime,x^\prime,y^\prime)$ at $t^\prime = %g$' %(t*dt*tconv), fontsize=fs_large )
    if ( arg == 'delta_rho' ):
        clmin = min_delta_rho
        clmax = max_delta_rho
        clevels = np.linspace( clmin, clmax, 101 )
        im = ax.contourf( x1, x2, rho[t]-rho_eq, clevels, cmap=cmap )
        for c in im.collections:
            c.set_edgecolor('face')
        if ( nc > 0 ):
           for ic in range(nc):
               xc = point_charges[t][ic,0]
               yc = point_charges[t][ic,1]
               ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'$\rho_1(t,x,y)$ at $t = %g$' %(t*dt), fontsize=fs_large )
    if ( arg == 'phi_eq' ):
        clmin = min_phi_eq
        clmax = max_phi_eq
        clevels = np.linspace( clmin, clmax, 101 )
        im = ax.contourf( x1, x2, phi_eq, clevels, cmap=cmap )
        for c in im.collections:
            c.set_edgecolor('face')
        if ( nc > 0 ):
           for ic in range(nc):
               xc = point_charges[t][ic,0]
               yc = point_charges[t][ic,1]
               ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'$\Phi_0(x,y)$', fontsize=fs_large )
    if ( arg == 'phi' ):
        clmin = min_phi
        clmax = max_phi
        clevels = np.linspace( clmin, clmax, 101 )
        im = ax.contourf( x1, x2, phi[t], clevels, cmap=cmap )
        for c in im.collections:
            c.set_edgecolor('face')
        if ( nc > 0 ):
           for ic in range(nc):
               xc = point_charges[t][ic,0]
               yc = point_charges[t][ic,1]
               ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'$\Phi(t,x,y)$ at $t = %g$' %(t*dt), fontsize=fs_large )
    if ( arg == 'Em_squared' ):
        clmin = min_Em
        clmax = max_Em
        clevels = np.linspace( clmin, clmax, 101 )
        im = ax.contourf( x1, x2, Em[t], clevels, cmap=cmap )
        for c in im.collections:
            c.set_edgecolor('face')
        if ( nc > 0 ):
           for ic in range(nc):
               xc = point_charges[t][ic,0]
               yc = point_charges[t][ic,1]
               ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'$|\boldsymbol{E}(t,x,y)|^2$ at $t = %g$' %(t*dt), fontsize=fs_large )
    ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=lw )
    ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=lw )
    ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=lw )
    ax.set_xlabel( r'$x$', fontsize=fs_large )
    ax.set_ylabel( r'$y$', fontsize=fs_large, rotation=0 )
    ax.tick_params( axis='both', labelsize=fs_small )
    ax.set_aspect( 'equal' )
    cb = fg.colorbar( im, ax=ax )
    #tk = np.linspace( clmin, clmax, 10 )
    #cb.set_ticks( tk )
    cb.ax.tick_params( labelsize=fs_small )
    cb.formatter.set_powerlimits((0,0))
    cb.update_ticks()
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )

# animated plots
def plot_time_evolution( arg, tmax ):

    fg = plt.figure()
    ax = fg.add_subplot(111)
    cax = make_axes_locatable(ax).append_axes( 'right', size='5%', pad='5%' )

    for t in (ti for ti in tt if ti <= tmax):
        ax.clear()
        cax.clear()
        if ( arg == 'rho' ):
            clmin = min_rho
            clmax = max_rho
            clevels = np.linspace( clmin, clmax, 101 )
            im = ax.contourf( x1, x2, rho[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            if ( nc > 0 ):
                for ic in range(nc):
                    xc = point_charges[t][ic,0]
                    yc = point_charges[t][ic,1]
                    ax.plot( xc, yc, marker='o', color='k', markersize=5. )
            ax.set_title( r'$\rho(t,x,y)$ at $t = %g$' %(t*dt), fontsize=fs_large )
        if ( arg == 'rho_rot' ):
            clmin = min_rho
            clmax = max_rho
            clevels = np.linspace( clmin, clmax, 101 )
            x1r, x2r = rotate( x1, x2, -omega*(t*dt) )
            im = ax.contourf( x1r, x2r, rho[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            if ( nc > 0 ):
                for ic in range(nc):
                    xc = point_charges[t][ic,0]
                    yc = point_charges[t][ic,1]
                    xr, yr = rotate( xc, yc, -omega*(t*dt) )
                    ax.plot( xr, yr, marker='o', color='k', markersize=5. )
            ax.set_xlim([-1.,1.])
            ax.set_ylim([-1.,1.])
            ax.set_title( r'$\rho(t^\prime,x^\prime,y^\prime)$ at $t^\prime = %g$' %(t*dt*tconv), fontsize=fs_large )
        if ( arg == 'delta_rho' ):
            clmin = min_delta_rho
            clmax = max_delta_rho
            clevels = np.linspace( clmin, clmax, 101 )
            im = ax.contourf( x1, x2, rho[t]-rho_eq, clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            ax.set_title( r'$\rho_1(t,x,y)$ at $t = %g$' %(t*dt), fontsize=fs_large )
        if ( arg == 'phi' ):
            clmin = min_phi
            clmax = max_phi
            clevels = np.linspace( clmin, clmax, 101 )
            im = ax.contourf( x1, x2, phi[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            if ( nc > 0 ):
                for ic in range(nc):
                    xc = point_charges[t][ic,0]
                    yc = point_charges[t][ic,1]
                    ax.plot( xc, yc, marker='o', color='k', markersize=5. )
            ax.set_title( r'$\Phi(t,x,y)$ at $t = %g$' %(t*dt), fontsize=fs_large )
        if ( arg == 'Em_squared' ):
            clmin = min_Em
            clmax = max_Em
            clevels = np.linspace( clmin, clmax, 101 )
            im = ax.contourf( x1, x2, Em[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            if ( nc > 0 ):
               for ic in range(nc):
                   xc = point_charges[t][ic,0]
                   yc = point_charges[t][ic,1]
                   ax.plot( xc, yc, marker='o', color='k', markersize=5. )
            ax.set_title( r'$|\boldsymbol{E}(t,x,y)|^2$ at $t = %g$' %(t*dt), fontsize=fs_large )
        ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=lw )
        ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=lw )
        ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=lw )
        ax.set_xlabel( r'$x$', fontsize=fs_large )
        ax.set_ylabel( r'$y$', fontsize=fs_large, rotation=0 )
        ax.tick_params( axis='both', labelsize=fs_small )
        ax.set_aspect( 'equal' )
        cb = fg.colorbar( im, cax=cax )
        cb.formatter.set_powerlimits((0,0))
        cb.update_ticks()
        fg.canvas.draw()
        plt.pause(1e-3)

## Limit t in some range
#mask = ( xt >= 0. )
#xt = xt[mask]
#y1 = y1[mask]
#ya = ya[mask]

# time evolution of scalar diagnostics
def plot_scalar_data( arg, save=False ):

    fg = plt.figure()
    ax = fg.add_subplot(111)
    lw = 1.0
    if ( arg == 'mass' ):
        # optimal size [6.4,2.4]
        ax.plot( xt, (y1[0]-y1[:])/y1[0], '-r' , lw=lw )
        ax.set_ylabel( r'$\delta M$', rotation=0, fontsize=fs_large )
    if ( arg == 'energy' ):
        # optimal size [6.4,2.4]
        ax.plot( xt, (y2[0]-y2[:])/y2[0], '-r' , lw=lw )
        ax.set_ylabel( r'$\delta W$', rotation=0, fontsize=fs_large )
    if ( arg == 'l2_norm_phi' ):
        # optimal size [9.5,6.0]
        ax.plot( xt, y3, '-r' , lw=lw, label=r'Numerical instability (mode $m=%d$)' %(l) )
        ax.plot( xt, ya, '--k', lw=lw, label=r'Analytical growth rate: Im$(\omega\approx %1.2g$' %(gamma) )
        ax.set_ylabel( r'$||\Phi-\Phi_0||_{L^2}$', fontsize=fs_large )
        ax.set_yscale( 'log' )
        ymin = ax.get_ylim()[0]
        ymax = 2.*max(y3)
        ax.set_ylim( [ymin,ymax] )
        ax.legend( loc='best', fontsize=fs_small )
    ax.set_xlabel( r'$t$', fontsize=fs_large )
    ax.tick_params( axis='both', labelsize=fs_small )
    #ax.ticklabel_format( axis='y', style='sci', scilimits=(0,0), useMathText=True )
    ax.grid()
    fg.tight_layout()
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )

# vorticity, shear, expansion
def plot_advection_field( arg, tmax ):

    fg = plt.figure()
    ax = fg.add_subplot(111)
    cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )

    for t in (ti for ti in tt if ti <= tmax):
        ax.clear()
        cax.clear()
        if ( arg == 'vorticity' ):
            clevels = np.linspace( min_vorticity, max_vorticity, 101 )
            im = ax.contourf( x_meshgrid, y_meshgrid, vorticity[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            ax.set_title( r'Vorticity of advection field at $t^\prime = %g$' %(t*dt*tconv), fontsize=fs_large )
        if ( arg == 'shear' ):
            clevels = np.linspace( min_shearrate, max_shearrate, 101 )
            im = ax.contourf( x_meshgrid, y_meshgrid, shearrate[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            ax.set_title( r'Shear rate of advection field at $t^\prime = %g$' %(t*dt*tconv), fontsize=fs_large )
        if ( arg == 'expansion' ):
            clevels = np.linspace( min_expansion, max_expansion, 101 )
            im = ax.contourf( x_meshgrid, y_meshgrid, expansion[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            ax.set_title( r'Expansion rate of advection field at $t^\prime = %g$' %(t*dt*tconv), fontsize=fs_large )
        nx = 16
        ax.plot( x_meshgrid[:,::nx], y_meshgrid[:,::nx], color='lightgrey', lw=lw )
        ax.plot( x_meshgrid.transpose()[:,::nx], y_meshgrid.transpose()[:,::nx], color='lightgrey', lw=lw )
        ax.plot( x1[:,::nr], x2[:,::nr], color='k', lw=lw )
        ax.plot( x1[:,n1-1], x2[:,n1-1], color='k', lw=lw )
        ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='k', lw=lw )
        ax.set_xlim( [ -1., 1. ] )
        ax.set_ylim( [ -1., 1. ] )
        ax.set_xlabel( r'$x$', fontsize=fs_large )
        ax.set_ylabel( r'$y$', fontsize=fs_large, rotation=0 )
        ax.tick_params( axis='both', labelsize=fs_small )
        ax.set_aspect( 'equal' )
        cb = fg.colorbar( im, cax=cax )
        cb.formatter.set_powerlimits((0,0))
        cb.update_ticks()
        fg.canvas.draw()
        plt.pause(0.1)

# stream lines (works for circle only)

index_x1 = return_index( x1_cart,  0.3, 0.5 )
index_x2 = return_index( x2_cart, -0.1, 0.1 )

index = tuple( ( index_x2, index_x1 ) )

def plot_stream_lines_global( frame, save=False ):
    t = 0
    fg = plt.figure()
    ax = fg.add_subplot(111)
    if ( frame == 'inertial' ):
        ax.streamplot( x1_mesh, x2_mesh, -Ey_cart[t], Ex_cart[t], density=6. )
        if ( nc > 0 ):
            for ic in range(nc):
                xc = point_charges[t][ic,0]
                yc = point_charges[t][ic,1]
                ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'Stream lines of advection field at $t = %g$' %(t*dt), fontsize=fs_large )
    if ( frame == 'rotating' ):
        ax.streamplot( x1_mesh, x2_mesh, Ax_rot, Ay_rot, density=6. )
        if ( nc > 0 ):
            for ic in range(nc):
                xc = point_charges[t][ic,0]
                yc = point_charges[t][ic,1]
                ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'Stream lines of advection field at $t^\prime = %g$' %(t*dt*tconv), fontsize=fs_large )
    ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=lw )
    ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=lw )
    ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=lw )
    ax.set_xlabel( r'$x$', fontsize=fs_large )
    ax.set_ylabel( r'$y$', fontsize=fs_large, rotation=0 )
    ax.tick_params( axis='both', labelsize=fs_small )
    ax.set_aspect( 'equal' )
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )

def plot_stream_lines_local( frame, save=False ):
    t = 0
    d = 1.5
    fg = plt.figure()
    ax = fg.add_subplot(111)
    if ( frame == 'inertial' ):
        ax.streamplot( x1_mesh[index], x2_mesh[index], -Ey_cart[t][index], Ex_cart[t][index], density=d )
        if ( nc > 0 ):
            for ic in range(nc):
                xc = point_charges[t][ic,0]
                yc = point_charges[t][ic,1]
                ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'Stream lines of advection field at $t = %g$' %(t*dt), fontsize=fs_large )
    if ( frame == 'rotating' ):
        ax.streamplot( x1_mesh[index], x2_mesh[index], Ax_rot[index], Ay_rot[index], density=d )
        if ( nc > 0 ):
            for ic in range(nc):
                xc = point_charges[t][ic,0]
                yc = point_charges[t][ic,1]
                ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'Stream lines of advection field at $t^\prime = %g$' %(t*dt*tconv), fontsize=fs_large )
    ax.grid()
    ax.set_xlabel( r'$x$', fontsize=fs_large )
    ax.set_ylabel( r'$y$', fontsize=fs_large, rotation=0 )
    ax.tick_params( axis='both', labelsize=fs_small )
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )
