# run interactively ("-i") after running sim_bsl_gc_2d0v_smooth_polar_splines.py

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

# to plot fine grid
nr = n1 // 8

# color map for plots
cmap = 'jet'

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

# equilibrium density
def plot_equilibrium_density( save=False ):
    fg = plt.figure(figsize=[9.0,9.0])
    ax = fg.add_subplot(111)
    #cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )
    clevels = np.linspace( min_rho_eq, max_rho_eq, 101 )
    im = ax.contourf( x1, x2, rho_eq, clevels, cmap=cmap )
    for c in im.collections:
        c.set_edgecolor('face')
    ax.set_title( r'Equilibrium density $\rho_{eq}$' )
    ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
    ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=0.5 )
    ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
    ax.set_xlabel( r'$x$' )
    ax.set_ylabel( r'$y$', rotation=0 )
    ax.set_aspect( 'equal' )
    cb = fg.colorbar( im, ax=ax )
    cb.formatter.set_powerlimits((0,0))
    cb.formatter.useMathText = True
    cb.update_ticks()
    fg.tight_layout()
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )

# point charge trajectory
def plot_point_charge_trajectory( frame, tmax, save=False ):
    fg = plt.figure(figsize=[9.0,9.0])
    ax = fg.add_subplot(111)
    if ( frame == 'inertial' ):
        xc = np.fromiter( ( point_charges[t][0,0] for t in range(tmax+1) ), dtype=float )
        yc = np.fromiter( ( point_charges[t][0,1] for t in range(tmax+1) ), dtype=float )
        ax.plot( xc, yc, '-', color='k', markersize=1. )
        ax.set_title( 'Point charge trajectory (inertial frame)' )
    if ( frame == 'rotating' ):
        xc = np.fromiter( ( point_charges[t][0,0] for t in range(tmax+1) ), dtype=float )
        yc = np.fromiter( ( point_charges[t][0,1] for t in range(tmax+1) ), dtype=float )
        th = np.fromiter( ( -omega*(t*dt) for t in range(tmax+1) ), dtype=float )
        xr, yr = rotate( xc, yc, th )
        ax.plot( xr, yr, '-', color='k', markersize=1. )
        ax.set_title( 'Point charge trajectory (rotating frame)' )
    ax.plot( x1[:,:], x2[:,:], color='lightgrey', lw=0.5 )
    ax.plot( x1.transpose()[:,:], x2.transpose()[:,:], color='lightgrey', lw=0.5 )
    ax.set_xlabel( r'$x$' )
    ax.set_ylabel( r'$y$', rotation=0 )
    ax.set_aspect( 'equal' )
    fg.tight_layout()
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )

# time evolution of radial position (rotating frame)
def plot_point_charge_coordinates( tmax, save=False ):
    fg = plt.figure(figsize=[9.0,9.0])
    ax1 = fg.add_subplot(121)
    ax2 = fg.add_subplot(122)
    xc = np.fromiter( ( point_charges[t][0,0] for t in range(tmax+1) ), dtype=float )
    yc = np.fromiter( ( point_charges[t][0,1] for t in range(tmax+1) ), dtype=float )
    th = np.fromiter( ( -omega*(t*dt) for t in range(tmax+1) ), dtype=float )
    xr, yr = rotate( xc, yc, th )
    rr = np.sqrt( xr**2 + yr**2 )
    tr = np.arctan2( yr, xr ) % ( 2.0 * np.pi )
    ax1.plot( th*(-tconv/omega), rr, '-', color='k', markersize=1. )
    ax2.plot( th*(-tconv/omega), tr, '-', color='k', markersize=1. )
    ax1.set_xlabel( r'$T$' )
    ax1.set_ylabel( r'$r$', rotation=0 )
    ax1.set_title( 'Radial position of point charge (rotating frame)' )
    ax2.set_xlabel( r'$T$' )
    ax2.set_ylabel( r'$\theta$', rotation=0 )
    ax2.set_title( 'Angular position of point charge (rotating frame)' )
    fg.tight_layout()
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )

# time shot of density and point charge
def plot_point_charge_time_shot( frame, t, save=False ):
    fg = plt.figure(figsize=[9.0,9.0])
    ax = fg.add_subplot(111)
    #cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )
    clevels = np.linspace( min_rho, max_rho, 101 )
    if ( frame == 'inertial' ):
        im = ax.contourf( x1, x2, rho[t], clevels, cmap=cmap )
        for c in im.collections:
            c.set_edgecolor('face')
        if ( nc > 0 ):
            for ic in range(nc):
                xc = point_charges[t][ic,0]
                yc = point_charges[t][ic,1]
                ax.plot( xc, yc, marker='o', color='w', markersize=5. )
        ax.set_title( r'Density $\rho$ at $T = %g$ (inertial frame)' %(t*dt*tconv) )
    if ( frame == 'rotating' ):
        x1r, x2r = rotate( x1, x2, -omega*(t*dt) )
        im = ax.contourf( x1r, x2r, rho[t], clevels, cmap=cmap )
        for c in im.collections:
            c.set_edgecolor('face')
        if ( nc > 0 ):
            for ic in range(nc):
                xc = point_charges[t][ic,0]
                yc = point_charges[t][ic,1]
                xr, yr = rotate( xc, yc, -omega*(t*dt) )
                ax.plot( xr, yr, marker='o', color='w', markersize=5. )
        ax.set_title( r'Density $\rho$ at $T = %g$ (rotating frame)' %(t*dt*tconv) )
    ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
    ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=0.5 )
    ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
    ax.set_xlabel( r'$x$' )
    ax.set_ylabel( r'$y$', rotation=0 )
    ax.set_aspect( 'equal' )
    cb = fg.colorbar( im, ax=ax )
    cb.formatter.set_powerlimits((0,0))
    cb.formatter.useMathText = True
    cb.update_ticks()
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )

# animated plots
def plot_time_evolution( arg, tmax ):

    fg = plt.figure(figsize=[9.0,9.0])
    ax = fg.add_subplot(111)
    #cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )

    for t in (ti for ti in tt if ti <= tmax):
        ax.clear()
        #cax.clear()
        if ( arg == 'rho' ):
            clevels = np.linspace( min_rho, max_rho, 101 )
            im = ax.contourf( x1, x2, rho[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            if ( nc > 0 ):
                for ic in range(nc):
                    xc = point_charges[t][ic,0]
                    yc = point_charges[t][ic,1]
                    ax.plot( xc, yc, marker='o', color='k', markersize=5. )
            ax.set_title( r'Density $\rho$ at $T = %g$ (inertial frame)' %(t*dt*tconv) )
        if ( arg == 'rho_rot' ):
            clevels = np.linspace( min_rho, max_rho, 101 )
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
            ax.set_title( r'Density $\rho$ at $T = %g$ (rotating frame)' %(t*dt*tconv) )
        if ( arg == 'delta_rho' ):
            clevels = np.linspace( min_delta_rho, max_delta_rho, 101 )
            im = ax.contourf( x1, x2, rho[t]-rho_eq, clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            ax.set_title( r'Density $\rho$ at $T = %g$' %(t*dt*tconv) )
        if ( arg == 'phi' ):
            clevels = np.linspace( min_phi, max_phi, 101 )
            im = ax.contourf( x1, x2, phi[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            ax.set_title( r'Potential $\phi$ at $T = %g$' %(t*dt*tconv) )
        if ( arg == 'Ex' ):
            clevels = np.linspace( min_Ex, max_Ex, 101 )
            im = ax.contourf( x1, x2, Ex[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            ax.set_title( r'Electric field component $E_x$ at $T = %g$' %(t*dt*tconv) )
        if ( arg == 'Ey' ):
            clevels = np.linspace( min_Ey, max_Ey, 101 )
            im = ax.contourf( x1, x2, Ey[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            ax.set_title( r'Electric field component $E_y$ at $T = %g$' %(t*dt*tconv) )
        if ( arg == 'Em' ):
            clevels = np.linspace( min_Em, max_Em, 101 )
            im = ax.contourf( x1, x2, Em[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            ax.set_title( r'Electric field magnitude $|E|$ at $T = %g$' %(t*dt*tconv) )
        ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
        ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=0.5 )
        ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
        ax.set_xlabel( r'$x$' )
        ax.set_ylabel( r'$y$', rotation=0 )
        ax.set_aspect( 'equal' )
        cb = fg.colorbar( im, ax=ax )
        cb.formatter.set_powerlimits((0,0))
        cb.formatter.useMathText = True
        cb.update_ticks()
        fg.canvas.draw()
        plt.pause(0.1)

## Limit t in some range
#mask = ( xt >= 0. )
#xt = xt[mask]
#y1 = y1[mask]
#ya = ya[mask]

# time evolution of scalar diagnostics
def plot_scalar_data( arg, save=False ):

    fg = plt.figure()
    ax = fg.add_subplot(111)
    fs = 12
    if ( arg == 'mass' ):
        # optimal size [6.4,2.4]
        ax.plot( xt, (y1[0]-y1[:])/y1[0], '-r' , lw=1.5 )
        ax.set_ylabel( r'$\mathcal{E}_\mathcal{M}$', rotation=0, fontsize=fs )
    if ( arg == 'energy' ):
        # optimal size [6.4,2.4]
        ax.plot( xt, (y2[0]-y2[:])/y2[0], '-r' , lw=1.5 )
        ax.set_ylabel( r'$\mathcal{E}_\mathcal{W}$', rotation=0, fontsize=fs )
    if ( arg == 'l2_norm' ):
        # optimal size [9.5,6.0]
        ax.plot( xt, y3, '-r' , lw=1.5, label=r'Numerical instability (mode $m=%d$)' %(l) )
        ax.plot( xt, ya, '--k', lw=1.5, label=r'Analytical growth rate: Im$(\omega)\approx %1.2g$' %(gamma) )
        ax.set_ylabel( r'$||\phi-\phi_0||_{L^2}$', fontsize=fs )
        ax.set_yscale( 'log' )
        ymin = ax.get_ylim()[0]
        ymax = 2.*max(y3)
        ax.set_ylim( [ymin,ymax] )
        ax.legend( loc='best', fontsize=fs )
    ax.set_xlabel( r'$t$', fontsize=fs )
    ax.ticklabel_format( axis='y', style='sci', scilimits=(0,0), useMathText=True )
    ax.tick_params( labelsize=fs )
    ax.grid()
    fg.tight_layout()
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )

# vorticity, shear, expansion
def plot_advection_field( arg, tmax ):

    fg = plt.figure(figsize=[9.0,9.0])
    ax = fg.add_subplot(111)
    #cax = make_axes_locatable(ax).append_axes( 'right', size='8%', pad='5%' )

    for t in (ti for ti in tt if ti <= tmax):
        ax.clear()
        #cax.clear()
        if ( arg == 'vorticity' ):
            clevels = np.linspace( min_vorticity, max_vorticity, 101 )
            im = ax.contourf( x_meshgrid, y_meshgrid, vorticity[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            ax.set_title( r'Vorticity of advection field at $T = %g$' %(t*dt*tconv) )
        if ( arg == 'shear' ):
            clevels = np.linspace( min_shearrate, max_shearrate, 101 )
            im = ax.contourf( x_meshgrid, y_meshgrid, shearrate[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            ax.set_title( r'Shear rate of advection field at $T = %g$' %(t*dt*tconv) )
        if ( arg == 'expansion' ):
            clevels = np.linspace( min_expansion, max_expansion, 101 )
            im = ax.contourf( x_meshgrid, y_meshgrid, expansion[t], clevels, cmap=cmap )
            for c in im.collections:
                c.set_edgecolor('face')
            ax.set_title( r'Expansion rate of advection field at $T = %g$' %(t*dt*tconv) )
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
        cb = fg.colorbar( im, ax=ax )
        cb.formatter.set_powerlimits((0,0))
        cb.formatter.useMathText = True
        cb.update_ticks()
        fg.canvas.draw()
        plt.pause(0.1)

# stream lines (works for circle only)

index_x1 = return_index( x1_cart,  0.3, 0.5 )
index_x2 = return_index( x2_cart, -0.1, 0.1 )

index = tuple( ( index_x2, index_x1 ) )

def plot_stream_lines_global( frame, save=False ):
    t = 0
    fg = plt.figure(figsize=[9.0,9.0])
    ax = fg.add_subplot(111)
    if ( frame == 'inertial' ):
        ax.streamplot( x1_mesh, x2_mesh, -Ey_cart[t], Ex_cart[t], density=6. )
        if ( nc > 0 ):
            for ic in range(nc):
                xc = point_charges[t][ic,0]
                yc = point_charges[t][ic,1]
                ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'Stream lines of advection field at $T = %g$ (inertial frame)' %(t*dt*tconv) )
    if ( frame == 'rotating' ):
        ax.streamplot( x1_mesh, x2_mesh, Ax_rot, Ay_rot, density=6. )
        if ( nc > 0 ):
            for ic in range(nc):
                xc = point_charges[t][ic,0]
                yc = point_charges[t][ic,1]
                ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'Stream lines of advection field at $T = %g$ (rotating frame)' %(t*dt*tconv) )
    ax.plot( x1[:,::nr], x2[:,::nr], color='lightgrey', lw=0.5 )
    ax.plot( x1[:,n1-1], x2[:,n1-1], color='lightgrey', lw=0.5 )
    ax.plot( x1.transpose()[:,::nr], x2.transpose()[:,::nr], color='lightgrey', lw=0.5 )
    ax.set_xlabel( r'$x$' )
    ax.set_ylabel( r'$y$', rotation=0 )
    ax.set_aspect( 'equal' )
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )

def plot_stream_lines_local( frame, save=False ):
    t = 0
    fg = plt.figure(figsize=[9.0,9.0])
    ax = fg.add_subplot(111)
    if ( frame == 'inertial' ):
        ax.streamplot( x1_mesh[index], x2_mesh[index], -Ey_cart[t][index], Ex_cart[t][index], density=6. )
        if ( nc > 0 ):
            for ic in range(nc):
                xc = point_charges[t][ic,0]
                yc = point_charges[t][ic,1]
                ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'Stream lines of advection field at $T = %g$ (inertial frame)' %(t*dt*tconv) )
    if ( frame == 'rotating' ):
        ax.streamplot( x1_mesh[index], x2_mesh[index], Ax_rot[index], Ay_rot[index], density=6. )
        if ( nc > 0 ):
            for ic in range(nc):
                xc = point_charges[t][ic,0]
                yc = point_charges[t][ic,1]
                ax.plot( xc, yc, marker='o', color='k', markersize=5. )
        ax.set_title( r'Stream lines of advection field at $T = %g$ (rotating frame)' %(t*dt*tconv) )
    ax.grid()
    ax.set_xlabel( r'$x$' )
    ax.set_ylabel( r'$y$', rotation=0 )
    fg.show()
    if ( save ):
       fg_name = 'figure.pdf'
       fg.savefig( fg_name, dpi=300 )
