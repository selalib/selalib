import numpy as np
import matplotlib.pyplot as plt

#==============================================================================
# INPUT
#==============================================================================

# Domain for f(x,v)
x = np.linspace( 0,13,300 ) 
v = np.linspace( -4.5,4.5,200 )

shading = 'flat'
#shading = 'goraud'

# Number of snapshots
Nfiles = 200

#==============================================================================
# FUNCTIONS
#==============================================================================

def get_file_name( j ):
    s = 'distribution__{:4d}.xmf'.format( j ).replace( ' ','0' )
    return s

def get_fig_title( j ):
    n = '{:4d}'.format( j ).replace( ' ','0' )
    s = 'f(x,v) from dataset #' + n
    return s

def update_plot( quad, j, shading='flat' ):
    file_name = get_file_name( j )
    fig_title = get_fig_title( j )
    f = np.genfromtxt( file_name, skip_header=16, skip_footer=5 )
    if shading == 'flat':
        f = f.reshape( [Nv,Nx] )[:-1,:-1].ravel()
    elif shading == 'gouraud':
        pass
    else:
        raise( NotImplementedError )
    quad.set_array( f )
    quad.get_axes().set_title( fig_title )
    quad.get_figure().canvas.draw()
    plt.pause( 0.01 )

#==============================================================================
# SCRIPT FUNCTIONALITY
#==============================================================================

Nx = len( x )
Nv = len( v )

# Create empty figure
fig = plt.figure()
ax  = fig.add_subplot(111)
ax.set_xlabel('x')
ax.set_ylabel('v')
ax.set_xlim( x[0],x[-1] )
ax.set_ylim( v[0],v[-1] )
fig.show()

# Create first plot (initial conditions?)
file_name = get_file_name( 1 )
fig_title = get_fig_title( 1 )
f = np.genfromtxt( file_name,skip_header=16,skip_footer=5 ).reshape( [Nv,Nx] )
quad = ax.pcolormesh( x,v,f )
ax.set_title( fig_title )
plt.colorbar( quad )
fig.canvas.draw()
plt.pause( 0.05 )

# Cycle over all files
for j in range( 2,Nfiles ):
    update_plot( quad, j )

# Keep window open
plt.show( block=False )

