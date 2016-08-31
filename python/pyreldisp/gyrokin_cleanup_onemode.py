from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from zealpy_gyrokin_cleanup_anal import zealpy_gyrokin_cleanup_anal

#------------------------------------
# INPUT DATAS (defined as GYSELA)
#------------------------------------
Zi           = 1
#--> For geometry definition
NNr          = 100
rhomin       = 0.006896551724 #0.00001
rhomax       = 1.
minor_radius = 14.5
aspect_ratio = 16.53849335 #16.5
#--> For analytical density and temperature profiles definition 
kappaTi  = 66.
kappaTe  = 66.
kappan   = 13.2
deltarTi = 0.1
deltarTe = 0.1
deltarn  = 0.2
#--> Choice of the box where zeros will be search
xmin = -1.
xmax = 1.
ymin = 0.0001
ymax = 0.01 #0.1
iota = 0
B0   = -1.

#------------------------------------------------------------
# Normalization to obtain the same profiles than in GYSELA
#------------------------------------------------------------
#--> Geometry definition
R0   = minor_radius*aspect_ratio
rmin = rhomin*minor_radius
rmax = rhomax*minor_radius
Lr   = np.abs(rmax-rmin)
Lz   = 2.*np.pi*R0
#--> Same normalization than GYSELA for profiles 
deltarTi = deltarTi*Lr
deltarTe = deltarTe*Lr
deltarn  = deltarn*Lr
invLTi   = kappaTi/R0
invLTe   = kappaTe/R0
invLn    = kappan/R0
iota     = iota/R0

def print_func( name, val ):
    print( '{:8s} = {:.15g}'.format( name, val ) )

print_func( 'R0'  , R0 )
print_func( 'rmin', rmin )
print_func( 'rmax', rmax )
print_func( 'Lr'  , Lr )
print_func( 'Lz'  , Lz )
print_func( 'deltarTi', deltarTi )
print_func( 'deltarTe', deltarTe )
print_func( 'deltarn' , deltarn )
print_func( 'invLTi'  , invLTi )
print_func( 'invLTe'  , invLTe )
print_func( 'invLn'   , invLn )
print_func( 'iota', iota )
print_func( 'B0'  , B0 )

#------------------------------------------------------------
# Dispersion relation solving
#------------------------------------------------------------
#--> Constructor for using dispersion relation for 
#-->   gyrokinetic with analytic n and T profiles
zp = zealpy_gyrokin_cleanup_anal( Zi=Zi,
                NNr=NNr, rmin=rmin, Lr=Lr,
                invLTi=invLTi,deltarTi=deltarTi,
                invLTe=invLTe,deltarTe=deltarTe,
                invLn0=invLn, deltarn0=deltarn,
                iota=iota, B0=B0 )

#--> Choice of the mode (m,n)
print( '\n------------------' )
mm_choice = raw_input( 'Poloidal mode m = ' )
nn_choice = raw_input( 'Toroidal mode n = ' )

ifig = 1
while mm_choice and nn_choice:
    zp.zeros = []
    mm  = int(mm_choice)
    nn  = int(nn_choice)

    print( '\nSearch zeros of mode (m,n) = (%d,%d)'% (mm,nn) )
    kpar    = float(nn) * 2*np.pi/Lz
    ktheta  = float(mm)

    #---> Initialisation of the dispersion relation for (m,n) mode
    zp.init_rel_disp( ktheta=ktheta, kpar=kpar )

    #---> Find all the zeros contained is the box defined before
    zp.get_zeros( xmin, xmax, ymin, ymax )

    #---> Select the most unstable modes
    if len(zp.zeros) != 0:
        imag_zeros   = np.imag( zp.zeros )
        mostunstable = zp.zeros[np.argmax(imag_zeros)]
        print( 'Zero with largest imaginary part: ', mostunstable )
        print( 'All zeros found:' )
        print( *zp.zeros, sep='\n', end='\n\n' )

        #---> Computation of the eigenvector (phi) corresponding to the
        #--->  most unstable mode
        zp.function_input_module.kernel( mostunstable )
        #-----> Rescaling of the eigenvector to get phi from phi
        rmesh = zp.params.dic['rmesh']
        n0    = zp.params.dic['n0']
        phi_eigvect = zp.function_input_module.vector[:,0] / np.sqrt( rmesh*n0 )

        #--> Figure 1: zeros of dispersion relation
        fig = plt.figure( ifig )
        ax  = fig.add_subplot( 1,1,1 )
        ax.plot( np.real(zp.zeros), np.imag(zp.zeros), '.',
                  label = 'm={} n={}'.format( mm, nn ) )
        ax.axis( [xmin,xmax,ymin,ymax] )
        ax.set_title( 'Zeros of dispersion relation' )
        ax.set_xlabel( r'Re($\omega$)', size=14 )
        ax.set_ylabel( r'Im($\omega$)', size=14, rotation='horizontal' )
        ax.legend()
        ax.grid()
        fig.show()
        ifig += 1

        #--> Figure 2: Plot radial profile of the most unstable mode
        fig = plt.figure( ifig )
        ax  = fig.add_subplot( 1,1,1 )
        ax.plot( rmesh, abs(phi_eigvect) )
        ax.set_title( 'Profile of most unstable mode' )
        ax.set_xlabel( 'r' )
        ax.grid()
        fig.show()
        ifig += 1

        #---> Display figures at screen
        if plt.rcParams['backend'] in ['TkAgg','WXAgg']:
            plt.show()
    else:
        print( '---> BECAREFUL: Zeros not found for (m,n) = (%d,%d)\n' % (mm,nn) )

    #--> Choice of the mode (m,n)
    print( '\n------------------' )
    mm_choice = raw_input( 'Poloidal mode m = ' )
    nn_choice = raw_input( 'Toroidal mode n = ' )

raw_input( "\nPress enter to quit:\n" )
