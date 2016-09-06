from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from zealpy_gyrokin_cleanup_anal import zealpy_gyrokin_cleanup_anal

try :
    raw_input
except NameError:
    raw_input = input

#===============================================================================
# INITIALIZATION
#===============================================================================

# Input data defined as GYSELA
Zi = 1

# Geometry
#---------
rhomin       = 0.006896551724 #0.00001
rhomax       = 1.
minor_radius = 14.5
aspect_ratio = 16.53849335 #16.5

# Magnetic field
#---------------
# B0     : Magnetic field at axis of cylinder
# iota(r): Radial profile of rotational transform
B0   = -1.
iota = lambda r: 0.8

# Numerical parameters
#---------------------
# NNr                  : Number of cells in radial direction
# [xmin,xmax,ymin,ymax]: Box in complex plane where zeros will be searched
NNr  = 100
xmin = -1.
xmax = 1.
ymin = 0.0001
ymax = 0.01 #0.1

# Analytic density and temperature profiles
#------------------------------------------
kappaTi  = 66.
kappaTe  = 66.
kappan0  = 13.2
deltarTi = 0.1
deltarTe = 0.1
deltarn0 = 0.2

# Normalization to obtain the same profiles of GYSELA
R0   = minor_radius*aspect_ratio
rmin = rhomin*minor_radius
rmax = rhomax*minor_radius
Lr   = np.abs(rmax-rmin)
Lz   = 2.*np.pi*R0
deltarTi = deltarTi*Lr
deltarTe = deltarTe*Lr
deltarn0 = deltarn0*Lr
invLTi   = kappaTi/R0
invLTe   = kappaTe/R0
invLn0   = kappan0/R0

# Print geometry and equilibrium parameters
def print_params():

    print( '\nParameters:\n' )
  
    def print_func( name, val ):
        print( '{:8s} = {:.14g}'.format( name, val ) )
     
    print_func( 'R0', R0 )
    print_func( 'rmin', rmin )
    print_func( 'rmax', rmax )
    print_func( 'Lr', Lr )
    print_func( 'Lz', Lz )

    print()
    print_func( 'deltarTi', deltarTi )
    print_func( 'deltarTe', deltarTe )
    print_func( 'deltarn0', deltarn0 )
    print_func( 'invLTi', invLTi )
    print_func( 'invLTe', invLTe )
    print_func( 'invLn0', invLn0 )
    print_func( 'B0', B0 )


#===============================================================================
# DISPERSION RELATION SOLVER
#===============================================================================

def analyse_mode( zp, mm, nn, verbose=False ):
    zp.mm    = mm
    zp.nn    = nn
    zp.zeros = []

    ktheta = float(mm)
    kzeta  = float(nn)*(2.0*np.pi/Lz)

    # Initialization of dispersion relation for (m,n) mode
    zp.init_rel_disp( ktheta=ktheta, kzeta=kzeta )

    # Find all zeros contained in the box defined before
    zp.get_zeros( xmin, xmax, ymin, ymax )

    # Select the most unstable modes
    if len(zp.zeros) != 0:
        imag_zeros   = np.imag( zp.zeros )
        zp.mostunstable = zp.zeros[np.argmax(imag_zeros)]
        if verbose:
            print( 'All zeros found:' )
            print( *zp.zeros, sep='\n', end='\n\n' )
            print( 'Zero with largest imaginary part: ')
            print( 'omega =', zp.mostunstable, '\n' )

        # Computation of the eigenvector (psi) corresponding
        # to the most unstable mode
        zp.function_input_module.kernel( zp.mostunstable )
        # Rescaling of the eigenvector to get phi from psi
        # (see initial comment in 'fun_gyrokin.f90')
        rmesh = zp.params.dic['rmesh']
        n0    = zp.params.dic['n0']
        zp.phi_eigvect = zp.function_input_module.vector[:,0]/np.sqrt(rmesh*n0)

    else:
        if verbose:
            print( '---> BECAREFUL <---' )
            print( 'Zeros not found for  (m,n) = (%d,%d)\n' % (mm,nn) )


def manual_scan( zp ):
    zeros_dict = OrderedDict()
    while True:
        # Choice of the mode (m,n)
        print( '\n-------------' )
        mm_choice = raw_input( 'Poloidal mode m = ' )
        nn_choice = raw_input( 'Toroidal mode n = ' )
        print()
        if mm_choice and nn_choice:
            mm = int(mm_choice)
            nn = int(nn_choice)
            analyse_mode( zp, mm, nn, verbose=True )
            zeros_dict[(mm,nn)] = np.array( zp.zeros )
        else:
            break
    return zeros_dict


def automatic_scan( zp, mm_list, nn_list ):
    zeros_dict = OrderedDict()
    for nn in nn_list:
        for mm in mm_list:
            print( '\n-------------' )
            print( 'Poloidal mode m = %d' % mm )
            print( 'Toroidal mode n = %d' % nn )
            print()
            analyse_mode( zp, mm, nn, verbose=True )
            zeros_dict[(mm,nn)] = np.array( zp.zeros )
    print()
    return zeros_dict


def find_most_unstable_mode( zp, zeros_dict ):

    max_mode = (None,None)
    max_real = 0.0
    max_imag = 0.0
    for (mm,nn), zeros in zeros_dict.items():
        max_imag_loc = max( zeros.imag ) 
        if max_imag_loc > max_imag:
            max_mode = (mm,nn)
            max_real = zeros.real[np.argmax(zeros.imag)]
            max_imag = max_imag_loc

    return max_mode


#===============================================================================
# PLOTS
#===============================================================================

def plot_profiles( axes, zp ):

    rmesh = zp.params.dic['rmesh']
    Ti    = zp.params.dic['Ti']
    Te    = zp.params.dic['Te']
    n0    = zp.params.dic['n0']

    axes[0].plot( rmesh, Ti  )
    axes[0].set_title( 'Equilibrium ion temperature profile' )
    axes[1].plot( rmesh, Te )
    axes[1].set_title( 'Equilibrium electron temperature profile' )
    axes[2].plot( rmesh, n0 )
    axes[2].set_title( 'Equilibrium density profile (not normalized)' )
    for ax in axes:
        ax.set_xlabel('r')
        ax.grid()


def plot_zeros( ax, zeros_dict ):

    for (mm,nn), zeros in zeros_dict.items():
        ax.plot( zeros.real, zeros.imag, '.',
                      label = r'$m={}$, $n={}$'.format( mm, nn ) )

    ax.axis( [xmin,xmax,ymin,ymax] )
    ax.set_title( 'Zeros of dispersion relation' )
    ax.set_xlabel( r'Re($\omega$)', size=14 )
    ax.set_ylabel( r'Im($\omega$)', size=14, rotation='horizontal' )
    ax.legend()
    ax.grid()


def plot_mode( ax, zp ):
    ax.plot( zp.params.dic['rmesh'], abs( zp.phi_eigvect ) )
    ax.set_title(
      r'Profile of most unstable mode ({},{}): $\omega={:.6g}$' \
              .format( zp.mm, zp.nn, zp.mostunstable ) )
    ax.set_xlabel( 'r' )
    ax.grid()


#===============================================================================
# PARSER
#===============================================================================

def parse_input():

  import argparse, sys

  parser = argparse.ArgumentParser (
      prog        = 'python3 '+ sys.argv[0],
      description = 'Compute growth rates for ITG instability in \
                     4D screw pinch model in cylindrical geometry.',
      epilog      = ' ',
      formatter_class = argparse.ArgumentDefaultsHelpFormatter )

  subparsers = parser.add_subparsers( dest='scan', help='' )
  
  # 'Manual' scan
  parser_a = subparsers.add_parser( 'manual',
                                    help='scan poloidal and toroidal modes \
                                          manually inserted by the user.' )

  # 'Automatic' scan
  parser_b = subparsers.add_parser( 'auto',
                                    help='scan poloidal and toroidal modes     \
                                          automatically within the range given \
                                          by the user.                         \
                                          Type auto -h for more help.' )
  parser_b.add_argument( '-m', '--m_range',
          type     = int,
          nargs    = 2,
          required = True,
          metavar  = ('m_min','m_max'),
          help     = 'range for m mode (along theta)' )

  parser_b.add_argument( '-n', '--n_range',
          type     = int,
          nargs    = 2,
          required = True,
          metavar  = ('n_min','n_max'),
          help     = 'range for n mode (along zeta )' )

  return parser.parse_args()


#===============================================================================
# MAIN
#===============================================================================

def main():

    args = parse_input()

    print_params()

    # Constructor for drift-kinetic dispersion relation
    # with analytic temperature and density profiles
    zp = zealpy_gyrokin_cleanup_anal( Zi=Zi,NNr=NNr,rmin=rmin,Lr=Lr,R0=R0,
                                      invLTi=invLTi,deltarTi=deltarTi,
                                      invLTe=invLTe,deltarTe=deltarTe,
                                      invLn0=invLn0,deltarn0=deltarn0,
                                      B0=B0, iota=iota )

    if args.scan == 'manual':
        zeros_dict = manual_scan( zp )
    elif args.scan == 'auto':
        m_min, m_max = args.m_range
        n_min, n_max = args.n_range
        mm_list = range( m_min, m_max+1 )
        nn_list = range( n_min, n_max+1 )
        zeros_dict = automatic_scan( zp, mm_list, nn_list )
    else:
        raise ValueError( args.scan )

    mm, nn = find_most_unstable_mode( zp, zeros_dict )
    print( '----------------------------------' )
    print( 'Results for the most unstable mode (m,n)=(%d,%d): \n' % (mm,nn) )
    analyse_mode( zp, mm, nn )

    # Plot equilibrium radial profiles of Ti, Te and n0
    fig1, axes1 = plt.subplots( nrows=1, ncols=3, sharex=True )
    plot_profiles( axes1, zp )
    fig1.show()

    # Plot zeros of dispersion relation and 
    # radial profile of the most unstable mode
    fig2, axes2 = plt.subplots( nrows=1, ncols=2 )
    plot_zeros( axes2[0], zeros_dict )
    plot_mode ( axes2[1], zp )
    fig2.show()

    # Wait for ENTER for stopping execution
    msg = '\nPress ENTER to quit'
    try:
        raw_input( msg )
    except NameError:
        input( msg )
    plt.close('all')


if __name__ == '__main__':
    main()
