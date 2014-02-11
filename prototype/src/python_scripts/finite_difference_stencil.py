#!/usr/bin/env python3
"""
SYNOPSIS

    finite_difference_stencil [-h,--help] [-v,--verbose] [--version]

DESCRIPTION

    This script computes the coefficients needed to estimate the value of
    a discrete function or any of its derivatives by providing: the degree of
    the derivative desired, the order of accuracy desired (power of the cell-
    spacing) and a list of indices of the points to be used, where the index
    zero represents the point of interest for the evaluation.

EXAMPLES

    ./finite_difference_stencil.py list-of-points order-of-accuracy
    ./finite_difference_stencil.py 2 2 -1 0 1  would compute the SECOND
    derivative, with an accuracy of h^2, using the left neighbor (-1), the 
    point of interest (0) and the right neighbor (+1). The result is the
    well-known [ 1, -2, 1 ] stencil.

EXIT STATUS

    TODO:

AUTHOR

    Edwin CHACON-GOLCHER <chacongolcher@math.unistra.fr>

LICENSE

    

VERSION

    1.0
"""

import sys, os, traceback, argparse
import time
import re
import math as m
import numpy

def main ():
    global options, parser
    parser   = argparse.ArgumentParser( )
    rank     = 0
    i        = 0
    coeffs   = [ ]
    det      = 0
    dists    = [ ]
    parser.add_argument('derivative_degree',
                        help="desired derivative degree at reference point",
                        type=int)
    parser.add_argument('order_desired',
                        help="highest order of the Taylor expansion terms kept",
                        type=int)
    parser.add_argument('points_for_stencil', 
                        help="list of points for which the coefficients are needed, given in terms of their relative distance from the reference point (where the answer is sought) in units of the cell spacing. Example: -1 0 1 would request that the coefficients be given for the neighboring points and the reference point itself.",
                        nargs='+', type=float)

    args = parser.parse_args()
    print('arguments passed:')
    print(args)
    rank = args.order_desired + 1
    deg  = args.derivative_degree
    dists = args.points_for_stencil
    # build linear system Ax = b
    # build the 'b' vector:
    b = [ 0 if i != deg else 1 for i in range(0,rank)]
    # build the 'A' matrix:
    A = [ [(i**row)/m.factorial(row) for i in dists] for row in range(0,rank)]
    # print(A)
    coeffs = numpy.linalg.solve(A,b)
    print(coeffs)


if __name__ == '__main__':
    try:
        start_time = time.time()
        main()
        sys.exit(0)
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
    except SystemExit as e: # sys.exit()
        raise e
    except Exception as e:
        print( 'ERROR, UNEXPECTED EXCEPTION')
        print( str(e))
        traceback.print_exc()
        os._exit(1)


