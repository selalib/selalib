# coding: utf8
"""
Recursively search directory tree for UTF-8 errors.

Modules required
----------------
  * Built-in  : os, sys, argparse
  * Library   : maintenance_tools

"""
#
# Author: Yaman Güçlü, Aug 2016 - IPP Garching
#
# Last revision: 31 Aug 2016
#

__all__ = ['test_utf8_errors','main']
__docformat__ = 'reStructuredText'

#==============================================================================
# PARAMETERS.. CHANGE THIS
#==============================================================================

def ignore_dir( d ):
    """ Return True if subdirectory should be ignored.
    """
    return d in ['__pycache__']

def select_file( f ):
    """ Return True if a file should be selected for processing.
    """
    import os
    ext = os.path.splitext( f )[1]
    return ext in ['.py','.cgf','.h','.f','.f90','.F90']

#==============================================================================
# MAIN FUNCTION
#==============================================================================

def test_utf8_errors( root, verbose=False ):
    """ Recursively search 'root' directory tree for UTF8 errors in files.
    """
    import os
    import maintenance_tools as mtools
    counter = 0
    print( "Function 'test_utf8_errors': walking over directory '%s'.." \
            % os.path.abspath( root ) )
    for f in mtools.recursive_file_search( root, ignore_dir, select_file ):
        if verbose: print( 'Read: '+f )
        n = mtools.find_utf8_error( f )
        if n > 0:
            print( "  . Error in line %d of file %s" % (n,f) )
            counter += 1
    if verbose: print( 'DONE' )
    if counter == 0: print( '\nPASSED: All files analysed are UTF-8' )
    else           : print( '\nWARNING: %s files are not UTF-8' % counter )

#==============================================================================
# PARSER
#==============================================================================

def parse_input():

  import argparse, sys

  parser = argparse.ArgumentParser (
      prog        = 'python3 '+ sys.argv[0],
      description = 'Recursively search for UTF-8 errors in the given directory',
      epilog      = ' ',
      formatter_class = argparse.ArgumentDefaultsHelpFormatter,
      )

  parser.add_argument( metavar = 'ROOT',
                       dest    = 'root',
                       help    = 'relative path of the root directory' )

  parser.add_argument( '-v', '--verbose',
                       action = 'store_true',
                       help   = 'increase output verbosity' )

  return parser.parse_args()

#==============================================================================
# SCRIPT FUNCTIONALITY
#==============================================================================

def main():

    # Parse input arguments
    print('')
    args = parse_input()
    print(args)
    print('')

    # Walk directory tree and only search for UTF8 errors
    test_utf8_errors( args.root, args.verbose )

#------------------------------------------------------------------------------
if __name__ == '__main__':
    # Run as main program
    main()
