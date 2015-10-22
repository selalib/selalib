# coding: utf8
"""
Find and replace a full name in a Fortran library

Modules required
----------------
  * Built-in  : os, sys, argparse
  * Library   : maintenance_tools

"""
#
# Author: Yaman Güçlü, Oct 2015 - IPP Garching
#
# Last revision: 19 Oct 2015
#

__all__ = ['find_and_replace','main']
__docformat__ = 'reStructuredText'

#==============================================================================
# PARAMETERS.. CHANGE THIS
#==============================================================================

def ignore_dir( d ):
    """ Return True if subdirectory should be ignored.
    """
    return d in []


def select_file( f ):
    """ Return True if filename should be selected for processing.
    """
    is_fortran = f.endswith( '.F90' )
    is_header  = f.endswith( '.h' )
    is_cmake   = (f == 'CMakeLists.txt')
    return is_fortran or is_header or is_cmake

#==============================================================================
# MAIN FUNCTION
#==============================================================================

def find_and_replace( root, old_name, new_name ):
    """
    """
    import os, shutil
    import maintenance_tools as mtools

    for dirpath, dirnames, filenames in os.walk( root ):

        # Prune subdirectories that are not of interest
        for d in dirnames:
            if ignore_dir( d ):
                dirnames.remove( d )

        # Create greedy algorithm to perform name substitutions
        engine = mtools.RegexSubEngine( [old_name], [new_name] )

        # Search all files and replace full names
        for filename in filenames:

            # Ignore files
            if not select_file( filename ):
                continue

            filepath = os.path.join( dirpath, filename )

            # Search if file needs to be changed
            replace = False
            with open( filepath, 'r' ) as f:
                for line_no, old_line in enumerate( f ):
                    new_line = engine.process( old_line )
                    if new_line != old_line:
                        print( "  Found match in file '%s', at line %s" \
                                % (filepath, line_no+1) )
                        replace = True

            # Make necessary substitutions
            if replace:
                text = open( filepath, 'r' ).read()
                open( filepath, 'w' ).write( engine.process( text ) )

    # DONE
    print( "DONE" )

    # Return dictionary with local variables (useful for interactive work)
    return locals()

#==============================================================================
# PARSER
#==============================================================================

def parse_input():

  import argparse, sys

  parser = argparse.ArgumentParser (
      prog        = 'python3 '+ sys.argv[0],
      description = 'Find/replace full name in a Fortran library.',
      epilog      = ' ',
      formatter_class = argparse.ArgumentDefaultsHelpFormatter,
      )

  parser.add_argument( metavar = 'ROOT',
                       dest    = 'root',
                       help    = 'relative path of the root directory' )

  parser.add_argument( metavar = 'OLD_NAME',
                       dest    = 'old_name',
                       help    = 'full name to be searched for' )

  parser.add_argument( metavar = 'NEW_NAME',
                       dest    = 'new_name',
                       help    = 'replacement text' )

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

    # Walk directory tree and change names of all library modules
    find_and_replace( args.root, args.old_name, args.new_name )

#------------------------------------------------------------------------------
if __name__ == '__main__':
    # Run as main program
    main()
