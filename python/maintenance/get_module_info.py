# coding: utf8
"""
Extract info from pre-processed Fortran module file.

Modules required
----------------
  * Built-in  : os, sys, argparse
  * Library   : sll2py
  * 3rd-party : f2py.fparser

"""
#
# Author: Yaman Güçlü, Oct 2015 - IPP Garching
#
# Last revision: 20 Nov 2015
#
from __future__ import print_function

__all__ = ['get_module_info','main']
__docformat__ = 'reStructuredText'

#==============================================================================
# PARAMETERS.. CHANGE THIS
#==============================================================================

def select_file( f ):
    """ Return True if filename should be selected for processing.
    """
    is_preprocessed_fortran = f.endswith( '.f90' )
    return is_preprocessed_fortran

#==============================================================================
# UTILITY FUNCTIONS (needed for Python2, as 'maintenance_tools' is for Python3)
#==============================================================================

def find_module_def( filename ):
    mod_def_list = []
    with open( filename, 'r' ) as f:
        for line in f:
            spl = line.partition( '!' )[0].split()
            if len( spl ) > 1:
                if spl[0] == 'module' and spl[1] != 'procedure':
                    mod_def_list.append( spl[1] )
    return mod_def_list

def find_program_name( filename ):
    program_list = []
    with open( filename, 'r' ) as f:
        for line in f:
            spl = line.partition( '!' )[0].split()
            if len( spl ) > 1:
                if spl[0] == 'program':
                    program_list.append( spl[1] )
    return program_list

#==============================================================================
# MAIN FUNCTION
#==============================================================================

def get_module_info( *filepaths ):
    """ Extract info from pre-processed Fortran module file.
    """
    import os
    from sll2py.fortran_module import FortranModule
    from sll2py.fparser.api    import parse

    for filepath in filepaths:

        # Split path
        dirpath, filename = os.path.split( filepath )

        # Ignore files
        if not select_file( filename ):
            continue

        # Determine if only one module is present, otherwise ignore file
        num_mods = len( find_module_def( filepath ) )
        if len( find_program_name( filepath ) ) > 0:
            skip_file = True
            print( "WARNING: Fortran program. ", end='' )
        elif num_mods > 1:
            skip_file = True
            print( "WARNING: multiple module definitions. ", end='' )
        elif num_mods == 0:
            skip_file = True
            print( "WARNING: no modules or programs in file. ", end='' )
        else:
            skip_file = False
        if skip_file:
            print( "Skipping file '%s'" % filename )
            continue

        # Create fparser module object and print info
        tree = parse( filepath, analyze=False )
        fmod = tree.content[0]
        print( "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" )
        print( "Module '%s':" % fmod.name )
        print( "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" )
        print()
        print( "=========" )
        print( "STRUCTURE" )
        print( "=========" )
        for item in fmod.content:
            print( "%s" % type( item ) )
        print()

        # Create 'my' module object and print info
        mmod = FortranModule( filepath, fmod )
        mmod.show()

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
      prog        = 'python '+ sys.argv[0],
      description = 'Extract info from pre-processed Fortran module files.',
      epilog      = ' ',
      formatter_class = argparse.ArgumentDefaultsHelpFormatter,
      )

  parser.add_argument( metavar = 'FILE',
                       nargs   = '+',
                       dest    = 'filepaths',
                       help    = 'relative path of the Fortran module files' )

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
    get_module_info( *args.filepaths )

#------------------------------------------------------------------------------
if __name__ == '__main__':
    # Run as main program
    main()
