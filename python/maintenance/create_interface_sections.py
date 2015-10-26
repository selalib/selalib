# coding: utf8
"""
Create interface sections for all modules (and programs?) in a Fortran library.

Modules required
----------------
  * Built-in  : os, sys, argparse
  * Library   : maintenance_tools

"""
#
# Author: Yaman Güçlü, Oct 2015 - IPP Garching
#
# Last revision: 23 Oct 2015
#
from __future__ import print_function

__all__ = ['create_interface_sections','main']
__docformat__ = 'reStructuredText'

#==============================================================================
# PARAMETERS.. CHANGE THIS
#==============================================================================

def ignore_dir( d ):
    """ Return True if subdirectory should be ignored.
    """
    non_source_dirs = ['bin','CMakeFiles','doc','include','modules','Testing']
    useless_subdirs = ['ctags','no_gfortran']
    ignore = (d in non_source_dirs) or (d in useless_subdirs)
    return ignore

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

def contains_exactly_1_module( filepath, verbose=False ):
    """ Determine if only one module is present, otherwise ignore file.
    """
    num_mods = len( find_module_def( filepath ) )
    if len( find_program_name( filepath ) ) > 0:
        skip_file = True
        if verbose: print( "WARNING: Fortran program. ", end='' )
    elif num_mods > 1:
        skip_file = True
        if verbose: print( "WARNING: multiple module definitions. ", end='' )
    elif num_mods == 0:
        skip_file = True
        if verbose: print( "WARNING: no modules or programs in file. ", end='' )
    else:
        skip_file = False
    if skip_file:
        if verbose: print( "Skipping file '%s'" % filepath )
    return (not skip_file)

#==============================================================================
# MAIN FUNCTION
#==============================================================================

def create_interface_sections( root ):
    """
    Create interface sections for all modules (and programs?)
    in a Fortran library.

    Parameters
    ----------
    root : str
      Relative path to root of directory tree

    """
    from maintenance_tools import recursive_file_search
    from fortran_module    import FortranModule, populate_exported_symbols
    from fparser.api       import parse

    # Walk directory and store FortranModule objects
    modules = {}
    for fpath in recursive_file_search( root, ignore_dir, select_file ):
        if contains_exactly_1_module( fpath, verbose=True ):
            # Create fparser module object
            tree = parse( fpath, analyze=False )
            fmod = tree.content[0]
            # Create 'my' module object and store it in dictionary
            modules[fmod.name] = FortranModule( fpath, fmod )

    # [0] Link modules against used ones, creating a graph
    for name,mmod in modules.items():
        mmod.link_used_modules( *modules.values() )

    # [1] Update use statements (recursively search symbols in used modules)
    for name,mmod in modules.items():
        mmod.update_use_statements()

    #########################
    return        # STOP HERE
    #########################

    # [2] Populate exported symbols
    populate_exported_symbols( *modules.values() )

    # [3] Generate interface sections
    for name,mmod in modules.items():
        interface = mmod.generate_interface_section()
        filepath  = mmod.filepath[:-4] + '-interface.txt'
        with open( filepath, 'w' ) as f:
            print( interface, file=f )

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
      description = 'Create interface sections for all modules (and programs?)'
                    ' in a Fortran library.',
      epilog      = ' ',
      formatter_class = argparse.ArgumentDefaultsHelpFormatter,
      )

  parser.add_argument( metavar = 'ROOT',
                       dest    = 'root',
                       help    = 'relative path of the root directory' )

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

    # Walk directory tree and create interface sections
    create_interface_sections( args.root )

#------------------------------------------------------------------------------
if __name__ == '__main__':
    # Run as main program
    main()
