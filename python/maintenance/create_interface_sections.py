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
# Last revision: 18 Nov 2015
#
from __future__ import print_function

__all__ = ['create_interface_sections','main']
__docformat__ = 'reStructuredText'

#==============================================================================
# PARAMETERS.. CHANGE THIS
#==============================================================================

ignored_symbols = [ \
        'mudpack_curvilinear_cof',
        'mudpack_curvilinear_cofcr',
        'mudpack_curvilinear_bndcr',
        'sol',
        'mulku',
        ]

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

def create_interface_sections( root, src='src', interfaces='src/interfaces' ):
    """
    Create interface sections for all modules (and programs?)
    in a Fortran library.

    Parameters
    ----------
    root : str
      Relative path to root of directory tree

    """
    import os
    from maintenance_tools import recursive_file_search
    from fortran_module    import FortranModule
    from fortran_external  import external_modules, find_external_library
    from fparser.api       import parse

    # Walk library tree and store FortranModule objects
    print( "================================================================" )
    print( "Processing library modules")
    print( "================================================================" )
    src_modules = {}
    src_root    = os.path.join( root, src )
    for fpath in recursive_file_search( src_root, ignore_dir, select_file ):
        if interfaces in fpath:
            print( "WARNING: Interface. Skipping file %s" % fpath )
            continue
        if contains_exactly_1_module( fpath, verbose=True ):
            # Create fparser module object
            tree = parse( fpath, analyze=False )
            fmod = tree.content[0]
            # Create 'my' module object and store it in dictionary
            src_modules[fmod.name] = FortranModule( fpath, fmod )

    # Additional source directory (ad-hoc)
    src_root = os.path.join( root, 'external/burkardt' )
    for fpath in recursive_file_search( src_root, ignore_dir, select_file ):
        if contains_exactly_1_module( fpath, verbose=True ):
            # Create fparser module object
            tree = parse( fpath, analyze=False )
            fmod = tree.content[0]
            # Create 'my' module object and store it in dictionary
            src_modules[fmod.name] = FortranModule( fpath, fmod )

    # Interface modules
    print( "================================================================" )
    print( "Processing interface modules" )
    print( "================================================================" )
    int_modules = {}
    int_root    = os.path.join( root, interfaces )
    for fpath in recursive_file_search( int_root, ignore_dir, select_file ):
        if contains_exactly_1_module( fpath, verbose=True ):
            # Create fparser module object
            tree = parse( fpath, analyze=False )
            fmod = tree.content[0]
            # Create 'my' module object and store it in dictionary
            int_modules[fmod.name] = FortranModule( fpath, fmod )

    # Dictionary with all the modules
    all_modules = {}
    all_modules.update( src_modules )
    all_modules.update( int_modules )

    print( "Source modules:" )
    for src_mod in src_modules.values():
        print( src_mod.name )
    print()

    print( "Interface modules:" )
    for int_mod in int_modules.values():
        print( int_mod.name )
    print()

    print( "================================================================" )
    print( "Link library modules against used ones" )
    print( "================================================================" )
    # [0] Link modules against used ones, creating a graph
    for i,(name,mmod) in enumerate( src_modules.items() ):
        print("  - link module %3d: %s" % (i+1,name) )
        mmod.link_used_modules( all_modules.values(), externals=external_modules )

    print( "================================================================" )
    print( "Search symbols in used modules" )
    print( "================================================================" )
    # [1] Update use statements (recursively search symbols in used modules)
    for i,(name,mmod) in enumerate( src_modules.items() ):
        print("  - update module %3d: %s" % (i+1,name) )
        mmod.update_use_statements( find_external_library, ignored_symbols )

    print( "================================================================" )
    print( "Cleanup imported symbols lists" )
    print( "================================================================" )
    # [2] Cleanup use statements (remove duplicate symbols and useless modules)
    for i,(name,mmod) in enumerate( src_modules.items() ):
        print("  - cleanup module %3d: %s" % (i+1,name) )
        mmod.cleanup_use_statements()

    print( "================================================================" )
    print( "Scatter imported symbols" )
    print( "================================================================" )
    # [3] Populate exported symbols
    for i,(name,mmod) in enumerate( src_modules.items() ):
        print("  - scatter from module %3d: %s" % (i+1,name) )
        mmod.scatter_exported_symbols()

    print( "================================================================" )
    print( "Generate interface sections" )
    print( "================================================================" )
    # [4] Generate interface sections
    for i,(name,mmod) in enumerate( src_modules.items() ):
        print("  - interface for module %3d: %s" % (i+1,name) )
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

