# coding: utf8
"""
Create interface sections for all modules (and programs?) in a Fortran library.

Modules required
----------------
  * Built-in  : os, sys, argparse
  * Library   : maintenance_tools, sll2py

"""
#
# Author: Yaman Güçlü, Oct 2015 - IPP Garching
#
# Last revision: 03 Dec 2015
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
# Overwrite instance methods (special cases)
#==============================================================================

def add_exported_symbols_permissive( self, *symbols ):
    for s in symbols:
        if not self.defines_symbol( s ):
            origin = self.find_symbol_def( s )
            if origin:
                print( "WARNING processing file '%s':" % self.filepath )
                print( "  symbol '%s' is imported but not defined here" % s )
                print( "  original definition in module '%s'" % origin )
            else:
                print( "ERROR processing file '%s':" % self.filepath )
                print( "  symbol '%' is neither defined nor imported here" % s )
                raise SystemExit()
        self._exported_symbols.add( s )

permissive_modules = ['sll_m_hdf5_io_parallel']

def make_modules_permissive( *modules ):
    from types import MethodType
    for m in modules:
        if m.name in permissive_modules:
            m.add_exported_symbols = \
                    MethodType( add_exported_symbols_permissive, m )

#==============================================================================
# HELPER FUNCTION
#==============================================================================

def parse_file_and_create_unit( fpath, module_dict, program_dict ):

    from sll2py.fparser.api      import parse
    from sll2py.fortran_module   import NewFortranModule
    from sll2py.fortran_program  import FortranProgram

    im = len(  module_dict )
    ip = len( program_dict )

    mod_names = find_module_def  ( fpath )
    prg_names = find_program_name( fpath )

    if len( mod_names ) == 1 and len( prg_names ) == 0:
        # Create fparser module object
        tree = parse( fpath, analyze=False )
        fmod = tree.content[0]
        # Create 'my' module object and store it in dictionary
        print("  - read module  %3d: %s" % (im+1, fmod.name ) );  im += 1
        module_dict[fmod.name] = NewFortranModule( fpath, fmod )
    elif len( mod_names ) == 0 and len( prg_names ) == 1:
        # Create fparser program object
        tree = parse( fpath, analyze=False )
        fprg = tree.content[0]
        # Create 'my' program object and store it in dictionary
        print("  - read program %3d: %s" % (ip+1, fprg.name ) );  ip += 1
        program_dict[fprg.name] = FortranProgram( fpath, fprg )
    elif len( mod_names ) == 0 and len( prg_names ) == 0:
        print( "ERROR: No modules or programs in file %s" % fpath )
        raise SystemExit()
    else:
        print( "ERROR: Multiple modules/programs in file %s" % fpath )
        if mod_names: print( "  Modules  = [%s]" % ', '.join( mod_names ) )
        if prg_names: print( "  Programs = [%s]" % ', '.join( prg_names ) )
        raise SystemExit()

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
    from maintenance_tools        import recursive_file_search
    from sll2py.fortran_external  import external_modules, find_external_library

    print( "================================================================" )
    print( "[1] Processing library modules and programs")
    print( "================================================================" )
    src_modules  = {}
    src_programs = {}

    # Walk library tree and store FortranModule (or FortranProgram) objects
    src_root = os.path.join( root, src )
    for fpath in recursive_file_search( src_root, ignore_dir, select_file ):
        if interfaces in fpath:
            print( "WARNING: Interface. Skipping file %s" % fpath )
            continue
        parse_file_and_create_unit( fpath, src_modules, src_programs )

    # Additional source directory (ad-hoc)
    src_root = os.path.join( root, 'external/burkardt' )
    for fpath in recursive_file_search( src_root, ignore_dir, select_file ):
        parse_file_and_create_unit( fpath, src_modules, src_programs )

    # Interface modules
    print( "================================================================" )
    print( "[2] Processing interface modules and programs" )
    print( "================================================================" )
    int_modules  = {}
    int_programs = {}
    int_root     = os.path.join( root, interfaces )
    for fpath in recursive_file_search( int_root, ignore_dir, select_file ):
        parse_file_and_create_unit( fpath, int_modules, int_programs )

    # Dictionary with all the modules
    all_modules = {}
    all_modules.update( src_modules )
    all_modules.update( int_modules )

    print( "================================================================" )
    print( "[3] Library modules/programs: Link against used modules" )
    print( "================================================================" )
    # Link library modules/programs against used modules, creating a graph
    for i,(name,mmod) in enumerate( src_modules.items() ):
        print("  - link module %3d: %s" % (i+1,name) )
        mmod.link_used_modules( all_modules.values(), externals=external_modules )

    # TODO: what about interface programs?
    print( "----------------------------------------------------------------" )
    for i,(name,mprg) in enumerate( src_programs.items() ):
        print("  - link program %3d: %s" % (i+1,name) )
        mprg.link_used_modules( all_modules.values(), externals=external_modules )

    print( "================================================================" )
    print( "[4] Library modules/programs: Search symbols in used modules" )
    print( "================================================================" )
    # Update use statements (recursively search symbols in used modules)
    for i,(name,mmod) in enumerate( src_modules.items() ):
        print("  - update module %3d: %s" % (i+1,name) )
        mmod.update_use_statements( find_external_library, ignored_symbols )

    print( "----------------------------------------------------------------" )
    for i,(name,mprg) in enumerate( src_programs.items() ):
        print("  - update program %3d: %s" % (i+1,name) )
        mprg.update_use_statements( find_external_library, ignored_symbols )

    print( "================================================================" )
    print( "[5] Library modules/programs: Cleanup use statements" )
    print( "================================================================" )
    # Cleanup use statements (remove duplicate symbols and useless modules)
    for i,(name,mmod) in enumerate( src_modules.items() ):
        print("  - cleanup module %3d: %s" % (i+1,name) )
        mmod.cleanup_use_statements()

    print( "----------------------------------------------------------------" )
    for i,(name,mprg) in enumerate( src_programs.items() ):
        print("  - cleanup program %3d: %s" % (i+1,name) )
        mprg.cleanup_use_statements()

    print( "================================================================" )
    print( "[6] Library modules/programs: Scatter imported symbols" )
    print( "================================================================" )

    # Some modules must be made permissive
    make_modules_permissive( *src_modules.values() )

    for i,(name,mmod) in enumerate( src_modules.items() ):
        print("  - scatter from module %3d: %s" % (i+1,name) )
        mmod.scatter_imported_symbols()

    print( "----------------------------------------------------------------" )
    for i,(name,mprg) in enumerate( src_programs.items() ):
        print("  - scatter from program %3d: %s" % (i+1,name) )
        mprg.scatter_imported_symbols()

    print( "================================================================" )
    print( "[7] Library modules/programs: Generate interface sections" )
    print( "================================================================" )
    for i,(name,mmod) in enumerate( src_modules.items() ):
        print("  - interface for module %3d: %s" % (i+1,name) )
        interface = mmod.generate_interface_section()
        if interface:
            filepath  = mmod.filepath[:-4] + '-interface.txt'
            with open( filepath, 'w' ) as f:
                print( interface, file=f )

    print( "----------------------------------------------------------------" )
    for i,(name,mprg) in enumerate( src_programs.items() ):
        print("  - interface for program %3d: %s" % (i+1,name) )
        interface = mprg.generate_interface_section()
        filepath  = mprg.filepath[:-4] + '-interface.txt'
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

