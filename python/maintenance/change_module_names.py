# coding: utf8
"""
Recursively search directory tree for Fortran module definitions, and perform
name substitutions according to an external function.

Modules required
----------------
  * Built-in  : os, sys, argparse, itertools
  * Library   : maintenance_tools, renaming

"""
#
# Author: Yaman Güçlü, Oct 2015 - IPP Garching
#
# Last revision: 15 Oct 2015
#

__all__ = ['test_utf8_errors','change_module_names','main']
__docformat__ = 'reStructuredText'

#==============================================================================
# PARAMETERS.. CHANGE THIS
#==============================================================================

def ignore_dir( d ):
    """ Return True if subdirectory should be ignored.
    """
    return d in ['ctags', 'no_gfortran']


def select_file( f ):
    """ Return True if a file should be selected for processing.
    """
    is_fortran = f.endswith( '.F90' )
    is_header  = f.endswith( '.h' )
    return is_fortran or is_header

# File domains where module definitions should be changed
new_module_domains = ['library','simulations']

#==============================================================================
# SIMPLE TESTS
#==============================================================================

def test():
    """ Test: walk directory and print filenames.
    """
    from maintenance_tools import recursive_file_search
    directory = 'tests/selalib/src/meshes'
    print( 'Walking on directory ' + directory + ':' )
    for f in recursive_file_search( directory, ignore_dir, select_file ):
        print( '  . ' + f )
    print( 'DONE' )


def test_utf8_errors( root ):
    """ Recursively search 'root' directory tree for UTF8 errors in files.
    """
    import maintenance_tools as mtools
    print( "Function 'test_utf8_errors': walking over directory '%s'.."% root )
    for f in mtools.recursive_file_search( root, ignore_dir, select_file ):
        n = mtools.find_utf8_error( f )
        if n > 0: print( "  . Error in line %d of file %s" % (n,f) )
    print( 'DONE' )

#==============================================================================
# MAIN FUNCTION
#==============================================================================

def change_module_names( root, dry_run=True ):
    """
    Recursively search directory tree for Fortran module definitions,
    and perform name substitutions according to an external function.

    Parameters
    ----------
    root : str
      Relative path to root of directory tree

    dry_run : bool
      If True, do not substitute names in the original files

    """
    import os, itertools
    import maintenance_tools as mtools
    from renaming import renaming as convert_module_name

    # Walk directory and store fileinfo
    file_metadata_list = []
    for f in mtools.recursive_file_search( root, ignore_dir, select_file ):
        md = mtools.get_file_metadata( f )
        file_metadata_list.append( md )

    # Iterate over all metadata, and determine file type (library, testing, sims)
    for md in file_metadata_list:
        md['domain'] = mtools.compute_file_domain( md['filename'] )

    # Collect all module definitions from library domain
    module_def_list = [md['module_def'] for md in file_metadata_list \
                                         if md['domain'] in new_module_domains]
    module_def_list = sorted( itertools.chain( *module_def_list ) )

    # Verify if there are repeated names
    if not mtools.all_different( *module_def_list ):
        print( "\nWarning: repeated old module names:" )
        for (name,count) in mtools.get_repetition_count( *module_def_list ):
            print( "  %d occurrencies for %s:" % (count,name) )
            for md in file_metadata_list:
                if name in md['module_def']:
                    print( "    - " + md['filename'] )
        stop = mtools.user_select( "Do you wish to interrupt execution?", \
                                   ['Y','N'] )
        if stop == 'Y': raise SystemExit()

    # Convert all names of the defined modules to new values
    new_module_names = [convert_module_name( w ) for w in module_def_list ]

    # Iterate over all metadata, and add new names to file dictionary
    for md in file_metadata_list:
        if md['domain'] in new_module_domains:
            md['new_module_def'] = []
            for old_name in md['module_def']:
                idx = module_def_list.index( old_name )
                new_name = new_module_names[idx]
                md['new_module_def'].append( new_name )
        else:
            md['new_module_def'] = list( md['module_def'] )

    # Again, verify that there are no repeated names
    if not mtools.all_different( *new_module_names ):
        print( "\nWarning: repeated new module names:" )
        for (name,count) in mtools.get_repetition_count( *new_module_names ):
            print( "  %d occurrencies for %s:" % (count,name) )
            for md in file_metadata_list:
                if name in md['new_module_def']:
                    print( "    - " + md['filename'] )
        stop = mtools.user_select( "Do you wish to interrupt execution?", \
                                   ['Y','N'] )
        if stop == 'Y': raise SystemExit()

    # Cycle through all files, and substitute all names (greedy algorithm)
    greedy_engine = mtools.RegexSubEngine( module_def_list, new_module_names )
    for md in file_metadata_list:
        print( "Processing " + md['filename'] )
        old_file  = md['filename']
        new_file  = md['filename'] + '.new_names'
        with open( old_file, 'r' ) as f1,\
             open( new_file, 'w' ) as f2:
            for old_line in f1:
                if mtools.is_preprocessor_line( old_line ) or \
                   mtools.is_fortran_include  ( old_line ):
                    new_line = old_line
                else:
                    new_line = greedy_engine.process( old_line )
                print( new_line, file=f2, end='' )

        # Substitute old file with new file
        if not dry_run:
            os.rename( new_file, old_file )

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
      description = 'Change all module names in a Fortran library, according '
                    'to a predefined algorithm.',
      epilog      = ' ',
      formatter_class = argparse.ArgumentDefaultsHelpFormatter,
      )

  parser.add_argument( metavar = 'ROOT',
                       dest    = 'root',
                       help    = 'relative path of the root directory' )

  parser.add_argument( '-t', '--test_utf8',
                       action  = 'store_true',
                       help    = 'do not run, check for utf8 errors instead' )

  parser.add_argument( '-d', '--dry_run',
                       action  = 'store_true',
                       help    = 'do not substitute old files with new ones' )

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

    if args.test_utf8:
        # Walk directory tree and only search for UTF8 errors
        test_utf8_errors( args.root )
    else:
        # Walk directory tree and change names of all library modules
        change_module_names( args.root, args.dry_run )

#------------------------------------------------------------------------------
if __name__ == '__main__':
    # Run as main program
    main()
