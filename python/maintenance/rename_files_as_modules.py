# coding: utf8
"""
Recursively search directory tree for Fortran module definitions,
change file name to match Fortran module, and make appropriate
substitutions within all 'CMakeLists.txt' files.

Modules required
----------------
  * Built-in  : os, sys, argparse
  * Library   : maintenance_tools

"""
#
# Author: Yaman Güçlü, Oct 2015 - IPP Garching
#
# Last revision: 15 Oct 2015
#

__all__ = ['match_files_to_modules','main']
__docformat__ = 'reStructuredText'

#==============================================================================
# PARAMETERS.. CHANGE THIS
#==============================================================================

def ignore_dir( d ):
    """ Return True if subdirectory should be ignored.
    """
    return d in ['ctags', 'no_gfortran', 'testing', 'simulations']


def is_source( f ):
    """ Return True if filename should be selected for processing.
    """
    is_fortran = f.endswith( '.F90' )
    is_doc     = f.endswith( '_doc.F90' )
    is_test    = f.startswith( 'test' ) or f.startswith( 'unit_test' )
    return is_fortran and (not is_doc) and (not is_test)

#==============================================================================
# MAIN FUNCTION
#==============================================================================

def rename_files_as_modules( root, dry_run=True, additional=None ):
    """
    Recursively search directory tree for Fortran module definitions,
    change file name to match Fortran module, and make appropriate
    substitutions within all 'CMakeLists.txt' files.

    Parameters
    ----------
    root : str
      Relative path to root of directory tree

    dry_run : bool
      If True, do not change Fortran file names and do not modify cmake files

    additional : str
      Relative path to additional file where substitutions are needed

    """
    import os, shutil
    import maintenance_tools as mtools

    # If needed, create temporary for additional file
    if additional is not None:
        if dry_run:
            tmp1 = additional + '.new'
            shutil.copy( additional, tmp1 )
        else:
            tmp1 = additional

    for dirpath, dirnames, filenames in os.walk( root ):
        # Give info
        print( "Walking directory '%s'" % dirpath )

        # Prune subdirectories that are not of interest
        for d in dirnames:
            if ignore_dir( d ):
                dirnames.remove( d )

        # If there is no cmake file, skip current directory!
        if not 'CMakeLists.txt' in filenames:
            print( "  WARNING: no 'CMakeLists.txt' file, skipping directory" )
            print()
            continue

        # Select all Fortran source files, excluding doc files and unit tests
        source_filenames = sorted( f for f in filenames if is_source( f ) )

        # Select files containing EXACTLY ONE module and NO programs
        for f in list( source_filenames ):
            filepath  = os.path.join( dirpath, f )
            num_mods  = len( mtools.find_module_def( filepath ) )

            if len( mtools.find_program_name( filepath ) ) > 0:
                skip_file = True
                print( "  WARNING: Fortran program. ", end='' )
            elif num_mods > 1:
                skip_file = True
                print( "  WARNING: multiple module definitions. ", end='' )
            elif num_mods == 0:
                skip_file = True
                print( "  WARNING: no modules or programs in file. ", end='' )
            else:
                skip_file = False

            if skip_file:
                source_filenames.remove( f )
                print( "Skipping file '%s'" % f )

        # Collect Fortran module names
        modules = [mtools.find_module_def( os.path.join( dirpath, f ) )[0] \
                   for f in source_filenames]

        # Remove repeated modules from list, and skip files that contain them
        if not mtools.all_different( *modules ):
            print( "  WARNING: repeated module names in directory:" )
            for (name,count) in mtools.get_repetition_count( *modules ):
                print( "           %d occurrencies for '%s'" % (count,name) )
                src_skip = \
                   [s for m,s in zip( modules, source_filenames ) if m == name]
                for s in src_skip:
                    modules.remove( name )
                    source_filenames.remove( s )
                    print("             . skipping file '%s'" % s )

        # Sanity check
        assert( len( modules ) == len( source_filenames ) )

        # Create list of (old_name,new_name) substitution pairs
        filename_subs = []
        for m,f in zip( modules, source_filenames ):
            fnew = m + '.F90'
            if fnew == f:
                print( "  File is OK: %s" % f )
            else:
                filename_subs.append( (f,fnew) )

        # If no substitutions are needed, move to next directory
        if len( filename_subs ) == 0:
            print( "  Nothing to do here", end='\n\n' )
            continue

        # Change file names
        maxlen   = max( [len( f ) for f,fnew in filename_subs] )
        template = "  Renaming  {{:<{:d}s}}  >>  {{:s}}".format( maxlen )
        for old_filename, new_filename in filename_subs:
            print( template.format( old_filename, new_filename ) )
            if not dry_run:
                old_filepath = os.path.join( dirpath, old_filename )
                new_filepath = os.path.join( dirpath, new_filename )
                os.rename( old_filepath, new_filepath )

        # Create greedy algorithm to perform name substitutions
        engine = mtools.RegexSubEngine( *zip( *filename_subs ) )

        # Substitute names in all CMakeLists.txt files
        print( "  Processing 'CMakeLists.txt'" )
        old_filepath = os.path.join( dirpath, 'CMakeLists.txt' )
        new_filepath = old_filepath + '.new'
        with open( old_filepath, 'r' ) as f1, \
             open( new_filepath, 'w' ) as f2:
            for old_line in f1:
                new_line = engine.process( old_line )
                print( new_line, file=f2, end='' )
        if not dry_run:
            os.rename( new_filepath, old_filepath )

        # Update names in additional file, if given
        if additional is not None:
            print( "  Processing additional file: " + additional )
            tmp2 = tmp1 + '.new'
            with open( tmp1, 'r' ) as f1, \
                 open( tmp2, 'w' ) as f2:
                for old_line in f1:
                    new_line = engine.process( old_line )
                    print( new_line, file=f2, end='' )
            os.rename( tmp2, tmp1 )

        # Print empty line
        print()

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
      description = 'Change file names in a Fortran library, '
                    'in order to match the module names',
      epilog      = ' ',
      formatter_class = argparse.ArgumentDefaultsHelpFormatter,
      )

  parser.add_argument( metavar = 'ROOT',
                       dest    = 'root',
                       help    = 'relative path of the root directory' )

  parser.add_argument( '-d', '--dry_run',
                       action  = 'store_true',
                       help    = 'do not substitute old files with new ones' )

  parser.add_argument( '-a', '--additional',
                       metavar = 'FILE',
                       default = None,
                       help    = 'change names in additional file' )

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
    rename_files_as_modules( args.root, args.dry_run, args.additional )

#------------------------------------------------------------------------------
if __name__ == '__main__':
    # Run as main program
    main()
