# coding: utf8
"""
Post-process interface sections:
    1. change some use statements with header include
    2. add some preprocessor #ifdef
    3. include additional headers

Modules required
----------------
  * Built-in  : re, os, sys, argparse
  * Library   : maintenance_tools

"""
#
# Author: Yaman Güçlü, Nov 2015 - IPP Garching
#
# Last revision: 27 Nov 2015
#
import re

__all__ = ['main']
__docformat__ = 'reStructuredText'

#==============================================================================
# PARAMETERS.. CHANGE THIS
#==============================================================================

module_to_header_subs = { \
'sll_m_assert'            : 'sll_assert.h',
'sll_m_errors'            : 'sll_errors.h',
'sll_m_memory'            : 'sll_memory.h',
'sll_m_working_precision' : 'sll_working_precision.h',
}

ignored_headers = { \
'sll_assert.h',
'sll_errors.h',
'sll_memory.h',
'sll_working_precision.h',
'sll_fft.h',
#'sll_fftw.h',
}

interface_ending = '-interface.txt' 

def is_interface_file( fpath ):
    return fpath.endswith( interface_ending )

use_pattern = r'\s*use\s+([a-zA-Z]\w*)'
use_regex   = re.compile( use_pattern, re.I )

include_pattern = r'#include\s+"([a-zA-Z]\w*\.h)"'
include_regex   = re.compile( include_pattern, re.I )

ifdef_pattern = r'#ifn?def\s+([\w^\d]\w*)'
ifdef_regex   = re.compile( ifdef_pattern, re.I )

#==============================================================================
# UTILITY FUNCTIONS
#==============================================================================

def get_original_file( fpath, original_root, preproc_root ):
    import os.path

    n = len( interface_ending )
    dname, fname = os.path.split( fpath )
    dname_orig = dname.replace( preproc_root, original_root )
    fname_orig = fname[:-n] + '.F90'

    # If original file exists, return filepath
    fpath_orig = os.path.join( dname_orig, fname_orig )
    if os.path.isfile( fpath_orig ):
        return fpath_orig

    # Otherwise, search file with '.f90' extension
    fpath_orig = os.path.join( dname_orig, fname_orig[:-4] + '.f90' )
    if os.path.isfile( fpath_orig ):
        return fpath_orig

    # Otherwise, search in parent directory
    fpath_orig = os.path.join( os.path.dirname( dname_orig ), fname_orig )
    if os.path.isfile( fpath_orig ):
        return fpath_orig

    # Otherwise, raise error
    print( "ERROR: original file not found" )
    print( "  interface file: %s" % fpath )
    raise SystemExit()

#------------------------------------------------------------------------------
def convert_modules_to_headers( f ):
    headers  = []
    lines    = []
    add_line = True
    for line in f:
        if add_line:
            match = use_regex.match( line )
            if match:
                mod_name = match.group( 1 )
                if mod_name in module_to_header_subs.keys():
                    add_line = False
                    headers.append( module_to_header_subs[mod_name] )
            if add_line:
                lines.append( line.rstrip( '\n' ) )
        else:
            # If line is empty, next line should be added in general
            if not line.strip():
                add_line = True

    return headers, lines

#------------------------------------------------------------------------------
def extract_use_blocks( f, *module_names ):
    blocks  = {}
    lines   = []
    extract = False
    for line in f:
        line = line.rstrip('\n')
        if extract:
            if not line.strip():
                extract = False
            else:
                block.append( line )
        else:
            match = use_regex.match( line )
            if match:
                mod_name = match.group( 1 )
                if mod_name in module_names:
                    block = [line]
                    blocks[mod_name] = block
                    extract = True
            if not extract:
                lines.append( line )

    return lines, blocks

#------------------------------------------------------------------------------
def expand_use_statements( f, **blocks ):
    new_lines = []
    for line in f:
        match = use_regex.match( line )
        mod_name = match.group( 1 ) if match else None
        if mod_name in blocks.keys():
            new_lines.extend( blocks[mod_name] )
            new_lines.append( '' )
        else:
            new_lines.append( line )
    return new_lines

#------------------------------------------------------------------------------
def find_missing_headers( f ):
    missing_headers = []
    for line in f:
        match = include_regex.match( line )
        if match:
            header = match.group( 1 )
            if header not in ignored_headers:
                missing_headers.append( header )
    return missing_headers

#------------------------------------------------------------------------------
def extract_interface_section( f ):
    # Identify the meaningful lines to be processed
    lines = []
    skip  = True
    for i,line in enumerate( f ):
        line = line.lstrip().rstrip('\n')
        if line.startswith( 'module' ) or line.startswith( 'program' ):
            skip = False

        elif line.startswith( 'implicit' ) or line.startswith( 'contains' ):
            break
        if not skip:
            lines.append( line )
    return lines

#------------------------------------------------------------------------------
def find_ifdef_clauses( lines ):
    # Extract #ifdef and #ifndef clauses from the lines above
    clauses  = []
    outer_if = None
    inclause = False
    modules  = set()
    for line in lines:
        if line.startswith('#if'):
            inclause   = True
            new_clause = []
            if not ifdef_regex.match( line ):
                print( "ERROR: not a #ifdef nor a #ifndef clause")
                raise SystemExit()

        if inclause:
            if line.startswith('#'):
                new_clause.append( line )
            else:
                match = use_regex.match( line )
                if match:
                    new_clause.append( match.group( 0 ) )
                    modules.add( match.group( 1 ) )

        if line.startswith('#endif'):
            inclause = False
            clauses.append( new_clause )

    # Last clause is open and requires special treatment
    if inclause:
        if clauses:
            # If this is not the only clause, return an open clause
            clauses.append( new_clause )
        else:
            # If this is the only clause, treat as an outer if
            outer_if = new_clause[0]

    # Remove empty clauses
    for c in list( clauses ):
        num_use = sum( bool( use_regex.match( line ) ) for line in c )
        num_pp  = sum( line.startswith('#') for line in c if not \
          (line.startswith('#ifdef') or line.startswith('#ifndef') \
                                     or line.startswith('#else') ) )
        if (num_use + num_pp) == 0:
            clauses.remove( c )
#            print("REMOVE")

    # Return outer-if (if any) and list of clauses
    return outer_if, clauses

#------------------------------------------------------------------------------
def find_all_used_modules( lines ):
    modules = []
    for line in lines:
        match = use_regex.match( line )
        if match:
            modules.append( match.group( 1 ) )
    return modules

#------------------------------------------------------------------------------
#def add_openmp_clause( f ):
#    lines = []
#    ifdef_block = None
#    for line in f:
#        match = use_regex.match( line )
#        if match:
#            mod_name = match.group( 1 )
#            if mod_name == 'omp_lib':
#                lines.append( "#ifdef _OPENMP" )
#                ifdef_block = "_OPENMP"
#                print()
#                print( "OPENMP LIB" )
#                print()
#        if ifdef_block and not line.strip():
#            lines.append( "#endif // %s" % ifdef_block )
#            ifdef_block = None
#        lines.append( line )
#    return lines

#==============================================================================
# MAIN FUNCTION
#==============================================================================

def process_interface_sections( original_root, preproc_root ):

    from maintenance_tools import recursive_file_search

    original_root = original_root.rstrip('/')
    preproc_root  =  preproc_root.rstrip('/')

#    # Counters
#    num_ifclauses  = 0
#    num_ifmods     = 0

    for fpath in recursive_file_search( 
            preproc_root, select_file = is_interface_file ):

        fpath_orig = get_original_file( fpath, original_root, preproc_root )
#        fpath_prpr = fpath.rpartition( '-interface.txt' )[0] + '.f90'
        print()
        print( "Post-processing interface: %s" % fpath )
#        print( "     . pre-processed file: %s" % fpath_prpr )
        print( "   . original source file: %s" % fpath_orig )

        #----------------------------------------------------------------------
        # Convert some specific modules to headers
        with open( fpath, 'r' ) as f:
            headers, lines = convert_modules_to_headers( f )

        # Collect all headers that should not be ignored
        with open( fpath_orig, 'r' ) as f_orig:
            headers += find_missing_headers( f_orig )

        # Add header lines at beginning of file
        if headers:
            lines = ['#include "%s"' % h for h in headers] + [''] + lines

        # Find ifdef clauses and contained use statements
        with open( fpath_orig, 'r' ) as f_orig:
            interface_section_orig = extract_interface_section( f_orig )
        outer_if, ifclauses = find_ifdef_clauses( interface_section_orig )
        ifmods = set(m for c in ifclauses for m in find_all_used_modules( c ))

        # Add outer ifdef at very beginning of file
        if outer_if:
            lines = [outer_if,''] + lines

#        #----------------------------------------------------------------------
#        # TEST
#        for c in ifclauses:
#            print()
#            print( '\n'.join( c ) )
#            print()
#        if ifmods:
#            print( ifmods )
#            print()
#        # Update counters
#        num_ifclauses  += len( ifclauses  )
#        num_ifmods     += len( ifmods )
#        #----------------------------------------------------------------------

        # Extract use statements that belong to ifdef clauses
        lines, blocks = extract_use_blocks( lines, *ifmods )

        # Expand ifdef clauses using the 'blocks' dictionary of use statements
        ifclauses = [expand_use_statements( c, **blocks ) for c in ifclauses]

        # Concatenate all lines from ifdef clauses, separated by a blank line
        iflines = []
        for c in ifclauses:
            iflines.extend( c )
            iflines.append('')
        if iflines:
            iflines.pop()

        # Insert all ifdef lines before the 'implicit none' statement
        for i,l in enumerate( lines ):
            if l.lstrip().startswith('implicit'):
                break
        lines = lines[:i] + iflines + lines[i:]

        # Save interface to a new file
        new_fpath = fpath[:-4] + '2.txt'
        with open( new_fpath, 'w' ) as new_f:
            print( '\n'.join( lines ), file=new_f )

        #----------------------------------------------------------------------

#    # Some statistics
#    print( "No. of #ifdef/#ifndef clauses = %d" % num_ifclauses  )
#    print( "No. of effected modules       = %d" % num_ifmods )

#==============================================================================
# PARSER
#==============================================================================

def parse_input():

  import argparse, sys

  parser = argparse.ArgumentParser (
      prog        = 'python3 '+ sys.argv[0],
      description = 'Post-process generated interface sections.',
      epilog      = ' ',
      formatter_class = argparse.ArgumentDefaultsHelpFormatter,
      )

  parser.add_argument(
          metavar = 'ROOT-ORIGINAL',
          dest    = 'original_root',
          help    = 'relative path of root directory with library files' )

  parser.add_argument(
          metavar = 'ROOT-PREPROC',
          dest    = 'preproc_root',
          help    = 'relative path of root directory with preprocessed files' )

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

    # Walk directory tree and post-process all interface sections
    process_interface_sections( args.original_root, args.preproc_root )

#------------------------------------------------------------------------------
if __name__ == '__main__':
    # Run as main program
    main()

