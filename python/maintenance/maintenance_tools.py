# coding: utf8
"""
Collection of low-level functions and classes, useful for performing standard
and exceptional maintenance operations in Selalib.

Selalib developers should not use directly this module, but rather refer to
various high-level specialized scripts.

Modules required
----------------
  * Built-in  : re, os

"""
#
# Author: Yaman Güçlü, Oct 2015 - IPP Garching
#
# Last revision: 13 Oct 2015
#

import re        # Regular expressions
import os        # Miscellaneous operating system operations

#==============================================================================
# FUNCTIONS
#==============================================================================

def iterate_over_file( filename, func, *args, **kwargs ):
    """ Iterate over file with a given function
    """
    with open( filename ) as f:
        for line in f:
            func( line, *args, **kwargs )


def recursive_file_search( path,
        ignore_dir  = lambda d: False,
        select_file = lambda f: True  ):
    """
    Iterate over a directory tree, ignoring directories and selecting files
    according to the given functions, and return relative file paths.

    Parameters
    ----------
    path : str
      Root directory to be searched

    ignore_dir : callable
      Function that returns True if the current directory should be ignored

    select_file: callable
      Function that returns True if the current file should be selected

    """
    for dirpath, dirnames, filenames in os.walk( path ):
        # Prune subdirectories that are not of interest
        for d in dirnames:
            if ignore_dir( d ):
                dirnames.remove( d )
        # Apply function to each file in directory
        for f in filenames:
            if select_file( f ):
                yield( os.path.join( dirpath, f ) )


def strip_comments( line ):
    """ Strip comments, trailing whitespaces and return character from line.
    """ 
    return line.rstrip( '\n' ).partition( '!' )[0].rstrip()


def is_preprocessor_line( line ):
    """ Return True if line is preprocessor instruction.
    """
    return line.lstrip().startswith('#')


def is_fortran_include( line ):
    """ Return True if line is Fortran include instruction.
    """
    return line.lstrip().startswith('include ')


def find_module_name( line, mod_def_list, mod_use_list ):
    """ Find module declaration or import, and add name to proper list.
    """
    spl = line.partition( '!' )[0].split( maxsplit=2 )
    if len( spl ) > 1:
        if spl[0] == 'module':
            next_word = spl[1]
            if next_word != 'procedure':
                mod_def_list.append( next_word )
        elif spl[0] == 'use':
            next_word = spl[1].partition( ',' )[0]
            mod_use_list.append( next_word )


def find_module_def( filename ):
    mod_def_list = []
    with open( filename, 'r' ) as f:
        for line in f:
            spl = line.partition( '!' )[0].split( maxsplit=2 )
            if len( spl ) > 1:
                if spl[0] == 'module' and spl[1] != 'procedure':
                    mod_def_list.append( spl[1] )
    return mod_def_list


def find_program_name( filename ):
    program_list = []
    with open( filename, 'r' ) as f:
        for line in f:
            spl = line.partition( '!' )[0].split( maxsplit=2 )
            if len( spl ) > 1:
                if spl[0] == 'program':
                    program_list.append( spl[1] )
    return program_list


def get_file_metadata( filename ):
    """ Determine various info about Fortran source file, and store it in
        dictionary.
    """
    # Extract module names from file
    module_def = []
    module_use = []
    iterate_over_file( filename,
            find_module_name, module_def, module_use )

    # Give error if there are multiple definitions of the same module!
    if not all_different( *module_def ):
        err_msg = "Error in file " + filename + \
                  ": multiple definitions of the same module"
        raise SystemExit( err_msg )

    # Remove multiple use statements
    if not all_different( *module_use ):
        module_use = sorted( set( module_use ) )

    # Return dictionary with all information
    md = dict()
    md['filename']   = filename
    md['module_def'] = module_def
    md['module_use'] = module_use
    return md


def compute_file_domain( path ):
    """
    Based on file and directory names, determine if a certain path belongs to
    one of the following domains in Selalib: library, testing or simulations.

    """
    dirname,filename = os.path.split( path )
    split_dir = dirname.split( '/' )
    if 'testing' in split_dir or \
            filename.startswith( 'test_' ) or \
            filename.startswith( 'unit_test_' ):
        return 'testing'
    if 'simulations' in split_dir:
        return 'simulations'
    return 'library'


def all_different( *items ):
    """ Return True if no equal items exist.
    """
    return len( items ) == len( set( items ) )


def get_repetition_count( *items ):
    """ Generator that provides (item,count) pairs in case item is repeated.
    """
    smallest_set = sorted( set( items ) )
    for item in smallest_set:
        count = items.count( item )
        if count > 1:
            yield (item,count)


def convert_module_name( name ):
    """ Dummy function that converts module names
    """
    return "NEW_" + name


def user_select( prompt, choices ):
    """ Ask user to choose one of the possible choices.
    """
    choices_str = [ str( c ) for c in choices]
    message = prompt + " [" + ",".join( choices_str ) + "]: "
    s = input( message )
    while s not in choices_str:
        s = input( "Cannot understand input, try again: " )
    return choices[choices_str.index( s )]


def find_utf8_error( filename ):
    """ Search for first non-utf8 character in file, and return line number.
        Returns 0 if no error is found.
    """
    with open( filename, 'r' ) as f:
        try:
            n = 1
            for line in f:
                n += 1
        except UnicodeDecodeError:
            return n
        else:
            return 0

#==============================================================================
# CLASSES
#==============================================================================

class RegexSubEngine( object ):
    """
    Engine for substituting 'old names' with 'new names' in a string.

    Parameters
    ----------
    oldnames: list
      Names to be changed (matched as full words)

    newnames: list
      New names that replace the old ones

    """
    __slots__ = ['_oldnames', '_newnames', '_re_obj']

    def __init__( self, oldnames, newnames, lookbehind='', lookahead='' ):
       pattern_code = lookbehind + r'\b%s\b' + lookahead
       pattern = '|'.join( pattern_code % re.escape( w ) for w in oldnames )
       re_obj  = re.compile( pattern )
       self._oldnames = oldnames 
       self._newnames = newnames
       self._re_obj   = re_obj

    def _replace_func( self, match ):
        w = match.group( 0 )
        return self._newnames[ self._oldnames.index( w ) ]

    def process( self, string ):
        """ Apply all substitutions to given string, and return new string.
        """
        return re.sub( self._re_obj, self._replace_func, string )

    @property
    def oldnames( self ): return self._oldnames

    @property
    def newnames( self ): return self._newnames

#==============================================================================
