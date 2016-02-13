# coding: utf8
"""
Simple Python script for generating a list of external symbols for the MPI
library.

Modules required
----------------
  * Built-in  : sys, argparse, re

"""
#
# Author: Yaman Güçlü, Dec 2015 - IPP Garching
#
# Last revision: 17 Dec 2015
#
from __future__ import print_function
import re

__all__ = ['create_external_list','main']
__docformat__ = 'reStructuredText'

#===============================================================================
# Parameters
#===============================================================================

tab = 2*' '

#===============================================================================
# MAIN FUNCTION
#===============================================================================

def create_external_list( err_file, ext_file ):

    err_text = open( err_file, 'r' ).read()

    if err_text:
        pattern = r'\s+(mpi_\w+)\b'
        missing = [w.lower() for w in re.findall( pattern, err_text )]
        missing = sorted( set( missing ) )
    else:
        return

    if missing:
        max_wlen = max( len( w ) for w in missing )
        line_len = len( tab ) + max_wlen + 3
        lines = ['external :: %s&' % ((line_len-13)*' ')]
        for w in missing:
            lines.append( '%s%s, %s&' % (tab, w, (max_wlen-len( w ))*' ') )
        lines[-1] = lines[-1].rstrip( ', &' )
        with open( ext_file, 'w' ) as f:
            print( *lines, file=f, sep='\n' )

#===============================================================================
# PARSER
#===============================================================================

def parse_input():

  import argparse, sys

  parser = argparse.ArgumentParser (
      prog        = 'python '+ sys.argv[0],
      description = 'Create list of external symbols from the compiler errors.',
      epilog      = ' ',
      formatter_class = argparse.ArgumentDefaultsHelpFormatter,
      )

  parser.add_argument(
          metavar = 'ERR_FILE',
          dest    = 'err_file',
          help    = 'standard error from compiler' )

  parser.add_argument(
          metavar = 'EXT_FILE',
          dest    = 'ext_file',
          help    = 'file to be generated, with all external symbols' )

  return parser.parse_args()

#==============================================================================
# SCRIPT FUNCTIONALITY
#==============================================================================

def main():

    # Parse input arguments
    args = parse_input()

    # Generate list of external symbols
    create_external_list( args.err_file, args.ext_file )

#------------------------------------------------------------------------------
if __name__ == '__main__':
    # Run as main program
    main()

