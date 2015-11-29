# coding: utf8
"""
Python 2 module exposing the FortranProgram class, which extracts information
from a Fortran program (partially) parsed by the F2Py library.

Modules required
----------------
  * built-in  : __future__
  * library   : content_extractor, fortran_units
  * 3rd-party : f2py.fparser

"""
#
# Author: Yaman Güçlü, Nov 2015 - IPP Garching
#
# Last revision: 24 Nov 2015
#
from __future__        import print_function
from content_extractor import compute_local_symbols, compute_external_symbols
from fortran_units     import FortranUnit
from fparser           import statements, block_statements

__all__ = ['FortranProgram']
__docformat__ = 'reStructuredText'

#==============================================================================
# CLASS: Fortran program
#==============================================================================

class FortranProgram( FortranUnit ):

    def __init__( self, filepath, program ):

        assert( isinstance( filepath, str ) )
        assert( isinstance(  program, block_statements.Program ) )
        self.filepath = filepath
        self._name    = program.name

        # Extract external symbols (not declared anywhere in program)
        self._imported_symbols = compute_external_symbols( program.content )

        # Used modules
        self._used_modules = {}
        for item in program.content:
            if type( item ) == statements.Use:
                self._used_modules[item.name] = {'isonly':item.isonly, \
                        'items':item.items, 'object':None }

    #--------------------------------------------------------------------------
    # Implementation of abstract interface
    #--------------------------------------------------------------------------
    @property
    def used_modules( self ): return self._used_modules

    @property
    def imported_symbols( self ): return self._imported_symbols

    #--------------------------------------------------------------------------
    def defines_symbol( self, symbol ):
        return False

    #--------------------------------------------------------------------------
    def add_exported_symbols( self, *symbols ):
        if symbols:
            print( "ERROR processing file '%s':" % self.filepath )
            for s in symbols:
                print( "  adding exported symbol '%s' to program" % s )
            raise SystemExit()

    #--------------------------------------------------------------------------
    def generate_interface_section( self ):
        """ Generate the interface section of a module file.
        """
        lines = []
        # Use only section
        for name,data in sorted( self._used_modules.items() ):
            tab = '! ' if name.startswith('F77_') else ''
            lines.append( tab + "use %s, only: &" % name )
            for m in sorted( data['items'] ):
                lines.append( tab + "  %s, &" % m )
            lines[-1] = lines[-1].rstrip( ', &' )
            lines.append( "" )
        # Implicit none statement
        lines.append( "implicit none" )
        lines.append( "" )
        # Concatenate line strings using newline characters
        return "\n".join( lines )

    #--------------------------------------------------------------------------
    # Additional properties/methods
    #--------------------------------------------------------------------------
    @property
    def name( self ): return self._name

    #--------------------------------------------------------------------------
    def show( self ):
        printer = lambda n: print( '  %s' % n )
        print( "================" )
        print( "IMPORTED SYMBOLS" )
        print( "================" )
        map( printer, sorted( self._imported_symbols ) )
        print()
        print( "============" )
        print( "USED MODULES" )
        print( "============" )
        for name,info in self._used_modules.items():
            if info['isonly']:
                print( '\n%s:' % name )
                map( printer, info['items'] )
            else:
                print( '\n%s: (all)' % name )
        print()
