# coding: utf8
"""
Python 2 module exposing the FortranModule class, which extracts information
from a Fortran module (partially) parsed by the F2Py library.

Modules required
----------------
  * built-in  : __future__
  * library   : content_extractor, fortran_units
  * 3rd-party : f2py.fparser

"""
#
# Author: Yaman Güçlü, Oct 2015 - IPP Garching
#
# Last revision: 24 Nov 2015
#
from __future__        import print_function
from content_extractor import compute_local_symbols, compute_external_symbols
from fortran_units     import FortranUnit
from fparser           import statements, block_statements

__all__ = ['FortranModule','FortranInterfaceModule','NewFortranModule']
__docformat__ = 'reStructuredText'

#==============================================================================
# CLASS: Fortran module
#==============================================================================

class FortranModule( FortranUnit ):

    def __init__( self, filepath, module ):

        assert( isinstance( filepath, str ) )
        assert( isinstance(   module, block_statements.Module ) )
        self.filepath = filepath
        self._name    = module.name

        # Extract module entities
        self._flocals = compute_local_symbols( module.content, return_dict=True )

        # Extract external symbols (not declared anywhere in module)
        self._imported_symbols = compute_external_symbols( module.content )

        # Set: exported symbols (empty for now)
        self._exported_symbols = set()

        # Used modules
        self._used_modules = {}
        for item in module.content:
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
        """ Return True if the symbol is defined in the module.
        """
        assert( isinstance( symbol, str ) )
        if   symbol in self.variables  :  return True
        elif symbol in self.types      :  return True
        elif symbol in self.classes    :  return True
        elif symbol in self.functions  :  return True
        elif symbol in self.subroutines:  return True
        elif symbol in self.interfaces :  return True
        elif symbol in self.abstract_i :  return True
        else:
            return False

    #--------------------------------------------------------------------------
    def add_exported_symbols( self, *symbols ):
        """ Notify module that symbols are used outside and should be public.
        """
        for s in symbols:
            if not self.defines_symbol( s ):
                print( "ERROR processing file '%s':" % self.filepath )
                print( "  adding symbol '%s' that is not defined here" % s )
                raise SystemExit()
            self._exported_symbols.add( s )

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
        # Public section
        if self._exported_symbols:
            lines.append( "public :: &" )
            for s in sorted( self._exported_symbols ):
                lines.append( "  %s, &" % s )
            lines[-1] = lines[-1].rstrip( ', &' )
            lines.append( "" )
        # Private blanket statement
        lines.append( "private" )
        # Concatenate line strings using newline characters
        return "\n".join( lines )

    #--------------------------------------------------------------------------
    # Additional properties
    #--------------------------------------------------------------------------
    @property
    def name( self ): return self._name

    @property
    def      locals( self ): return self._flocals
    @property
    def   variables( self ): return self._flocals['variables']
    @property
    def       types( self ): return self._flocals['types']
    @property
    def     classes( self ): return self._flocals['classes']
    @property
    def   functions( self ): return self._flocals['functions']
    @property
    def subroutines( self ): return self._flocals['subroutines']
    @property
    def  interfaces( self ): return self._flocals['interfaces']
    @property
    def  abstract_i( self ): return self._flocals['abstract_i']
    
    @property
    def exported_symbols( self ): return tuple( self._exported_symbols )

    #--------------------------------------------------------------------------
    # Additional methods
    #--------------------------------------------------------------------------
    def show( self ):
        printer = lambda n: print( '  %s' % n )
        print( "==================" )
        print( "MODULE DEFINITIONS" )
        print( "==================" )
        print()
        print( "VARIABLES" )
        print( "---------" )
        map( printer, sorted( self.variables ) )
        print()
        print( "TYPES" )
        print( "-----" )
        map( printer, sorted( self.types ) )
        print()
        print( "ABSTRACT TYPES" )
        print( "--------------" )
        map( printer, sorted( self.classes ) )
        print()
        print( "SUBROUTINES" )
        print( "-----------" )
        map( printer, sorted( self.subroutines ) )
        print()
        print( "FUNCTIONS" )
        print( "---------" )
        map( printer, sorted( self.functions ) )
        print()
        print( "INTERFACES" )
        print( "----------" )
        map( printer, sorted( self.interfaces ) )
        print()
        print( "ABSTRACT INTERFACES" )
        print( "-------------------" )
        map( printer, sorted( self.abstract_i ) )
        print()

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

#==============================================================================
# CLASS: Selalib library interface (= Fortran module with no definitions)
#==============================================================================

class FortranInterfaceModule( FortranModule ):

    def __init__( self, filepath, module ):
        FortranModule.__init__( self, filepath, module )
        self.used_symbols = []
        for name,data in self._used_modules.items():
            self.used_symbols.extend( data['items'] )

    #--------------------------------------------------------------------------
    def add_exported_symbols( self, *symbols ):
        """ Notify module that symbols are used outside and should be public.
        """
        for s in symbols:
            if s not in self.used_symbols:
                print( "ERROR processing file '%s':" % self.filepath )
                print( "  adding symbol '%s' that is not present here" % s )
                raise SystemExit()

    #--------------------------------------------------------------------------
    def generate_interface_section( self ):
        """ Generate the interface section of a module file.
        """
        print( "WARNING: no interface section for library interface module" )
        return ''

#==============================================================================
# FACTORY FUNCTION: returns FortranModule or LibraryInterfaceModule
#==============================================================================

def NewFortranModule( filepath, module ):
    name = module.name
    mmod = FortranModule( filepath, module )
    nsym = len( set.union( *mmod.locals.values() ) )
    if name.startswith('sll_') and not name.startswith('sll_m_') and nsym == 0:
        return FortranInterfaceModule( filepath, module )
    else:
        return mmod
