# coding: utf8
"""
Module exposing the FortranModule class, which extracts information from a
Fortran module (partially) parsed by the F2Py library.

Modules required
----------------
  * Built-in  : re
  * 3rd-party : f2py.fparser

"""
#
# Author: Yaman Güçlü, Oct 2015 - IPP Garching
#
# Last revision: 23 Oct 2015
#
from __future__ import print_function
from   fparser  import statements, typedecl_statements, block_statements

__all__ = ['FortranModule','populate_exported_symbols']
__docformat__ = 'reStructuredText'

#==============================================================================
# FUNCTIONS
#==============================================================================

variable_declaration_types = \
  ( typedecl_statements.Integer, typedecl_statements.Real,
    typedecl_statements.Complex, typedecl_statements.Logical,
    typedecl_statements.Type )

namespace_types = \
  ( block_statements.Type,
    block_statements.Function,
    block_statements.Subroutine )

#------------------------------------------------------------------------------
def get_all_used_types( *content ):
    """ Extract ALL used types from items in content.
    """
    for item in content:
        if isinstance( item, typedecl_statements.Type ):
            yield item.name
        elif isinstance( item, typedecl_statements.Class ):
            yield item.get_kind() # NOTE: this makes no sense, but...
        elif isinstance( item, namespace_types ):
            for t in get_all_used_types( *item.content ):
                yield t
            # In Python3.3+:
            # yield from get_all_used_types( *item.content )

#------------------------------------------------------------------------------
def get_all_calls( *content ):
    """ Extract all non type-bound calls from items in content.
    """
    import re
#    pattern = r"([A-Za-z]\w*) *\("
    pattern = r"(?<![%\s])\s*(\b[A-Za-z]\w*\b) *\("  # discard type-bound call
    regex   = re.compile( pattern )
    for item in content:
        # TODO: do not parse strings!!!!!
        # If assignment, parse r.h.s.
        if isinstance( item, statements.Assignment ):
            for match in regex.findall( item.expr ):
                yield match
        # If pointer assignment, parse r.h.s.
        elif isinstance( item, statements.PointerAssignment ):
            for match in regex.findall( item.expr ):
                yield match
        # If subroutine call, parse arguments
        elif isinstance( item, statements.Call ):
            for s in item.items:
                for match in regex.findall( s ):
                    yield match
        # If block, search recursively:
        elif hasattr( item, 'content' ):
            for c in get_all_calls( *item.content ):
                yield c
            # In Python3.3+:
            # yield from get_all_calls( *item.content )

#------------------------------------------------------------------------------
def get_all_subroutine_calls( *content ):
    """ Extract all non type-bound subroutine calls from items in content.
    """
    for item in content:
        # If subroutine call, return its name
        if isinstance( item, statements.Call ):
            if '%' not in item.designator: # discard type-bound subroutine call
                yield item.designator
        # If this is a block, search recursively:
        elif hasattr( item, 'content' ):
            for c in get_all_subroutine_calls( *item.content ):
                yield c
            # In Python3.3+:
            # yield from get_all_subroutine_calls( *item.content )

#------------------------------------------------------------------------------
# Extract names of declared variables
def get_declared_variables_names( *content ):
    for item in content:
        if is_variable_typedecl( item ):
            for var in item.entity_decls:
                yield var

#==============================================================================
# CLASSES
#==============================================================================

class FortranModule( object ):

    def __init__( self, filepath, module ):

        assert( isinstance( filepath, str ) )
        assert( isinstance(   module, block_statements.Module ) )
        self.filepath = filepath
        self._name    = module.name

        # Extract module entities
        self._variables   = []
        self._types       = []
        self._classes     = []
        self._functions   = []
        self._subroutines = []
        self._interfaces  = []
        for item in module.content:
            if isinstance( item, variable_declaration_types ):
                self._variables.extend( item.entity_decls )
            elif isinstance( item, block_statements.Type ):
                if 'abstract' in item.item.get_line(): # TODO: use regex
                    self._classes.append( item.name )
                else:
                    self._types.append( item.name )
            elif isinstance( item, block_statements.Function ):
                self._functions.append( item.name )
            elif isinstance( item, block_statements.Subroutine ):
                self._subroutines.append( item.name )
            else:
                # TODO: get interfaces
                pass
        # Create tuples, sorted
        self._variables   = tuple( sorted( self._variables   ) )
        self._types       = tuple( sorted( self._types       ) )
        self._classes     = tuple( sorted( self._classes     ) )
        self._functions   = tuple( sorted( self._functions   ) )
        self._subroutines = tuple( sorted( self._subroutines ) )
        self._interfaces  = tuple( sorted( self._interfaces  ) )

        # Extract external types (also abstract)
        ext_types = set( get_all_used_types( *module.content ) ) - \
                    set( self._types ) - set( self._classes )
        # Extract external calls (functions and arrays)
        ext_calls = set( get_all_calls( *module.content ) ) - \
                    set( self._functions )
        # Extract external subroutine calls
        ext_sub_calls = set( get_all_subroutine_calls( *module.content ) ) - \
                        set( self._subroutines )
        # TODO: extract external module variables

        # Get a list of all imported symbols
        self._imported_symbols = tuple( sorted( ext_types ) + \
                sorted( ext_calls ) + sorted( ext_sub_calls ) )

        # Set: exported symbols (empty for now)
        self._exported_symbols = set()

        # Used modules
        self._used_modules = {}
        for item in module.content:
            if type( item ) == statements.Use:
                self._used_modules[item.name] = {'isonly':item.isonly, \
                                                 'items' :item.items }

    #--------------------------------------------------------------------------
    def link_used_modules( self, *modules ):
        """ Given a list of library modules, link against the used ones.
        """
        for name,data in self._used_modules.items():
            objects = []
            for m in modules:
                if m.name == name:
                    objects.append( m )
            if len( objects ) == 0:
                print( "WARNING: missing link for module %s" % name )
                data['object'] = None
            elif len( objects ) == 1:
                data['object'] = m
            else:
                print( "ERROR: multiple modules with name %s" % name )
                raise SystemExit()

    #--------------------------------------------------------------------------
    def defines_symbol( self, symbol ):
        """ Return True if the symbol is defined in the module.
        """
        assert( isinstance( symbol, str ) )
        if   symbol in self._variables  :  return True
        elif symbol in self._types      :  return True
        elif symbol in self._classes    :  return True
        elif symbol in self._functions  :  return True
        elif symbol in self._subroutines:  return True
        elif symbol in self._interfaces :  return True
        else:
            return False

    #--------------------------------------------------------------------------
    def find_symbol_def( self, symbol ):
        """ Find name of module where symbol is defined (recursive search).
        """
        # Search in current module
        if self.defines_symbol( symbol ):
            # Symbol is found in module
            return self.name

        # Search in all used modules
        for name,data in self._used_modules.items():
            if data['isonly']:
                if symbol in data['items']:
                    # Symbol is found in "use [], only:" list
                    return name
                continue
            if data['object'] is None:
                print( "WARNING: missing link for module %s" % name )
                continue
            # Recursively call same function
            mod_name = data['object'].find_symbol_def( symbol )
            if mod_name is not None:
                # Symbol is found in some other module in the use hierarchy
                return mod_name

        # Nothing was found
        return None

    #--------------------------------------------------------------------------
    def update_use_statements( self ):
        """ Update all use statements.
        """
        unlocated_symbols = list( self._imported_symbols )

        for s in self._imported_symbols:
            mod_name = self.find_symbol_def( s )
            if mod_name is None:
                print( "ERROR: cannot locate symbol %s" % s )
                raise SystemExit()

            if mod_name in self._used_modules.keys():
                self._used_modules[mod_name]['items'].append( s )
                unlocated_symbols.remove( s )
            else:
                new_mod = { 'isonly': True, 'items': [s] }
                self._used_module['mod_name'] = new_mod
                unlocated_symbols.remove( s )

        assert( unlocated_symbols == [] )

        # Remove useless modules, remove duplicate items, add only
        # NOTE: modules from external libraries will give problems
        useless = []
        for name,data in self._used_modules.items():
            data['isonly'] = True
            data['items' ] = tuple( set( data['items'] ) )
            if data['items'] == ():
                useless.append( name )
        for m in useless:
            self._used_modules.pop( m )

    #--------------------------------------------------------------------------
    def add_exported_symbols( self, *symbols ):
        """ Notify module that symbols are used outside and should be public.
        """
        for s in symbols:
            if not self.defines_symbol( symbol ):
                print( "ERROR: adding symbol that is not defined in module." )
                raise SystemExit()
            self._exported_symbols.add( symbol )

    #--------------------------------------------------------------------------
    def generate_interface_section( self ):
        """ Generate the interface section of a module file.
        """
        lines = []
        # Use only section
        for name,data in self._used_modules.items():
            lines.append( "use %s, only: &" % name )
            for m in data['items']:
                lines.append( "  %s, &" % m )
            lines[-1] = lines[-1].rstrip( ', &' )
            lines.append( "" )
        # Implicit none statement
        lines.append( "implicit none" )
        lines.append( "" )
        # Public section
        lines.append( "public :: &" )
        for s in self._exported_symbols:
            lines.append( "  %s, &" )
        lines[-1] = lines[-1].rstrip( ', &' )
        # Private blanket statement
        lines.append( "private" )
        # Concatenate line strings using newline characters
        return "\n".join( lines )

    #--------------------------------------------------------------------------
    @property
    def name( self ): return self._name

    @property
    def   variables( self ): return self._variables
    @property
    def       types( self ): return self._types
    @property
    def     classes( self ): return self._classes
    @property
    def   functions( self ): return self._functions
    @property
    def subroutines( self ): return self._subroutines
    @property
    def  interfaces( self ): return self._interfaces
    
    @property
    def imported_symbols( self ): return self._imported_symbols

    @property
    def exported_symbols( self ): return tuple( self._exported_symbols )

    @property
    def used_modules( self ): return dict( self._used_modules )

    #--------------------------------------------------------------------------
    def show( self ):
        printer = lambda n: print( '  %s' % n )
        print( "==================" )
        print( "MODULE DEFINITIONS" )
        print( "==================" )
        print()
        print( "VARIABLES" )
        print( "---------" )
        map( printer, self._variables )
        print()
        print( "TYPES" )
        print( "-----" )
        map( printer, self._types )
        print()
        print( "ABSTRACT TYPES" )
        print( "--------------" )
        map( printer, self._classes )
        print()
        print( "SUBROUTINES" )
        print( "-----------" )
        map( printer, self._subroutines )
        print()
        print( "FUNCTIONS" )
        print( "---------" )
        map( printer, self._functions )
        print()
        print( "INTERFACES" )
        print( "----------" )
        map( printer, self._interfaces )
        print()

        print( "================" )
        print( "IMPORTED SYMBOLS" )
        print( "================" )
        map( printer, self._imported_symbols )
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
# Functions that operate on module objects
#==============================================================================

def populate_exported_symbols( *modules ):
    for m in modules:
        for name,data in m.used_modules:
            mod = data['object']
            mod.add_exported_symbols( data['items'] )

