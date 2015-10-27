# coding: utf8
"""
Module exposing the FortranModule class, which extracts information from a
Fortran module (partially) parsed by the F2Py library.

Modules required
----------------
  * Built-in  : re
  * Library   : fortran_intrinsics
  * 3rd-party : f2py.fparser

"""
#
# Author: Yaman Güçlü, Oct 2015 - IPP Garching
#
# Last revision: 26 Oct 2015
#
from __future__ import print_function
from fparser    import statements, typedecl_statements, block_statements
from fortran_intrinsics import intrinsic_procedures

__all__ = ['FortranModule','populate_exported_symbols']
__docformat__ = 'reStructuredText'

#==============================================================================
# FUNCTIONS
#==============================================================================

variable_declaration_types = \
  ( typedecl_statements.Integer,   typedecl_statements.Real,
    typedecl_statements.Complex,   typedecl_statements.Logical,
    typedecl_statements.Character, typedecl_statements.Type )

namespace_types = \
  ( block_statements.Type,
    block_statements.Function,
    block_statements.Subroutine )

#------------------------------------------------------------------------------
def is_fortran_string( text ):
    text  = text.strip()
    first = text[ 0]
    last  = text[-1]
    if (first == last) and (first in ["'",'"']):
        return True
    elif "//" in text:
        return True
    else:
        return False

#------------------------------------------------------------------------------
def get_external_symbols( content, fglobals=set() ):
    """ TODO: this should include
           1) all used variables
           2) all used types
           3) all used subroutines (NOT intrinsic)
           4) all calls (used functions and arrays, NOT intrinsic)
    """
    # Compute locally defined symbols
    flocals = compute_locals( content )
    # Compute all used symbols
    all_used_symbols = compute_all_used_symbols( content )
    # Externals symbols at this level (not from blocks)
    fexternals = all_used_symbols - flocals - fglobals - intrinsic_procedures
    # Recursion: search for external symbols
    for item in content:
        if hasattr( item, 'content' ):
            fexternals.update( get_external_symbols( item.content, \
                    fglobals = set.union( fglobals, flocals ) ) )
    # Return all external symbols
    return fexternals

#------------------------------------------------------------------------------
def compute_locals( content, return_dict=False ):
    """ Compute local symbols from:
          1. r.h.s. of type declarations (local variables)
          2. type definitions (local types)
          3. procedure definitions (local functions and subroutines)
          4. interface definitions (TODO)
    """
    variables   = []
    types       = []
    classes     = []
    functions   = []
    subroutines = []
    interfaces  = []
    # Search for local definitions and update lists
    for item in content:
        if isinstance( item, variable_declaration_types ):
            for v in item.entity_decls:
                v = v.split('(')[0].strip() # drop (:,:) from arrays
                variables.append( v )
        elif isinstance( item, block_statements.Type ):
            if 'abstract' in item.item.get_line(): # TODO: use regex
                classes.append( item.name )
            else:
                types.append( item.name )
        elif isinstance( item, block_statements.Function ):
            functions.append( item.name )
        elif isinstance( item, block_statements.Subroutine ):
            subroutines.append( item.name )
        else:
            # TODO: get interfaces
            pass
    # If required return a dictionary of sets, otherwise just one set
    if return_dict:
        return dict( variables   = set( variables   ),
                     types       = set( types       ),
                     classes     = set( classes     ),
                     functions   = set( functions   ),
                     subroutines = set( subroutines ),
                     interfaces  = set( interfaces  ), )
    else:
        return set( variables + types + classes + functions + \
                subroutines + interfaces )

#------------------------------------------------------------------------------
def compute_all_used_symbols( content ):
    """ Compute all used symbols from:
          1. l.h.s. of type declarations (used types)
          2. subroutine calls (used subroutines)
          3. r.h.s. of assignments (used functions and variables)
          4. l.h.s. of assignments (used variables)
          5. type blocks (entends...)
    """
    import re
    pattern_name     = r"\b([a-zA-Z]\w*)\b"
    pattern_variable = r"(?<![%\s])\s*" + pattern_name + r"\s*(?![\s\(])"
    pattern_call     = r"(?<![%\s])\s*" + pattern_name + r"\s*\("
    pattern_extends  = r"extends\s*\("  + pattern_name + r"\s*\)"

    types       = []
    subroutines = []
    variables   = []
    calls       = []

    for item in content:
        # l.h.s. of type declarations
        if isinstance( item, typedecl_statements.Type ):
            types.append( item.name )
        elif isinstance( item, typedecl_statements.Class ):
            types.append( item.get_kind() ) # NOTE: this makes no sense, but...
        # Subroutine calls (both caller and arguments)
        elif isinstance( item, statements.Call ):
            # caller
            sub_name = item.designator
            if '%' in sub_name:
                variables.append( sub_name.split('%')[0].strip() )
            else:
                subroutines.append( sub_name )
            # arguments
            for s in item.items:
                if not is_fortran_string( s ):
                    variables.extend( re.findall( pattern_variable, s ) )
                    calls    .extend( re.findall( pattern_call    , s ) )
        # Assignments (both sides)
        elif isinstance( item, (statements.Assignment, 
                                statements.PointerAssignment) ):
            # l.h.s.
            variables.append( re.findall( pattern_name, item.variable )[0] )
            # r.h.s.
            if not is_fortran_string( item.expr ):
                variables.extend( re.findall( pattern_variable, item.expr ) )
                calls    .extend( re.findall( pattern_call    , item.expr ) )
        # Type blocks
        elif isinstance( item, block_statements.Type ):
            if len( item.specs ) > 0:
                for s in item.specs:
                    if not is_fortran_string( s ):
                        types.extend( re.findall( pattern_extends, s, re.I ) )

    return set( types + subroutines + variables + calls )

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
        self._flocals = compute_locals( module.content, return_dict=True )

        # Extract external symbols (not declared anywhere in module)
        self._imported_symbols = get_external_symbols( module.content )

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
        if   symbol in self.variables  :  return True
        elif symbol in self.types      :  return True
        elif symbol in self.classes    :  return True
        elif symbol in self.functions  :  return True
        elif symbol in self.subroutines:  return True
        elif symbol in self.interfaces :  return True
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
        map( printer, self.variables )
        print()
        print( "TYPES" )
        print( "-----" )
        map( printer, self.types )
        print()
        print( "ABSTRACT TYPES" )
        print( "--------------" )
        map( printer, self.classes )
        print()
        print( "SUBROUTINES" )
        print( "-----------" )
        map( printer, self.subroutines )
        print()
        print( "FUNCTIONS" )
        print( "---------" )
        map( printer, self.functions )
        print()
        print( "INTERFACES" )
        print( "----------" )
        map( printer, self.interfaces )
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

