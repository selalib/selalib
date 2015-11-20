# coding: utf8
"""
Python 2 module exposing the FortranModule class, which extracts information
from a Fortran module (partially) parsed by the F2Py library.

Modules required
----------------
  * 3rd-party : f2py.fparser

"""
#
# Author: Yaman Güçlü, Oct 2015 - IPP Garching
#
# Last revision: 20 Nov 2015
#
from __future__        import print_function
from fparser           import statements, block_statements
from content_extractor import compute_local_symbols, compute_external_symbols

from fortran_units import FortranUnitBase

__all__ = ['FortranModule','LibraryInterfaceModule','NewFortranModule']
__docformat__ = 'reStructuredText'

#==============================================================================
# CLASS: Fortran module
#==============================================================================

class FortranModule( FortranUnitBase ):

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
    def link_used_modules( self, modules, externals={} ):
        """ Given a list of library modules, link against the used ones.
        """
        for name,data in self._used_modules.items():

            if name in externals:
                print( "WARNING: skipping external module %s" % name )
                data['object']   = None
                data['external'] = True
                continue

            objects = [m for m in modules if m.name == name]

            if len( objects ) == 0:
                print( "ERROR: missing link for module %s" % name )
                raise SystemExit()
            elif len( objects ) == 1:
                data['object']   = objects[0]
                data['external'] = False
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
        elif symbol in self.abstract_i :  return True
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

        # Search in used modules
        for name,data in self._used_modules.items():
            # Search in "use only" symbol list
            if data['isonly']:
                if symbol in data['items']:
                    # Symbol is found in "use [], only:" list
                    return name
                continue
            # Recursively call same function in "use all" modules (if linked!)
            if data['object']:
                mod_name = data['object'].find_symbol_def( symbol )
                if mod_name:
                    # Symbol is found in some other module in the use hierarchy
                    return mod_name

        # Nothing was found
        return None

    #--------------------------------------------------------------------------
    def update_use_statements( self,
            find_external_library = None,
            ignored_symbols       = [] ):
        """
        Update all use statements.

        """
        for s in self._imported_symbols:

            # Recursively search for symbol in used modules
            mod_name = self.find_symbol_def( s )

            # Symbol was NOT FOUND. If available, search in external libraries
            if (mod_name is None) and (find_external_library is not None):
                external_match = find_external_library( s )
                if external_match:
                    # Get library and module name
                    lib_name, mod_name = external_match
                    if mod_name == '':
                        # If module name is missing, F77 interface is used
                        mod_name = 'F77_' + lib_name

            # Symbol was NOT FOUND. If it should be ignored, print warning and
            # move to next symbol. Otherwise, raise error
            if mod_name is None:
                if s in ignored_symbols:
                    print( "WARNING: ignoring symbol '%s'" % s )
                    continue
                else:
                    print( "ERROR processing file '%s':" % self.filepath )
                    print( "  cannot locate symbol '%s'" % s )
                    raise SystemExit()

            # Symbol was FOUND in module 'mod_name'. If 'mod_name' is already
            # in the list of used modules, add symbol to the module item list.
            # Otherwise, add 'mod_name' to the list of used modules
            if mod_name in self._used_modules.keys():
                item_list = self._used_modules[mod_name]['items']
                if s not in item_list:
                    item_list.append( s )
            else:
                new_mod = { 'isonly': True, 'items': [s], 'object': None }
                self._used_modules[mod_name] = new_mod

    #--------------------------------------------------------------------------
    def cleanup_use_statements( self ):
        """ Remove useless modules, remove duplicate items, add only
        """
        # Set 'isonly=True' for all used modules, remove duplicate symbols,
        # and identify useless modules
        useless = []
        for name,data in self._used_modules.items():
            data['isonly'] = True
            data['items' ] = tuple( set( data['items'] ) )
            if data['items'] == ():
                useless.append( name )
        # Remove useless modules
        for m in useless:
            self._used_modules.pop( m )

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
    def scatter_imported_symbols( self ):
        """ Notify used modules that symbols here imported should be public.
        """
        for name,data in self._used_modules.items():
            m = data['object']
            if m is not None:
                m.add_exported_symbols( *data['items'] )

    #--------------------------------------------------------------------------
    def generate_interface_section( self ):
        """ Generate the interface section of a module file.
        """
        lines = []
        # Use only section
        for name,data in self._used_modules.items():
            tab = '! ' if name.startswith('F77_') else ''
            lines.append( tab + "use %s, only: &" % name )
            for m in data['items']:
                lines.append( tab + "  %s, &" % m )
            lines[-1] = lines[-1].rstrip( ', &' )
            lines.append( "" )
        # Implicit none statement
        lines.append( "implicit none" )
        lines.append( "" )
        # Public section
        if self._exported_symbols:
            lines.append( "public :: &" )
            for s in self._exported_symbols:
                lines.append( "  %s, &" % s )
            lines[-1] = lines[-1].rstrip( ', &' )
            lines.append( "" )
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
    def  abstract_i( self ): return self._flocals['abstract_i']
    
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

class LibraryInterfaceModule( FortranModule ):

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
        return LibraryInterfaceModule( filepath, module )
    else:
        return mmod
