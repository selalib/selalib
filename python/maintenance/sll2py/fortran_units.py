# coding: utf8
"""
Python 2 module exposing the FortranUnitBase abstract class, which should be
subclassed by FortranModule, FortranInterfaceModule, FortranProgram

Modules required
----------------
  * Built-in  : abc

"""
#
# Author: Yaman Güçlü, Nov 2015 - IPP Garching
#
# Last revision: 24 Nov 2015
#

from abc import ABCMeta, abstractmethod, abstractproperty

__all__ = ['FortranUnitBase','FortranUnit']

#==============================================================================
# ABSTRACT BASE CLASS
#==============================================================================

class FortranUnitBase( object ):

    __metaclass__ = ABCMeta 

    @abstractmethod
    def link_used_modules( self, modules, externals={} ):
        """ Given a list of library modules, link against the used ones. """

    @abstractmethod
    def update_use_statements( self, 
            find_external_library = None,
            ignored_symbols       = [] ):
        """ Update all use statements. """

    @abstractmethod
    def cleanup_use_statements( self ):
        """ Remove useless modules, remove duplicate items, add only. """

    @abstractmethod
    def scatter_imported_symbols( self ):
        """ Notify used modules that symbols here imported should be public."""

    @abstractmethod
    def generate_interface_section( self ):
        """ Generate the interface section for the file. """

#==============================================================================
# ABSTRACT SUBCLASS: Fortran unit
#==============================================================================

class FortranUnit( FortranUnitBase ):
    # NOTE: this is an abstract class that implements 4 of the 5 abstract
    # methods in "FortranUnitBase"

    #--------------------------------------------------------------------------
    # Interface for derived classes
    #--------------------------------------------------------------------------
    @abstractproperty
    def used_modules( self ):
        """ Dictionary with used modules. """

    @abstractproperty
    def imported_symbols( self ):
        """ List of all imported symbols (strings). """

    @abstractmethod
    def defines_symbol( self, symbol ):
        """ [Modules only] Return True if the symbol is defined locally. """

    @abstractmethod
    def add_exported_symbols( self, *symbols ):
        """ [Modules only]
            Notify module that symbols are used outside and should be public.
        """

    @abstractmethod
    def generate_interface_section( self ):
        """ Generate the interface section for the file. """

    #--------------------------------------------------------------------------
    # Implemented methods
    #--------------------------------------------------------------------------
    def find_symbol_def( self, symbol ):
        """ Find name of module where symbol is defined (recursive search).
        """
        # Search in current unit
        if self.defines_symbol( symbol ):
            # Symbol is defined locally
            return self.name

        # Search in used modules
        for name,data in self.used_modules.items():
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
    def link_used_modules( self, modules, externals={} ):
        """ Given a list of library modules, link against the used ones.
        """
        for name,data in self.used_modules.items():

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
    def update_use_statements( self,
            find_external_library = None,
            ignored_symbols       = [] ):
        """
        Update all use statements.

        """
        for s in self.imported_symbols:

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
            if mod_name in self.used_modules.keys():
                item_list = self.used_modules[mod_name]['items']
                if s not in item_list:
                    item_list.append( s )
            else:
                new_mod = { 'isonly': True, 'items': [s], 'object': None }
                self.used_modules[mod_name] = new_mod

    #--------------------------------------------------------------------------
    def cleanup_use_statements( self ):
        """ Remove useless modules, remove duplicate items, add only
        """
        # Set 'isonly=True' for all used modules, remove duplicate symbols,
        # and identify useless modules
        useless = []
        for name,data in self.used_modules.items():
            data['isonly'] = True
            data['items' ] = tuple( set( data['items'] ) )
            if data['items'] == ():
                useless.append( name )
        # Remove useless modules
        for m in useless:
            self.used_modules.pop( m )

    #--------------------------------------------------------------------------
    def scatter_imported_symbols( self ):
        """ Notify used modules that symbols here imported should be public.
        """
        for name,data in self.used_modules.items():
            m = data['object']
            if m is not None:
                m.add_exported_symbols( *data['items'] )
