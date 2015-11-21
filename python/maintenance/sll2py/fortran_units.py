# coding: utf8
"""
Python 2 module exposing the FortranUnitBase abstract class, which should be
subclassed by FortranModule, FortranModuleInterface, FortranProgram

Modules required
----------------
  * Built-in  : abc

"""
#
# Author: Yaman Güçlü, Nov 2015 - IPP Garching
#
# Last revision: 20 Nov 2015
#

from abc import ABCMeta, abstractmethod

__all__ = ['FortranUnitBase']

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
        """ Remove useless modules, remove duplicate items, add only
        """

    @abstractmethod
    def scatter_imported_symbols( self ):
        """ Notify used modules that symbols here imported should be public.
        """

    @abstractmethod
    def generate_interface_section( self ):
        """ Generate the interface section for the file.
        """

#==============================================================================

