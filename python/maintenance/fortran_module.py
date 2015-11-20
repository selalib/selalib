# coding: utf8
"""
Python 2 module exposing the FortranModule class, which extracts information
from a Fortran module (partially) parsed by the F2Py library.

Modules required
----------------
  * Built-in  : re
  * Library   : fortran_intrinsics
  * 3rd-party : f2py.fparser

"""
#
# Author: Yaman Güçlü, Oct 2015 - IPP Garching
#
# Last revision: 18 Nov 2015
#
from __future__ import print_function
import re
from fparser    import statements, typedecl_statements, block_statements
from fortran_intrinsics import (intrinsic_types, intrinsic_procedures,\
                                logical_operators, logical_constants, \
                                relational_operators)

__all__ = ['FortranModule','populate_exported_symbols']
__docformat__ = 'reStructuredText'

#==============================================================================
# FUNCTIONS
#==============================================================================

variable_declaration_types = \
  ( typedecl_statements.Integer,   typedecl_statements.Real,
    typedecl_statements.Complex,   typedecl_statements.Logical,
    typedecl_statements.Character, typedecl_statements.Type,
    typedecl_statements.Class,     typedecl_statements.DoublePrecision )

has_interface_types = \
( statements.ProcedureDeclaration, statements.SpecificBinding )

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
class regex_bin_oct_hex( object ):
    bin_patterns = [r"\bb'[0-1]+'"  , r'\bb"[0-1]+"'  ]
    oct_patterns = [r"\bo'[0-7]+'"  , r'\bo"[0-7]+"'  ]
    hex_patterns = [r"\bz'[a-f\d]+'", r'\bz"[a-f\d]+"']
    patterns = bin_patterns + oct_patterns + hex_patterns
    regex    = re.compile( r'({})'.format( '|'.join( patterns ) ) )

def remove_fortran_bin_oct_hex( text ):
    """ Remove Fortran literal constants in binary, octal or hexadecimal form.
    """
    return regex_bin_oct_hex.regex.sub( ' ', text )

#------------------------------------------------------------------------------
def remove_fortran_strings( text ):
    """ Remove any fortran string from the given text.
    """
    in_string = False
    delimiter = None
    new_text  = str()
    # Construct new text
    for c in text:
        # Are we in a string?
        if in_string:
            if c == delimiter:
                # Recognize that string has finished
                in_string = False
                delimiter = None
        else:
            if c in ['"',"'"]:
                # Recognize that new string has just began, and store delimiter
                in_string = True
                delimiter = c
            else:
                # Add character to new string
                new_text += c
    # Return new string
    return new_text

#------------------------------------------------------------------------------
def remove_fortran_logicals( text ):
    """
    Remove any fortran logical operator or constant from the given text.
    Also remove fortran relational operators.

    NOTE
    ----
    First, one should remove all fortran strings!

    """
    text = text.lower()
    for s in    logical_operators:  text = text.replace( s, ' ' )
    for s in    logical_constants:  text = text.replace( s, ' ' )
    for s in relational_operators:  text = text.replace( s, ' ' )
    return text

#------------------------------------------------------------------------------
def compute_external_symbols( content, fglobals=set() ):
    """
    Compute all external symbols in a 'content' list of items (usually a
    Fortran block), which represent a certain scope. This is done by searching
    all used symbols in the current scope, and substracting the globals symbols
    inherited from the parent scope. When one of the items is itself a block,
    this function is recursively called in the block, once the global variables
    of the child scope are updated with the local symbols at the current scope,
    as well as the 'associated' variables.

    Parameters
    ----------
    content : list of 'fparser' objects
      List of items (Fortran statements), usually representing a Fortran block

    fglobals : set
      Global symbols in the parent scope, which the current scope inherits

    Returns
    -------
    fexternals : set
      External symbols in the current scope (child blocks included)

    """
    # Compute locally defined symbols
    flocals = compute_locals( content )
    # Compute all used symbols
    all_used_symbols = compute_all_used_symbols( content )
    # Externals symbols at this level (not from blocks)
    fexternals = all_used_symbols - flocals - fglobals - intrinsic_procedures
    # Create new set of globals by adding the locals
    fglobals = set.union( fglobals, flocals )
    # Search inside blocks for additional external symbols
    for item in content:
        if hasattr( item, 'content' ):
            # Get associated symbols from 'associate' block
            if isinstance( item, block_statements.Associate ):
                assoc_syms = [a.associate_name for a in item.association_list]
                fglobals_block = set.union( fglobals, assoc_syms )
            # Get associated symbol (if any) from 'select type' block
            elif isinstance( item, block_statements.SelectType ):
                if item.associate_name:
                    assoc_syms = [item.associate_name]
                    fglobals_block = set.union( fglobals, assoc_syms )
            else:
                fglobals_block = fglobals
            # Recursion: find external symbols in block contents
            fexternals_block = compute_external_symbols( item.content, \
                    fglobals_block )
            # Update set of all external symbols
            fexternals.update( fexternals_block )
    # Return all external symbols
    return fexternals

#------------------------------------------------------------------------------
def compute_locals( content, return_dict=False ):
    """
    Compute local symbols from:
      1. r.h.s. of type declarations (local variables)
      2. type definitions (local types)
      3. procedure definitions (local functions and subroutines)
      4. interface definitions
      5. abstract interfaces

    Parameters
    ----------
    content : list of 'fparser' objects
      List of items (Fortran statements), usually representing a Fortran block

    return_dict : bool
      If True, a dictionary with 7 separate entries is returned

    Returns
    -------
    flocals : set | dict
      Local symbols (divided into 7 cathegories if return_dict is True)

    """
    variables   = []
    types       = []
    classes     = []
    functions   = []
    subroutines = []
    interfaces  = []
    abstract_i  = []
    # Search for local definitions and update lists
    for item in content:
        if isinstance( item, variable_declaration_types ):
            for v in item.entity_decls:
                v = v.split('(')[0].strip() # drop (:,:) from arrays
                v = v.split('=')[0].strip() # get parameter name
                variables.append( v )
        elif isinstance( item, statements.ProcedureDeclaration ):
            for v in item.proc_decls:
                v = v.split('=>')[0].strip() # drop '=> null()' from pointers
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
        elif isinstance( item, block_statements.Interface ):
            if item.isabstract:
                # Abstract interface
                abstract_i.extend( c.name for c in item.content[:-1] )
            else:
                # Non-abstract interface
                interfaces.append( item.name )

    # If required return a dictionary of sets, otherwise just one set
    if return_dict:
        return dict( variables   = set( variables   ),
                     types       = set( types       ),
                     classes     = set( classes     ),
                     functions   = set( functions   ),
                     subroutines = set( subroutines ),
                     interfaces  = set( interfaces  ),
                     abstract_i  = set( abstract_i  ), )
    else:
        return set( variables + types + classes + functions + \
                subroutines + interfaces + abstract_i )

#------------------------------------------------------------------------------
patterns = {}
patterns['name']      = r"(?<!\d\.)\b([a-zA-Z]\w*)\b"
patterns['variable']  = r"(?<![%\s])\s*" + patterns['name'] + r"\s*(?![\s\(])"
patterns['call']      = r"(?<![%\s])\s*" + patterns['name'] + r"\s*\("
patterns['extends']   = r"extends\s*\("  + patterns['name'] + r"\s*\)"
patterns['dimension'] = r"dimension\s*\("+ patterns['name'] + r"\s*\)"

re_engines = {}
for key,val in patterns.items():
    re_engines[key] = re.compile( val, re.I )

#------------------------------------------------------------------------------
def extract_expr_symbols( expr, strip=False ):
    """
    Extract all symbols that are used in a certain mathematical expression,
    distinguishing between calls and variables.

    Parameters
    ----------
    expr : str
      Mathematical expression to be parsed

    strip : bool
      If True, Fortran strings and logicals should be removed before proceeding

    Returns
    -------
    calls : list of str
      Names of functions and arrays (they cannot be distinguished)

    variables : list of str
      Names of variables (also arrays, if used as a whole)

    """
    from fparser.splitline import string_replace_map
    from fparser.utils     import specs_split_comma

    # If expression is not in lower case, convert it
    if not expr.islower():
        expr = expr.lower()

    # If required, remove Fortran strings, logicals and relationals
    if strip:
        expr = remove_fortran_bin_oct_hex( expr )
        expr = remove_fortran_strings    ( expr )
        expr = remove_fortran_logicals   ( expr )

    # Substitute argument lists with 'F2PY_EXPR_TUPLE' strings, and store map
    string, string_map = string_replace_map( expr )

    # Extract all symbols at this level
    calls     = re_engines['call'    ].findall( string )
    variables = re_engines['variable'].findall( string )
    variables = [v for v in variables if not v.startswith('F2PY_EXPR_TUPLE')]

    # Search for symbols in each of the 'mapped' strings
    for text in string_map.values():
        # Split argument lists
        for expr in specs_split_comma( text ):
            # Ignore keywords in function calls
            if '=' in expr:
                string, string_map = string_replace_map( expr )
                expr = string_map( string.rpartition('=')[2].lstrip() )
            # Recursion: extract symbols from new expression and update lists
            if expr:
                new_calls, new_variables = extract_expr_symbols( expr )
                calls    .extend( new_calls     )
                variables.extend( new_variables )

    # Return symbol lists
    return calls, variables

#------------------------------------------------------------------------------
def compute_all_used_symbols( content ):
    """
    Compute all used symbols in a certain scope level, without searching any
    child scopes.

    Parameters
    ----------
    content : list of 'fparser' objects
      List of items (Fortran statements), usually representing a Fortran block

    Returns
    -------
    used_symbols : set
      All symbols used in the current scope level

    """
    types       = []
    subroutines = []
    variables   = []
    calls       = []
    interfaces  = []

    for item in content:
        # Type declaration statements
        if isinstance( item, variable_declaration_types ):
            # Array dimensions on l.h.s.
            for s in item.attrspec:
                # TODO: what about "dimension(size(x))"?
                variables.extend( re_engines['dimension'].findall( s ) )
            # Array dimensions on r.h.s.
            for v in item.entity_decls:
                # TODO: cleanup this code
                v = remove_fortran_bin_oct_hex( v )
                v = remove_fortran_strings    ( v )  # remove strings
                v = remove_fortran_logicals   ( v )  # remove logicals
                v = v.replace('(',',')
                v = v.replace(')',',')
                syms = re_engines['variable'].findall( v )
                syms = (s for s in syms[1:] if s not in intrinsic_procedures)
                variables.extend( syms )
            # kind parameter in numerical type declarations
            if isinstance( item, (typedecl_statements.Integer,
                                  typedecl_statements.Real,
                                  typedecl_statements.Complex) ):
                kind_param = str( item.get_kind() )
                if not kind_param.isdigit():
                    variables.append( kind_param )
            # len parameter in character type declarations
            elif isinstance( item, typedecl_statements.Character ):
                len_param = str( item.get_length() )
                if (not len_param.isdigit()) and (len_param not in ['*',':']):
                    variables.append( len_param )
            # Type name in extended type declarations
            elif isinstance( item, typedecl_statements.Type ):
                types.append( item.name )
            # Type name in polymorphic extended type declarations
            elif isinstance( item, typedecl_statements.Class ):
                type_name = item.get_kind() # NOTE: this makes no sense, but...
                if type_name != '*':
                    types.append( type_name )
        # Procedure declaration statements or bindings
        elif isinstance( item, has_interface_types ):
            if item.iname:
                interfaces.append( item.iname )
        # ALLOCATE statement
        elif isinstance( item, statements.Allocate ):
            # Type name and type parameters
            type_params_rhs = []
            if item.type_spec:
                type_name, type_params_dict = item.parse_type_spec()
                if type_name not in intrinsic_types:
                    types.append( type_name )
                for key,val in type_params_dict.items():
                    type_params_rhs.append( val )
            # Optional arguments
            alloc_opt_rhs = [a.partition('=')[2] for a in item.alloc_opt_list]
            # Create unique string together with allocation list
            s = ','.join( type_params_rhs + item.allocation_list + alloc_opt_rhs )
            # Find all symbols in string
            s = remove_fortran_bin_oct_hex( s )
            s = remove_fortran_strings    ( s )
            s = remove_fortran_logicals   ( s )
            variables.extend( re_engines['variable'].findall( s ) )
            calls    .extend( re_engines['call'    ].findall( s ) )
        # TYPE GUARD statement
        elif isinstance( item, statements.TypeGuard ):
            if item.type_spec:
                types.append( item.type_spec )
        # Subroutine calls (both caller and arguments)
        elif isinstance( item, statements.Call ):
            # caller
            sub_name = item.designator
            if '%' in sub_name:
                new_calls, new_vars = extract_expr_symbols( sub_name, strip=True )
                variables.extend( new_vars  )
                calls    .extend( new_calls )
            else:
                subroutines.append( sub_name )
            # arguments
            for s in item.items:
                if not is_fortran_string( s ):
                    # Remove binary/octal/hexadecimal literal constants
                    s = remove_fortran_bin_oct_hex( s )
                    # Remove strings
                    s = remove_fortran_strings( s )
                    # Remove logicals
                    s = remove_fortran_logicals( s )
                    # Neglect argument keywords!
                    if '=' in s:
                        s = s.partition('=')[2].strip()
                    # Search for symbols
                    variables.extend( re_engines['variable'].findall( s ) )
                    calls    .extend( re_engines['call'    ].findall( s ) )
        # Assignments (both sides)
        elif isinstance( item, (statements.Assignment, 
                                statements.PointerAssignment) ):
            # l.h.s.
            variables.append( re_engines['name'].findall( item.variable )[0] )
            # r.h.s.
            new_calls, new_vars = extract_expr_symbols( item.expr, strip=True )
            variables.extend( new_vars  )
            calls    .extend( new_calls )
        # Type blocks
        elif isinstance( item, block_statements.Type ):
            if len( item.specs ) > 0:
                for s in item.specs:
                    if not is_fortran_string( s ):
                        types.extend( re_engines['extends'].findall( s ) )
        # Extract symbols from all sorts of expressions
        else:
            # IF-THEN-ELSE block
            if isinstance( item, (block_statements.If, \
                    block_statements.IfThen, statements.ElseIf) ):
                text = item.expr
            # DO loop
            elif isinstance( item, block_statements.Do ):
                text = item.item.apply_map( item.loopcontrol )
                if text.startswith( 'while' ): # DO WHILE
                    text = text[6:-1]
            # PRINT or WRITE statement
            elif isinstance( item, (statements.Print, statements.Write) ):
                text = ','.join( item.items )
            # FORALL statement
            elif isinstance( item, statements.Forall ):
                text = ','.join( item.specs[0] )
            # FORALL block
            elif isinstance( item, block_statements.Forall ):
                text = item.item.apply_map( item.specs )
            # WHERE statement
            elif isinstance( item, statements.Where ):
                text = item.expr
            # WHERE block
            elif isinstance( item, block_statements.Where ):
                text = item.item.apply_map( item.expr )
            # SELECT TYPE block
            elif isinstance( item, block_statements.SelectType ):
                text = item.selector
            # ASSOCIATE block
            elif isinstance( item, block_statements.Associate ):
                text = ','.join( a.selector for a in item.association_list )
            # Nothing of the above
            else:
                text = ''

            # Extract variables and calls from expression
            if text:
                new_calls, new_vars = extract_expr_symbols( text, strip=True )
                variables.extend( new_vars  )
                calls    .extend( new_calls )

            # Double-check that we do not have F2Py strings
            if ('f2py' in text) or ('F2PY' in text):
                print('ERROR: ', text )
                raise SystemExit()

    return set( types + subroutines + variables + calls + interfaces )

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
                self._used_modules[mod_name]['items'].append( s )
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
# Functions that operate on module objects
#==============================================================================

def populate_exported_symbols( *modules ):
    for m in modules:
        for name,data in m.used_modules:
            mod = data['object']
            mod.add_exported_symbols( data['items'] )

