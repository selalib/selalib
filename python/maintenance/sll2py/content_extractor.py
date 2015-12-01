# coding: utf8
"""
This Python 2 module exposes two functions, which extract information
from a Fortran module (partially) parsed by the F2Py library.

Modules required
----------------
  * Built-in  : re
  * Library   : fortran_intrinsics
  * 3rd-party : f2py.fparser

"""
#
# Author: Yaman Güçlü, Nov 2015 - IPP Garching
#
# Last revision: 01 Dec 2015
#
from __future__ import print_function
import re
from fparser    import statements, typedecl_statements, block_statements
from fortran_intrinsics import (intrinsic_types, intrinsic_procedures,\
                                logical_operators, logical_constants, \
                                relational_operators)

__all__ = ['compute_local_symbols','compute_external_symbols']
__docformat__ = 'reStructuredText'

#==============================================================================
# LOW-LEVEL FUNCTIONS (UTILITIES)
#==============================================================================

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

#==============================================================================
# MID-LEVEL FUNCTIONS (NOT EXPOSED)
#==============================================================================

# Lists of Fparser types
variable_declaration_types = \
  ( typedecl_statements.Integer,   typedecl_statements.Real,
    typedecl_statements.Complex,   typedecl_statements.Logical,
    typedecl_statements.Character, typedecl_statements.Type,
    typedecl_statements.Class,     typedecl_statements.DoublePrecision )

has_interface_types = \
( statements.ProcedureDeclaration, statements.SpecificBinding )

#------------------------------------------------------------------------------
# Regex patterns
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

    # Empty symbol lists, to be populated
    calls     = []
    variables = []

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

    # Split argument lists, if any
    for text in specs_split_comma( string ):
        # Ignore keywords in function calls
        if '=' in text:
            text = text.rpartition('=')[2].lstrip()
        # Extract all symbols at this level
        if text:
            new_calls     = re_engines['call'    ].findall( text )
            new_variables = re_engines['variable'].findall( text )
            new_variables = [v for v in new_variables \
                    if not v.startswith('F2PY_EXPR_TUPLE')]
            calls    .extend( new_calls     )
            variables.extend( new_variables )

    # Recursion: search for symbols in the 'mapped' strings, and update lists
    for expr in string_map.values():
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
                    text = text[5:].strip()[1:-1].strip()
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
            # SELECT CASE block
            elif isinstance( item, block_statements.Select ):
                text = item.item.apply_map( item.expr )
            # CASE statement
            elif isinstance( item, statements.Case ):
                text = ','.join( s for c in item.items for s in c )
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
# HIGH-LEVEL FUNCTIONS (EXPOSED)
#==============================================================================

def compute_local_symbols( content, return_dict=False ):
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
    flocals = compute_local_symbols( content )
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
