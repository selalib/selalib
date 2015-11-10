import re
import fortran_iso_c_binding
import fortran_iso_fortran_env
import fortran_deboor
import fortran_fftpack
import fortran_fishpack
import fortran_mudpack

#==============================================================================
# Fortran extension modules
#==============================================================================

iso_c_binding = {}
iso_c_binding['lib_name'] = 'iso_c_binding'
iso_c_binding['mod_name'] = 'iso_c_binding' 
iso_c_binding['all_syms'] = set.union( \
        fortran_iso_c_binding.intrinsic_procedures,
        fortran_iso_c_binding.named_constants,
        fortran_iso_c_binding.parameters,
        )
iso_c_binding['all_syms'] = set( s.lower() for s in iso_c_binding['all_syms'] )
iso_c_binding['match']    = lambda s : s.lower() in iso_c_binding['all_syms']

iso_fortran_env = {}
iso_fortran_env['lib_name'] = 'iso_fortran_env'
iso_fortran_env['mod_name'] = 'iso_fortran_env'
iso_fortran_env['all_syms'] = set.union( \
        fortran_iso_fortran_env.named_constants,
        fortran_iso_fortran_env.derived_types,
        fortran_iso_fortran_env.intrinsic_procedures,
        )
iso_fortran_env['all_syms'] = set( s.lower() for s in iso_fortran_env['all_syms'] )
iso_fortran_env['match']    = lambda s : s.lower() in iso_fortran_env['all_syms']

#==============================================================================
# Libraries shipped with SeLaLib (inside 'external' directory)
#==============================================================================

deboor = {}
deboor['lib_name'] = 'deboor'
deboor['mod_name'] = ''
deboor['all_syms'] = fortran_deboor.procedures
deboor['match'   ] = lambda s : s.lower() in deboor['all_syms']

fftpack = {}
fftpack['lib_name'] = 'fftpack'
fftpack['mod_name'] = ''
fftpack['all_syms'] = fortran_fftpack.procedures
fftpack['match'   ] = lambda s : s.lower() in fftpack['all_syms']

fishpack = {}
fishpack['lib_name'] = 'fishpack'
fishpack['mod_name'] = ''
fishpack['all_syms'] = fortran_fishpack.procedures
fishpack['match'   ] = lambda s : s.lower() in fishpack['all_syms']

mudpack = {}
mudpack['lib_name'] = 'mudpack'
mudpack['mod_name'] = ''
mudpack['all_syms'] = fortran_mudpack.procedures
mudpack['match'   ] = lambda s : s.lower() in mudpack['all_syms']

#==============================================================================
# Libraries installed on the system
#==============================================================================

fftw = {}
fftw['lib_name'] = 'fftw'
fftw['mod_name'] = 'sll_m_fftw3'
fftw['match'   ] = re.compile( '^d?fftw_\w+\Z', re.I ).match

mpi = {}
mpi['lib_name'] = 'mpi'
mpi['mod_name'] = 'mpi'
mpi['match']    = re.compile( '^mpi_\w+\Z', re.I ).match

hdf5 = {}
hdf5['lib_name'] = 'hdf5'
hdf5['mod_name'] = 'hdf5'
hdf5['match']    = re.compile( '^h5\w+_f\Z', re.I ).match

#==============================================================================
# Convenience objects
#==============================================================================

library_collection = [iso_c_binding, iso_fortran_env, deboor, fftpack,
        fishpack, mudpack, fftw, mpi, hdf5]

external_modules = \
        {lib['mod_name'] for lib in library_collection if lib['mod_name']}

external_f77 = \
        {lib['lib_name'] for lib in library_collection if not lib['mod_name']}

# Add 'match' function to dictionaries that do not have it
for lib in library_collection:
    if lib.has_key('all_syms') and not lib.has_key('match'):
        lib['match'] = lambda s : s.lower() in lib['all_syms']

#==============================================================================
# Convenience functions
#==============================================================================

def find_external_library( symbol ):
    """ Find name of library (and module) containing the given symbol
    """
    for lib in library_collection:
        if lib['match']( symbol ):
            return lib['lib_name'], lib['mod_name']
    return ()
