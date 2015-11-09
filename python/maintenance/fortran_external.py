import re
import fortran_iso_c_binding
import fortran_iso_fortran_env
import fortran_fftpack

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

fftpack = {}
fftpack['lib_name'] = 'fftpack'
fftpack['mod_name'] = ''
fftpack['all_syms'] = fortran_fftpack.procedures
fftpack['match'   ] = lambda s : s.lower() in fftpack['all_syms']

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

library_collection = [fftpack, iso_c_binding, iso_fortran_env, fftw, mpi, hdf5]

def find_external_library( symbol ):
    for lib in library_collection:
        if lib['match']( symbol ):
            return lib['lib_name'], lib['mod_name']
    return ()
