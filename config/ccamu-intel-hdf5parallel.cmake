SET(CMAKE_BUILD_TYPE "Release" CACHE STRING " " FORCE)
SET(CMAKE_Fortran_COMPILER "/softs/intel/13.0/composer_xe_2013.0.079/bin/intel64/ifort" CACHE FILEPATH " " FORCE)

SET(HDF5_ENABLED ON CACHE BOOL "" FORCE)
SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "" FORCE)
SET(HDF5_INCLUDE_DIRS "/softs/hdf5/intel/1.8.8/include" CACHE PATH "" FORCE)
SET(HDF5_C_LIBRARY "/softs/hdf5/intel/1.8.8/lib/libhdf5.a" CACHE FILEPATH "" FORCE) 
SET(HDF5_FORTRAN_LIBRARY "/softs/hdf5/intel/1.8.8/lib/libhdf5_fortran.a" CACHE FILEPATH "" FORCE)

SET(FFTW_ENABLED ON CACHE BOOL "" FORCE)
SET(FFTW_INCLUDE_DIRS "/softs/fftw/openmpi/intel/3.3/include" CACHE PATH " " FORCE)
SET(FFTW_LIBRARY "/softs/fftw/openmpi/intel/3.3/lib/libfftw3.a" CACHE FILEPATH " " FORCE)
SET(FFTW_MPI_INCLUDE_DIR "/softs/fftw/openmpi/intel/3.3/lib"  CACHE PATH " " FORCE)
