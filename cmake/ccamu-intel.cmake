SET(CMAKE_CXX_COMPILER     "/softs/intel/13.0/composer_xe_2013.0.079/bin/intel64/icpc" CACHE FILEPATH "")
SET(CMAKE_C_COMPILER       "/softs/intel/13.0/composer_xe_2013.0.079/bin/intel64/icc" CACHE FILEPATH "")
SET(CMAKE_Fortran_COMPILER "/softs/intel/13.0/composer_xe_2013.0.079/bin/intel64/ifort" CACHE FILEPATH "")

SET(HDF5_ENABLED ON CACHE BOOL "" FORCE)
SET(HDF5_PARALLEL_ENABLED OFF CACHE BOOL "" FORCE)
SET(HDF5_C_LIBRARY "/softs/hdf5-serial/intel/1.8.8/lib/libhdf5.a" CACHE FILEPATH "") 
SET(HDF5_FORTRAN_LIBRARY "/softs/hdf5-serial/intel/1.8.8/lib/libhdf5_fortran.a" CACHE FILEPATH "")
SET(HDF5_INCLUDE_DIRS "/softs/hdf5-serial/intel/1.8.8/include" CACHE PATH "")
SET(ZLIB_LIBRARIES "/home/cpasseron/zlib/lib/libz.a" CACHE FILEPATH "")

SET(FFTW_ENABLED ON CACHE BOOL "" FORCE)
SET(FFTW_INCLUDE_DIRS "/softs/fftw/intel/3.3/include" CACHE PATH "")
SET(FFTW_LIBRARY "/softs/fftw/intel/3.3/lib/libfftw3.a" CACHE FILEPATH "")
SET(FFTW_MPI_INCLUDE_DIR "/softs/fftw/intel/3.3/lib"  CACHE PATH "")
SET(FFTW_THREADS_LIBRARY "/softs/fftw/intel/3.3/lib/libfftw3_omp.a" CACHE FILEPATH "")
