SET(CMAKE_Fortran_COMPILER   "/opt/local/bin/gfortran-mp-4.7" CACHE FILEPATH " " FORCE )
SET(FFTW_LIBRARY             "/opt/local/lib/libfftw3.a" CACHE FILEPATH " " FORCE )
SET(CMAKE_BUILD_TYPE         Release)
SET(CMAKE_CXX_COMPILER       "/opt/local/bin/g++-mp-4.7" CACHE FILEPATH " " FORCE)
SET(CMAKE_C_COMPILER         "/opt/local/bin/gcc-mp-4.7" CACHE FILEPATH " " FORCE)
SET(FFTW_INCLUDE_DIRS        "/opt/local/include" CACHE PATH " "FORCE)
SET(FFTW_LIBRARY             "/opt/local/lib/libfftw3.a" CACHE FILEPATH " " FORCE)
SET(HDF5_PARALLEL_ENABLED    ON)
SET(HDF5_C_LIBRARY           "/opt/local/lib/libhdf5.a" CACHE FILEPATH " " FORCE)
SET(HDF5_FORTRAN_LIBRARY     "/opt/local/lib/libhdf5_fortran.a" CACHE FILEPATH " " FORCE)
SET(HDF5_INCLUDE_DIRS        "/opt/local/include" CACHE PATH " "FORCE)
SET(MPIEXEC                  "/opt/local/bin/openmpirun" CACHE FILEPATH " " FORCE)
SET(MPI_C_COMPILER           "/opt/local/bin/openmpicc" CACHE FILEPATH " " FORCE)
SET(MPI_C_INCLUDE_PATH       "/opt/local/include" CACHE PATH " "FORCE)
SET(MPI_CXX_COMPILER         "/opt/local/bin/openmpicxx" CACHE FILEPATH " " FORCE)
SET(MPI_CXX_INCLUDE_PATH     "/opt/local/include" CACHE PATH " "FORCE)
SET(MPI_Fortran_COMPILER     "/opt/local/bin/openmpif90" CACHE FILEPATH " " FORCE)
SET(MPI_Fortran_INCLUDE_PATH "/opt/local/lib" CACHE PATH " " FORCE)
SET(MPI_C_LIBRARIES          "/opt/local/lib/libmpi.dylib" CACHE FILEPATH " " FORCE)
SET(MPI_CXX_LIBRARIES        "/opt/local/lib/libmpi_cxx.dylib" CACHE FILEPATH " " FORCE)
SET(MPI_Fortran_LIBRARIES    "/opt/local/lib/libmpi_f90.dylib" CACHE FILEPATH " " FORCE)
