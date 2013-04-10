SET(CMAKE_BUILD_TYPE "Release" CACHE STRING " " FORCE)
SET(CMAKE_Fortran_FLAGS "-I/home/vgdgirard/local/include" CACHE STRING " " FORCE) 

SET(CMAKE_CXX_COMPILER "/softs/intel/13.0/bin/icpc" CACHE FILEPATH " " FORCE)
SET(CMAKE_C_COMPILER "/softs/intel/13.0/bin/icc" CACHE FILEPATH " " FORCE)
SET(CMAKE_Fortran_COMPILER "/softs/intel/13.0/bin/ifort" CACHE FILEPATH " " FORCE)

SET(HDF5_ENABLED ON CACHE BOOL "" FORCE)
SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "" FORCE)
SET(HDF5_INCLUDE_DIRS "/home/vgdgirard/local/include" CACHE PATH "" FORCE)
SET(HDF5_C_LIBRARY "/home/vgdgirard/local/lib/libhdf5.a" CACHE FILEPATH "" FORCE) 
SET(HDF5_FORTRAN_LIBRARY "/home/vgdgirard/local/lib/libhdf5_fortran.a" CACHE FILEPATH "" FORCE)

SET(FFTW_ENABLED ON CACHE BOOL "" FORCE)
SET(FFTW_INCLUDE_DIRS "/home/vgdgirard/local/include" CACHE PATH " " FORCE)
SET(FFTW_LIBRARY "/home/vgdgirard/local/lib/libfftw3.a" CACHE FILEPATH " " FORCE)
SET(FFTW_MPI_INCLUDE_DIR "/home/vgdgirard/local/include"  CACHE PATH " " FORCE)

SET(MPI_C_COMPILER "/home/vgdgirard/local/bin/mpicc" CACHE FILEPATH " " FORCE)
SET(MPI_CXX_COMPILER "/home/vgdgirard/local/bin/mpic++" CACHE FILEPATH " " FORCE)
SET(MPI_Fortran_COMPILER "/home/vgdgirard/local/bin/mpif90" CACHE FILEPATH " " FORCE)
SET(MPI_Fortran_INCLUDE_PATH "/home/vgdgirard/local/include" CACHE STRING " " FORCE)
SET(MPIEXEC "/home/vgdgirard/local/bin/mpirun" CACHE FILEPATH " " FORCE)

SET(CMAKE_EXE_LINKER_FLAGS "-Xlinker -rpath -Xlinker /home/vgdgirard/local/lib -L/home/vgdgirard/local/lib" CACHE STRING " " FORCE)
