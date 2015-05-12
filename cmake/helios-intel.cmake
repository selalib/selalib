#once compiling in helios, use then
#cmake ../selalib/src -DOPTIONS_FILE=../cmake/helios-intel.cmake 

#bashrc should be similar to that
## .bashrc
## Source global definitions
#if [ -f /etc/bashrc ]; then
#        . /etc/bashrc
#fi
## User specific aliases and functions
#module load cmake
#module load intel
#module load intelmpi
#module load hdf5_p
#module load fftw
#export FC=ifort
#export CC=icc
#export CXX=icpc
#export ARCH=helios


#Configuration that should be (is) independent of the module versions
SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "We chose Release" FORCE)
SET(CMAKE_Fortran_COMPILER "ifort")

SET(HDF5_ENABLED ON CACHE BOOL "" FORCE)
SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "" FORCE)
SET(HDF5_ROOT $ENV{HDF5_DIR})

SET(MPI_Fortran_COMPILER "mpiifort")
SET(MPI_C_COMPILER "mpicc")
SET(MPI_CXX_COMPILER "mpiicpc")

SET(USE_MKL ON CACHE BOOL "" FORCE)
SET(CMAKE_Fortran_FLAGS_RELEASE "-nowarn -O3 -xHost -ip -fpic -g" CACHE STRING "Enable -g to analyse with vtune " FORCE)

# SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "We chose Release" FORCE)
# SET(CMAKE_Fortran_COMPILER "/opt/intel/composer_xe_2013_sp1.4.211/bin/intel64/ifort")
# 
# SET(FFTW_ENABLED ON CACHE BOOL "" FORCE)
# SET(FFTW_INCLUDE_DIRS "/csc/softs/fftw/fftw-3.3.4/intel-14.0.2.144/intelmpi-4.1.3.049/default/include")
# SET(FFTW_LIBRARY "/csc/softs/fftw/fftw-3.3.4/intel-14.0.2.144/intelmpi-4.1.3.049/default/lib/libfftw3.a")
# 
# 
# SET(HDF5_ENABLED ON CACHE BOOL "" FORCE)
# SET(HDF5_PARALLEL_ENABLED ON CACHE BOOL "" FORCE)
# SET(HDF5_ROOT "/csc/softs/hdf/hdf5-1.8.13_parallel/intel-14.0.2.144/intelmpi-4.1.3.049/default")
# SET(HDF5_C_LIBRARY "/csc/softs/hdf/hdf5-1.8.13_parallel/intel-14.0.2.144/intelmpi-4.1.3.049/default/lib/libhdf5.a")
# SET(HDF5_FORTRAN_LIBRARY "/csc/softs/hdf/hdf5-1.8.13_parallel/intel-14.0.2.144/intelmpi-4.1.3.049/default/lib/libhdf5_fortran.a")
# SET(HDF5_INCLUDE_DIR "/csc/softs/hdf/hdf5-1.8.13_parallel/intel-14.0.2.144/intelmpi-4.1.3.049/default/include")
# SET(HDF5_INCLUDE_DIR_FORTRAN "/csc/softs/hdf/hdf5-1.8.13_parallel/intel-14.0.2.144/intelmpi-4.1.3.049/default/include")
# 
# SET(MPI_Fortran_COMPILER "/opt/intel/impi/4.1.3.049/intel64/bin/mpiifort")
# SET(MPI_C_COMPILER "/opt/intel/impi/4.1.3.049/intel64/bin/mpicc")
# SET(MPI_CXX_COMPILER "/opt/intel/impi/4.1.3.049/intel64/bin/mpiicpc")

#intel/14.0.3.174          intel/15.0.0.090
#intel/14.0.4.211(default) intel/15.0.1.133
#bullxmpi/1.2.4.3            intelmpi/4.1.1.036          intelmpi/5.0.2.044
#bullxmpi/1.2.8.2(default)   intelmpi/4.1.3.049(default)
#bullxmpi_gnu/1.2.8.2        intelmpi/5.0.1.035

#old configuration 
#SET(CMAKE_BUILD_TYPE "Release")
#SET(CMAKE_Fortran_COMPILER "/opt/intel/composer_xe_2013.5.192/bin/intel64/ifort")
#SET(FFTW_INCLUDE_DIRS "/csc/softs/fftw/fftw-3.3.3/intel-13.1.3.192/intelmpi-4.1.1.036/default/include")
#SET(FFTW_LIBRARY "/csc/softs/fftw/fftw-3.3.3/intel-13.1.3.192/intelmpi-4.1.1.036/default/lib/libfftw3.a")
##SET(FFTW_MPI_INCLUDE_DIR "/u/system/SLES11/soft/fftw/3.3.2/intel-12.1/mpi.ibm-1.2/include")
##SET(FFTW_MPI_LIBRARY "/u/system/SLES11/soft/fftw/3.3.2/intel-12.1/mpi.ibm-1.2/lib/libfftw3_mpi.a")
##SET(FFTW_THREADS_LIBRARY "/u/system/SLES11/soft/fftw/3.3.2/intel-12.1/mpi.ibm-1.2/lib/libfftw3_threads.a")
#SET(HDF5_C_LIBRARY "/csc/softs/hdf/hdf5-1.8.9_parallel/intel-12.1.1.256/intelmpi-4.0.3/default/lib/libhdf5.a")
#SET(HDF5_FORTRAN_LIBRARY "/csc/softs/hdf/hdf5-1.8.9_parallel/intel-12.1.1.256/intelmpi-4.0.3/default/lib/libhdf5_fortran.a")
#SET(HDF5_INCLUDE_DIR "/csc/softs/hdf/hdf5-1.8.9_parallel/intel-12.1.1.256/intelmpi-4.0.3/default/include")
##SET(HDF5_INCLUDE_DIRS "/csc/softs/hdf/hdf5-1.8.9_parallel/intel-12.1.1.256/intelmpi-4.0.3/default/include")
#SET(HDF5_INCLUDE_DIR_FORTRAN "/csc/softs/hdf/hdf5-1.8.9_parallel/intel-12.1.1.256/intelmpi-4.0.3/default/include")
#SET(MPI_Fortran_COMPILER "/opt/intel/impi/4.1.1.036/intel64/bin/mpiifort")
