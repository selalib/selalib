#!/bin/bash

# ---
# Build selalib on MPCDF systems (Hydra, Linux clusters) using the Intel toolchain.
# 2015, khr@mpcdf.mpg.de
# ---


# --- Determine the absolute location of the selalib source tree,
#     assuming that this script is invoked from the 'scripts' directory.
SLL_BASE=`pwd`/../
SLL_BASE=`readlink -f $SLL_BASE`
#SLL_BASE=/home/...
# --- build (object) directory
BUILD_DIR=$HOME/selalib/obj
DO_CONF=true
DO_BUILD=true
# --- number of processors to be used for a parallel build
JMAKE=1


# --- load environment modules
module purge
# ---
module load gcc
module load intel
module load mkl
module load impi || module load mpi.intel
module load hdf5-mpi
module load fftw
module load cmake/3.2
module load git
module load python33/python
# ---
module list


# --- set enviroment variables to actually use the Intel toolchain
export FC=ifort
export CC=icc
export CXX=icpc
export I_MPI_F90=$FC
export I_MPI_CC=$CC
export I_MPI_CXX=$CXX


mkdir -p $BUILD_DIR
cd $BUILD_DIR

if [ x"$DO_CONF" == x"true" ]
then
  cmake $SLL_BASE \
    -DCMAKE_BUILD_TYPE:STRING="Release" \
    -DOPENMP_ENABLED:BOOL=ON \
    -DUSE_MKL:BOOL=ON \
    -DHDF5_PARALLEL_ENABLED:BOOL=ON \
    -DCMAKE_CXX_FLAGS:STRING="-g" \
    -DCMAKE_C_FLAGS:STRING="-g" \
    -DCMAKE_Fortran_FLAGS:STRING="-g" \
    -DMPI_CXX_COMPILE_FLAGS:STRING="-g" \
    -DMPI_C_COMPILE_FLAGS:STRING="-g" \
    -DMPI_Fortran_COMPILE_FLAGS:STRING="-g"
else
  echo "skipping configure step ..."
fi

if [ x"$DO_BUILD" == x"true" ]
then
  gmake -j${JMAKE} VERBOSE=1
  #gmake Experimental VERBOSE=1
else
  echo "skipping make step ..."
fi

