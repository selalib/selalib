#!/bin/bash

# ---
# Example script to build selalib on MPCDF systems (Hydra, Linux clusters).
#
# By default, this script needs to be run from the ./scripts/ folder in the
# git checkout to be able to correctly determine the source location.  If you
# want to run it elsewhere, please set SLL_BASE accordingly (see below).
#
# 2016, khr@mpcdf.mpg.de
# ---


# --- OPTIONS ---
DO_CONF=true
DO_BUILD=true
BUILD_DIR=$HOME/selalib/obj
BUILD_TYPE="Release"
#BUILD_TYPE="Debug"
# number of processors to be used for a parallel build
JMAKE=1
# --- automatically determine the absolute location of the selalib source tree
SLL_BASE=`pwd`/../
SLL_BASE=`readlink -f $SLL_BASE`


# --- load environment modules
source /etc/profile.d/modules.sh
module purge
# ---
module load gcc
module load intel
module load mkl
if [ x"$CLUSTER" != x"" ]; then
  # on any Linux cluster, we have Intel MPI in the impi module
  module load impi
else
  # whereas on Hydra, the module is labeled differently
  module load mpi.intel
fi
module load fftw
module load hdf5-mpi
module load git
module load cmake/3.2
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
    -DCMAKE_BUILD_TYPE:STRING="${BUILD_TYPE}" \
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
