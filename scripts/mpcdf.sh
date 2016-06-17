#!/bin/bash

# ---
# Example script to build selalib on MPCDF systems (Hydra, Linux clusters).
#
# Copy this script to your build folder and run it to compile selalib.
# Note that you need to reset SLL_BASE if your source folder has a different 
# name.
#
# 2016, khr@mpcdf.mpg.de
# ---


# --- OPTIONS ---
DO_CONF=true
DO_BUILD=true
BUILD_DIR=`pwd`
BUILD_TYPE="Release"
#BUILD_TYPE="Debug"
# --- optimization flags that override the CMAKE Fortran Release defaults ---
INTEL_OPT_FLAGS="-O3 -ip -fpic -xSSE4.2 -nowarn"
# --- add custom preprocessor macros (here: toggle 32 bit halos)
#MACROS="-DUSE_HALO_REAL32"
# number of processors to be used for a parallel build
JMAKE=16
# --- automatically determine the absolute location of the selalib source tree
SLL_BASE=$HOME/selalib
#SLL_BASE=`readlink -f $SLL_BASE`


# --- load environment modules
source /etc/profile.d/modules.sh
module purge
# ---
module load intel/15.0
module load mkl/11.3
if [ x"$CLUSTER" != x"" ]; then
  # on any Linux cluster, we have Intel MPI in the impi module
  module load impi/5.0.3
else
  # whereas on Hydra, the module is labeled differently
  module load mpi.intel/5.0.3
fi
module load fftw
module load hdf5-mpi
module load cmake/3.2
module load python33/python
# ---
module list


# --- set enviroment variables to actually use the Intel toolchain
export FC=ifort
export CC=icc
export CXX=icpc
#export I_MPI_F90=$FC
#export I_MPI_CC=$CC
#export I_MPI_CXX=$CXX


mkdir -p $BUILD_DIR
cd $BUILD_DIR

if [ x"$DO_CONF" == x"true" ]
then
  cmake $SLL_BASE \
    -DOPENMP_ENABLED:BOOL=ON \
    -DUSE_MKL:BOOL=ON \
    -DHDF5_PARALLEL_ENABLED:BOOL=ON \
    -DCMAKE_BUILD_TYPE:STRING="${BUILD_TYPE}" \
    -DMPI_C_COMPILER:STRING="mpiicc" \
    -DMPI_CXX_COMPILER:STRING="mpiicpc" \
    -DMPI_Fortran_COMPILER:STRING="mpiifort" \
    -DFORCE_Fortran_FLAGS_RELEASE:STRING="$INTEL_OPT_FLAGS" \
    -DCMAKE_C_FLAGS:STRING="-g $MACROS $INTEL_OPT_FLAGS" \
    -DCMAKE_CXX_FLAGS:STRING="-g $MACROS $INTEL_OPT_FLAGS" \
    -DCMAKE_Fortran_FLAGS:STRING="-g $MACROS"
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

