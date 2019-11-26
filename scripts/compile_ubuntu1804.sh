#!/bin/bash

# ---
# Build selalib locally on an Ubuntu 16.04 machine.
# Run this script directly from within the "scripts" folder.
# ---


BUILD_DIR=${BUILD_DIR:=${HOME}/selalib_obj}
# run the cmake configure step
DO_CONF=${DO_CONF:=1}
# run the actual build
DO_BUILD=${DO_BUILD:=1}
# toggle openmp usage
USE_OPENMP=${USE_OPENMP:=ON}
# build the external package "fmempool"
USE_FMEMPOOL=${USE_FMEMPOOL:=OFF}
# number of processors to be used for a parallel build
JMAKE=${JMAKE:=1}
# define the CMAKE build type
CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:="Release"}

#export ZFP_ROOT=/opt/apps/zfp/current


SLL_BASE=`readlink -f ..`

mkdir -p $BUILD_DIR
cd $BUILD_DIR

export CC=gcc
export CXX=g++
export FC=gfortran


# HDF5, OpenMPI, FFTW, OpenBLAS can be installed using the package manager:
# $ apt-get install libopenmpi-dev libhdf5-openmpi-dev libhdf5-dev libfftw3-dev libopenblas-dev openmpi-bin


if [ x"$DO_CONF" == x"1" ]
then
  cmake $SLL_BASE \
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} \
    -DOPENMP_ENABLED:BOOL=${USE_OPENMP} \
    -DUSE_MKL:BOOL=OFF \
    -DFFTW_ENABLED:BOOL=ON \
    -DHDF5_PARALLEL_ENABLED:BOOL=ON \
    -DMPI_C_COMPILER:STRING="mpicc" \
    -DMPI_CXX_COMPILER:STRING="mpic++" \
    -DMPI_Fortran_COMPILER:STRING="mpif90" \
    -DCMAKE_C_FLAGS:STRING="-g" \
    -DCMAKE_CXX_FLAGS:STRING="-g" \
    -DCMAKE_Fortran_FLAGS:STRING="-g" \
    -DUSE_FMEMPOOL=${USE_FMEMPOOL}
else
  echo "skipping configure step ..."
fi

if [ x"$DO_BUILD" == x"1" ]
then
  make -j${JMAKE} VERBOSE=0
fi
