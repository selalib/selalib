#!/bin/bash

# ---
# Example script to build selalib on MPCDF systems (Draco, Hydra, Linux clusters).
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
BUILD_TYPE="Release"
#BUILD_TYPE="Debug"

# select the build directory
BUILD_DIR=`pwd`

# add custom preprocessor macros (e.g. toggle 32 bit halos)
#MACROS="-DUSE_HALO_REAL32"

# number of processors to be used for a parallel build
JMAKE=16
# automatically determine the absolute location of the selalib source tree
SLL_BASE=$HOME/selalib


# --- machine-dependent modules and optimization flags
source /etc/profile.d/modules.sh
module purge
# ---
if [ x"$CLUSTER" == x"DRACO" ]; then
    module load intel/16.0
    module load mkl/11.3
    module load impi/5.1.3
    module load fftw
    module load hdf5-mpi
    module load cmake/3.5
    module load anaconda/3
    # --- optimization flags
    INTEL_OPT_FLAGS="-O3 -xHost -nowarn -qopt-report"
else
    if [ x"$CLUSTER" != x"" ]; then
      # on any Linux cluster, we have Intel MPI in the impi module
      module load intel/16.0
      module load mkl/11.3
      module load impi/5.1.3
      module load fftw
      module load hdf5-mpi
      module load cmake/3.2
      module load anaconda/3
      # --- optimization flags
      #INTEL_OPT_FLAGS="-O3 -xSSE4.2 -nowarn -qopt-report"
      INTEL_OPT_FLAGS="-O3 -xAVX -nowarn -qopt-report"
    else
      # whereas on Hydra, the module is labeled differently
      module load intel/16.0
      module load mkl/11.3
      module load mpi.intel/5.1.3
      module load fftw
      module load hdf5-mpi
      module load cmake/3.2
      module load anaconda/3
      # --- optimization flags
      #INTEL_OPT_FLAGS="-O3 -xSSE4.2 -nowarn -qopt-report"
      INTEL_OPT_FLAGS="-O3 -xAVX -nowarn -qopt-report"
    fi
fi
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

