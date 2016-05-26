#!/bin/bash
source /opt/intel/bin/compilervars.sh intel64 -platform linux
source /opt/intel/mkl/bin/mklvars.sh intel64
export FC=ifort
export CC=icc
export CXX=icpc
export I_MPI_F90=ifort
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export HDF5_ROOT=/usr/local
rm -rf build
mkdir build
cd build; {
cmake -DCMAKE_BUILD_TYPE=Release -DHDF5_PARALLEL_ENABLED=ON ..
make
}; cd -
exit 0
