#!/bin/bash
source /etc/profile.d/modules.sh
module purge
module load mpi/mpich-x86_64
rm -rf build
mkdir build
cd build
cmake ../ -DHDF5_PARALLEL_ENABLED=ON -DMPI_Fortran_INCLUDE_PATH=${MPI_FORTRAN_MOD_DIR} -DHDF5_INCLUDE_DIR=/usr/include/mpich-x86_64 -DHDF5_INCLUDE_DIR_FORTRAN=/usr/include/mpich-x86_64 -DHDF5_ROOT=/usr/lib64/mpich -DCMAKE_BUILD_TYPE=Release
make 
