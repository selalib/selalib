#!/bin/bash
source /etc/profile.d/modules.sh
module purge
module load mpi/openmpi-x86_64
rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DHDF5_PARALLEL_ENABLED=ON ..
make 
