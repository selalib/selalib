#!/bin/bash
source /etc/profile.d/modules.sh
module purge
module load mpi/mpich-x86_64
rm -rf build
mkdir build
cd build; {
cmake -DHDF5_INCLUDE_DIR=/usr/include/mpich-x86_64 -DCMAKE_BUILD_TYPE=Release -DHDF5_PARALLEL_ENABLED=ON ..
make
make Experimental
}; cd -
exit 0
