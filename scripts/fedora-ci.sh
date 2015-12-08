#!/bin/bash
source /etc/profile.d/modules.sh
module load mpi
rm -rf build
mkdir build
cd build; {
cmake -DCMAKE_BUILD_TYPE=Release -DHDF5_PARALLEL_ENABLED=ON ..
make
make Experimental
}; cd -
exit 0