#!/bin/bash
rm -rf build
mkdir build
cd build; {
cmake -DCMAKE_BUILD_TYPE=Release -DHDF5_PARALLEL_ENABLED=ON ..
make
ctest
}; cd -
exit 0
