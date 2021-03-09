# SeLaLib

## Required compiler

We use some additions to the Fortran language described in the Fortran 2003
standard. 

## Building the library modules

Upon cloning the repository, you will see a collection of directories. Each
directory contains at least a library module, an CMakeLists.txt file and
testing directory. To configure, create a build directory and use cmake
command:

~~~~
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release  .. 
make 
~~~~

This procedure is tested only on ubuntu 20.04 with following dependencies:

- gfortran g++ cmake libopenmpi-dev openmpi-bin libhdf5-openmpi-dev libfftw3-dev liblapack-dev libopenblas-dev

For documentation install

- doxygen texlive graphviz
