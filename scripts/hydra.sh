#!/bin/bash
cd /ptmp/pin/all
module purge
module load intel
module load mkl
module load mpi.ibm
module load cmake
module load hdf5-mpi
module load fftw
module load python33/python
module load git
mkdir /ptmp/$USER/build_on_hydra
cd /ptmp/$USER/build_on_hydra
cmake $HOME/selalib/src -DUSE_MKL=ON -DHDF5_PARALLEL_ENABLED=ON
make -j8
make Experimental
