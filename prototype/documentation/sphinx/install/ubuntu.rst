.. toctree::
   :maxdepth: 2

========================================
Install dependencies on UBUNTU or DEBIAN
========================================

The prerequisites can be installed with a couple of commands on Ubuntu. The only choice to make is between openmpi and mpich2. Most of our testing is done with openmpi but mpich2 should also work.

For openmpi do:

 $ sudo apt-get install libopenmpi-dev openmpi-bin libhdf5-openmpi-dev

For mpich2 do:

 $ sudo apt-get install mpich2 libmpich2-dev libhdf5-mpich2-dev

Please make sure not to install mpich2 and openmpi together. When both openmpi and mpich2 are installed strange errors will occur and selalib will not work. If you see both installed please remove both and install one.

Install source code manager :

 $ sudo apt-get install git

Install build tools :

 $ sudo apt-get install cmake

Install compilers :

 $ sudo apt-get install gfortran

Install libraries (MPI, HDF5, FFTW3, LAPACK) :

 $ sudo apt-get install  libopenmpi-dev libfftw3-dev \
    liblapack-dev 

Install software for documentation:

 $ sudo apt-get install doxygen texlive-latex3

For some outputs you can install also:

 $ sudo apt-get install gnuplot
