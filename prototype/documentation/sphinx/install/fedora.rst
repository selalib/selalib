.. toctree::
   :maxdepth: 2

==============================
Install dependencies on FEDORA
==============================

Install source code manager ::

 $ sudo yum install git-all

Install build tools ::

 $ sudo yum install cmake

Install compilers ::

 $ sudo yum install gcc-gfortran

Install libraries (MPI, HDF5, FFTW3, LAPACK) ::

 $ sudo yum install hdf5-openmpi-devel openmpi-devel fftw-devel

Install software for documentation::

 $ sudo yum install doxygen texlive-latex

Don't forget to load the mpi module ::

 $ module load openmpi-x86_64

or::

 $ module load mpich2-x86_64
