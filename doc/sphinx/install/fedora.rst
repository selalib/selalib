.. toctree::
   :maxdepth: 2

==============================
Install dependencies on FEDORA
==============================

Install source code manager ::

 $ sudo dnf install git-all

Install build tools ::

 $ sudo dnf install cmake

Install compilers ::

 $ sudo dnf install gcc-gfortran

Install libraries (MPI, HDF5, FFTW3, LAPACK) ::

 $ sudo dnf install hdf5-openmpi-devel openmpi-devel fftw-devel

Install software for documentation::

 $ sudo dnf install doxygen texlive-latex

Don't forget to load the mpi module ::

 $ module load openmpi-x86_64

or::

 $ module load mpich2-x86_64

In this last case change openmpi by mpich2 for hdf5 and mpi packages
