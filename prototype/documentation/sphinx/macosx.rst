.. toctree::
   :maxdepth: 2

==============================
Install dependencies on MACOSX
==============================

Download gcc+gfortran on http://hpc.sourceforge.net ::

   $ sudo tar -zxvf gcc-bin.tar.gz -C /

Change PATH to use new compilers::

   $ export PATH="/usr/local/bin:${PATH}"

Fortran must support 2003 norm (version > 4.6)::

   $ gfortran -v

Download and install OPENMPI::

   $ curl -L -O http://www.open-mpi.org/software/ompi/v1.4/downloads/openmpi-1.4.5.tar.bz2
   $ tar jxvf openmpi-1.4.5.tar.bz2
   $ cd  openmpi-1.4.5
   $ ./configure --prefix=/usr/local
   $ make
   $ sudo make install

Download and install HDF5::

   $ curl -L -O http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.9/src/hdf5-1.8.9.tar.bz2
   $ tar jxvf hdf5-1.8.9.tar.bz2
   $ cd hdf5-1.8.9
   $ export PATH=/usr/local/bin:${PATH}
   $ ./configure --enable-fortran --prefix=/usr/local
   $ make
   $ sudo make install

Download and install FFTW3::

   $ curl -L -O http://www.fftw.org/fftw-3.3.1.tar.gz
   $ tar zxvf fftw-3.3.1.tar.gz
   $ cd fftw-3.3.1
   $ ./configure --prefix=/usr/local --enable-shared --enable-threads
   $ make
   $ sudo make install

