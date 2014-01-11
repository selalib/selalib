.. toctree::
   :maxdepth: 2

==============================
Install dependencies with Homebrew
==============================
 
Install Homebrew (not compatible with macports) ::

	$ ruby -e "$(curl -fsSL https://raw.github.com/mxcl/homebrew/go/install)"
	$ brew tap homebrew/science
      $ brew install git
      $ brew install cmake
	$ brew install hdf5 --enable-fortran --enable-parallel
	$ brew install pastix
	$ brew install fftw
	$ cd selalib/prototype/build
	$ cmake ../src -DPASTIX_ENABLED=ON  \
   		-DZLIB_LIBRARIES="/usr/lib/libz.dylib;/usr/local/lib/libsz.a" \
   		-DHDF5_PARALLEL_ENABLED=ON
	$ make

==============================
Install dependencies with macports
==============================

Install xcode and command-line tools ::

	$ xcode-select --install﻿
 
Install macports (http://www.macports.org/install.php) and ::

	$ sudo port install openmpi 
	$ sudo port install hdf5-18 +fortran+openmpi
	$ sudo port install fftw-3
	$ sudo port install cmake  +universal 
	$ sudo port install git-core +bash_completion
	$ sudo port install gnuplot doxygen texlive-latex


==============================
Install dependencies from scratch on MACOSX
==============================

Download gcc+gfortran on http://hpc.sourceforge.net ::

   $ sudo tar -zxvf gcc-bin.tar.gz -C /

Change PATH to use new compilers::

   $ export PATH="/usr/local/bin:${PATH}"
   $ export CC=/usr/local/bin/gcc
   $ export CXX=/usr/local/bin/g++
   $ export FC=/usr/local/bin/gfortran

Fortran must support 2003 norm (version > 4.6)::

   $ gfortran -v

Download and install OPENMPI::

   $ curl -L -O http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.3.tar.bz2
   $ tar jxvf openmpi-1.6.3.tar.bz2
   $ cd  openmpi-1.6.3
   $ ./configure --prefix=/usr/local --enable-mpi-thread-multiple
   $ make
   $ sudo make install
   $ cd -

Download and install HDF5::

   $ curl -L -O http://www.hdfgroup.org/ftp/lib-external/szip/2.1/src/szip-2.1.tar.gz
   $ tar jxvf szip-2.1.tar.gz
   $ cd szip-2.1
   $ ./configure --prefix=/usr/local 
   $ make
   $ sudo make install
   $ cd -
   $ curl -L -O http://zlib.net/zlib-1.2.7.tar.gz
   $ tar jxvf zlib-1.2.7.tar.gz
   $ cd zlib-1.2.7
   $ ./configure --prefix=/usr/local 
   $ make
   $ sudo make install
   $ cd -
   $ curl -L -O http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.10/src/hdf5-1.8.10.tar.bz2
   $ tar jxvf hdf5-1.8.10.tar.bz2
   $ cd hdf5-1.8.10
   $ export PATH=/usr/local/bin:${PATH}
   $ ./configure --enable-fortran --prefix=/usr/local --enable-parallel CC=mpicc FC=mpif90
   $ make
   $ sudo make install
   $ cd -

Download and install FFTW3::

   $ curl -L -O http://www.fftw.org/fftw-3.3.2.tar.gz
   $ tar zxvf fftw-3.3.2.tar.gz
   $ cd fftw-3.3.2
   $ ./configure --prefix=/usr/local --enable-shared --enable-threads --enable-mpi
   $ make
   $ sudo make install
   $ cd -

CMAKE, the most recent version is the best http://www.cmake.org/cmake/resources/software.html::

   $ tar zxvf cmake-2.8.10.1.tar.gz
   $ cd cmake-2.8.10.1
   $ ./bootstrap
   $ make
   $ sudo make install
   $ cd -

Install GIT
   http://code.google.com/p/git-osx-installer/
