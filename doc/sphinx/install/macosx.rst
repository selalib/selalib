.. toctree::
   :maxdepth: 2

Install xcode and command-line tools (OSX Mavericks) ::

	$ xcode-select --install﻿

==================================
Install dependencies with Homebrew
==================================

.. note::
 Homebrew is the best way to install selalib dependencies on a mac.
 
Install Homebrew (`not compatible with macports <https://guide.macports.org/chunked/installing.macports.uninstalling.html>`_) ::

	$ ruby -e "$(curl -fsSL https://raw.github.com/mxcl/homebrew/go/install)"
	$ brew tap homebrew/science
	$ brew install git
	$ brew install cmake
	$ brew install gcc
	$ brew install openmpi --enable-mpi-thread-multiple
	$ brew install hdf5 --with-fortran --with-mpi
	$ brew install pastix
	$ brew install fftw
	$ cd selalib/build
	$ cmake ../ -DPASTIX_ENABLED=ON  \
		-DZLIB_LIBRARIES="/usr/lib/libz.dylib;/usr/local/lib/libsz.a" \
		-DHDF5_PARALLEL_ENABLED=ON
	$ make

==================================
Install dependencies with macports
==================================

.. warning::
 Don't forget to set aliases for openmpi and python with 'port select'

Install macports (http://www.macports.org/install.php) and ::

	$ sudo port install cmake git
	$ sudo port install openmpi
	$ sudo port install hdf5 +gfortran+openmpi
	$ sudo port install fftw-3
	$ sudo port install doxygen

If you want to run tests that use nurbs ::

	$ sudo port install python34
	$ sudo port select --set python python34

