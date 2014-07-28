.. toctree::
   :maxdepth: 2

Install xcode and command-line tools (OSX Mavericks) ::

	$ xcode-select --install﻿

==============================
Install dependencies with Homebrew
==============================

.. note::
 Homebrew is the best way to install selalib dependencies on a mac.
 
Install Homebrew (`not compatible with macports <https://guide.macports.org/chunked/installing.macports.uninstalling.html>`_) ::

	$ ruby -e "$(curl -fsSL https://raw.github.com/mxcl/homebrew/go/install)"
	$ brew tap homebrew/science
	$ brew install git
	$ brew install cmake
	$ brew install gcc
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


.. warning::
 Since 2014 update of macports, some changes in openmpi could cause troubles.

Install macports (http://www.macports.org/install.php) and ::

	$ sudo port install gcc48 
	$ sudo port select gcc mp-gcc48
	$ sudo port install openmpi-default +threads
	$ sudo port select mpi openmpi-gcc48-fortran
	$ sudo port install hdf5-18 +fortran+gfortran+openmpi
	$ sudo port install fftw-3
	$ sudo port install cmake git-core 

If you want to run tests that use nurbs ::

	$ sudo port install python34
	$ sudo port select --set python python34

A special options file is available in selalib/cmake directory ::

	$ mkdir selalib/prototype/build
	$ cd selalib/prototype/build
	$ cmake ../src -DOPTIONS_FILE=../../cmake/macports_update_2014.cmake 
