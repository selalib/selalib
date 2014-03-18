.. toctree::
   :maxdepth: 2

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

.. warning::
 Since 2014 update of macports, some changes in openmpi could cause troubles.

Install macports (http://www.macports.org/install.php) and ::

	$ sudo port install gcc48 
	$ sudo port install hdf5-18 +fortran+gfortran+openmpi_devel
	$ sudo port install fftw-3
	$ sudo port install cmake git-core 
	$ mkdir selalib/prototype/build
	$ cd selalib/prototype/build
	$ cmake ../src -DOPTIONS_FILE=../../cmake/macports_update_2014.cmake 
