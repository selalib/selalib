
=================
Getting Started
=================

To download the developer version you need to use GIT.
The latest development version of SeLaLib is tracked in the 'prototype-devel' branch.

Documentation for Git is available `here <http://git-scm.com/>`_.

First time set up
-----------------

Set up your developer info::

 git config --global user.name "Your NAME"
 git config --global user.email "you@somewhere"
 git config --global core.editor emacs (or other) 
 git config --global merge.tool opendiff (or other)
 
Check configuration with::

 git config --list

For help::

 git help <command-you-want-help-about>

Get the prototype
-----------------
Developer GIT Access via SSH

`Edit keys on your INRIA gforge account <https://gforge.inria.fr/account/editsshkeys.php>`_ and follow instructions.
Only project developers can access the GIT tree via this method. Substitute with the proper values::

 git clone git+ssh://YOUR_INRIA_LOGIN@scm.gforge.inria.fr//gitroot//selalib/selalib.git
 cd selalib/

Display all branches with::

 git branch -a

Create the local branch and merge the remote prototype branch:: 

 git checkout -b prototype-devel origin/prototype-devel

More information available on document `An overview of GIT A short tutorial for SELALIB developers <https://gforge.inria.fr/docman/view.php/3042/7642/selalib_coding_guidelines.pdf>`_ by Edwin Chacon-Golcher.

Build the prototype
-------------------

Create a new directory, build the libraries within ::

 cd prototype
 mkdir build
 cd build
 cmake ../src/
 make
 make test

CMake is a system for managing the build process of software. It writes in a CMakeCache.txt file all parameters that it has detected on your system. CMake can to not detect libraries installed for several reasons. So we need to edit the CMake cache.
If you have already performed cmake you can run ::

 make edit_cache

otherwise run ::

 ccmake ../src/

As you can see there is now a user interface that allows to edit parameters.
Press [c] to configure up to the configuration is valid and that the line "Press [g] to generate and exit" appear.

Details about principle parameters

+------------------------+----------------+--------------------------------------+
|          KEY           | POSSIBLE VALUE |                DESCRIPTION           |
+========================+================+======================================+
| MPI_MODULE_ENABLED     | ON/OFF         | When set to ON disable all reference |
|                        |                | to mpi and run SeLaLib in sequential.|
+------------------------+----------------+--------------------------------------+
|  FFT_LIB               | SLLFFT         | By default SLLFFT. SeLaLib provide   |
|                        | FFTPACK        | a fast fourier transform module      |
|                        | FFTW           | around 3 libraries, fftpack, fftw and|
|                        |                | his own implementation (see XXX)     |
+------------------------+----------------+--------------------------------------+
|  HDF5_ROOT             | /usr/local     | For output we need to link the lib   |
|  FFTW_ROOT             | /opt/local     | with HDF5 and FFTW3. Parallel is     |
|                        | /usr           | better but serial works. With this   |
|                        | ...            | you can use your own implementation  |
+------------------------+----------------+--------------------------------------+
