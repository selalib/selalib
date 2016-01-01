REQUIRED COMPILER
-----------------

We use some additions to the Fortran language described in the Fortran 2003
standard. This standard is the default one applied when a user decides to use
gfortran, but not all of the features of the standard are supported. 
For example, the use of allocatable string is not 
recognized by gfortran in the versions previous to the 4.8 release. Hence, 
when building with gfortran, we presently require at least version 4.8.
For Intel Fortran compiler we need at least version 13.

BUILDING THE LIBRARY MODULES
----------------------------

Upon cloning the repository, you will see a collection of directories. Each
directory contains at least a library module, an CMakeLists.txt file and 
testing directory. To configure, create a build directory and use cmake 
command:
~~~~
mkdir build 
cd build
cmake -DCMAKE_BUILD_TYPE=Release                                    \
      -DSLL_PACKAGE=1                                               \
      -DCMAKE_INSTALL_PREFIX=<path where selalib will be installed> \
      <the path of the selalib directory>
make install
~~~~

Problems in building any of the modules should be related with how to 
locate the right libraries for your system, as explained in the next section.

For cmake configuration see [CMake Quickstart](CMakeQuickstart.md).

EXTERNAL LIBRARY DEPENDENCIES
-----------------------------

The prototype presently depends on:
  - mpi
  - hdf5
  - lapack

We want to offer a fine-grained control to the user in terms of the choice of
which version of any of these libraries to use, in case that the user can't or
does not wish to use the default options available in the system. 

To offer this fine-grained control, we are presently configured the build
process to use shell variables. Thus the user is responsible for going into
.bashrc (or .tcshrc or some other) and define a few variables. Presently:

- HDF5_ROOT/include: directory where hdf5.h is located, 
- HDF5_ROOT/lib: directory where libhdf5.a is located, 
- HDF5_ROOT/include: directory where hdf5.mod is located

In .bashrc, this would be something like:

  export HDF5_ROOT=/usr/local

In .tcshrc, this would be something like:

  setenv HDF5_ROOT /usr/local

To develop in Selalib please read :
   - GitQuickstart.md
   - CMakeQuickstart.md
   - CONTRIBUTING.md

Selalib compilation, testing 
----------------------------

~~~~
mkdir build
cd build
cmake <the path of this directory>
make Experimental
~~~~
(the test result goes to http://cdash.inria.fr/CDash/index.php?project=Selalib)