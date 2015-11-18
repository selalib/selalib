REQUIRED COMPILER

The prototype is exploratory and forward-looking in character. We have decided
to test some additions to the Fortran language described in the Fortran 2003
standard. This standard is the default one applied when a user decides to use
gfortran, but not all of the features of the standard are supported. 
For example, the use of procedure pointers inside derived types is not 
recognized by gfortran in the versions previous to the 4.5 release. Hence, 
when building with gfortran, we presently require at least version 4.5.

This is a prototype, if this compiler requirement turned out to be too costly, 
we could re-implement the specific capabilities with external functions, but 
for the moment, we would like to see how far this approach can be taken.

BUILDING THE LIBRARY MODULES

Upon cloning the repository, you will see a collection of directories. Each
directory contains at least a library module, an SConstruct file and a unit 
test. You may go into any of these directories and run the 'scons' command
to locally build the library module and the corresponding unit test. To avoid
the need of building all of these modules manually, the top level of the 'src'
directory includes a Makefile whose sole function is to go to each library
directory in an appropriate order and run the 'scons' command. Effectively,
this tests the build of the existing modules. You can run either 'make clean'
or 'make' commands.

Problems in building any of the modules should be related with how to 
locate the right libraries for your system, as explained in the next section.

EXTERNAL LIBRARY DEPENDENCIES:

The prototype presently depends on:
- mpi
- hdf5
- lapack
- glibc: for high-resolution timer capabilities. This may pose a problem in
         systems that do not have a sufficiently recent version of glibc to
         use the high-resolution timer.

We want to offer a fine-grained control to the user in terms of the choice of
which version of any of these libraries to use, in case that the user can't or
does not wish to use the default options available in the system. 

To offer this fine-grained control, we are presently configured the build
process to use shell variables. Thus the user is responsible for going into
.bashrc (or .tcshrc or some other) and define a few variables. Presently:

HDF5_ROOT/include: directory where hdf5.h is located, 
HDF5_ROOT/lib: directory where libhdf5.a is located, 
HDF5_ROOT/include: directory where hdf5.mod is located

In .bashrc, this would be something like:
export HDF5_ROOT=/usr/local

In .tcshrc, this would be something like:
setenv HDF5_ROOT /usr/local

To develop in Selalib please read :
   - GitQuickstart.txt
   - CMakeQuickstart.txt

selalib compilation, installation
------------------------------------------

mkdir build 
cd build
cmake -DCMAKE_BUILD_TYPE=Release                                    \
      -DSLL_PACKAGE=1                                               \
      -DCMAKE_INSTALL_PREFIX=<path where selalib will be installed> \
      <the path of this directory>/src 
make install

selalib compilation, testing 
------------------------------------------

mkdir build
cd build
cmake <the path of this directory>
make Experimental
(the test result goes to http://cdash.inria.fr/CDash/index.php?project=Selalib)
