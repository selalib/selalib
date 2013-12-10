Installation
************

Linux packages
==============

Selalib is available as a package 

* as a `deb package </releases/selalib_1.0.0_amd64.deb>`_
* as a `rpm package </releases/selalib-1.0.0.fc18.x86_64.rpm>`_

Install the package on Ubuntu 12.10 with::

    sudo dpkg -i selalib_1.0.0_amd64.deb

or on fedora core 18 ::

    sudo rpm -i selalib-1.0.0.fc18.x86_64.rpm

Use the library::

    gfortran -I/usr/lib/selalib/include landau.F90 -L/usr/lib/selalib/lib -lselalib -ldfftpack

    ./a.out | gnuplot


From git repository
===================
       
Change the variable CMAKE_BUILD_TYPE to "Release".

Set the variable CMAKE_INSTALL_PREFIX to the path where you want to install selalib.

Fortran modules, header files and library will be installed to this path.

* header files in ${CMAKE_INSTALL_PREFIX}/include
* fortran module files in ${CMAKE_INSTALL_PREFIX}/include/fortran
* library archives in ${CMAKE_INSTALL_PREFIX}/lib

Just type::

   cmake -DCMAKE_BUILD_TYPE="Release" \
         -DCMAKE_INSTALL_PREFIX=<install_dir> \
         -DSLL_BUILD_PACKAGE=ON <prototype_src_dir>
   make 
   make install


Use selalib
===========

To use the sequential version of selalib just add::

   #include "selalib.h"

To use the parallel version of selalib just add::

   #include "selalib-mpi.h"

and link your program with flag *-lselalib*. 
You can also use one of selalib capabilities separately, headers file available are::

   #include "sll_working_precision.h"
   #include "sll_memory.h"
   #include "sll_assert.h"
   #include "sll_splines.h"
   #include "sll_constants.h"
   #include "sll_utilities.h"
   #include "sll_interpolators.h"

To use solvers::

   #include "sll_poisson_solvers.h"
   #include "sll_maxwell_solvers.h"

Be careful, you could have to link with HDF5, FFTW, LAPACK, etc ...

Il you want to use some macros like SLL_ALLOCATE or SLL_ASSERT with gfortran, just add
the compilation flag *-ffree-line-length-none*.

Find this example in directory : examples


CMakeLists.txt
==============

Find config files for selalib and fftpack in directory called "cmake"::

   PROJECT(Landau)
   ENABLE_LANGUAGE(Fortran)
   CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
   SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})
   FIND_PACKAGE(SeLaLib)
   FIND_PACKAGE(FFTPACK)
   MESSAGE(STATUS "SELALIB_INCLUDES:${SELALIB_INCLUDES}")
   MESSAGE(STATUS "SELALIB_LIBRARIES:${SELALIB_LIBRARIES}")
   INCLUDE_DIRECTORIES(${SELALIB_INCLUDES})
   INCLUDE_DIRECTORIES(${SELALIB_INCLUDES}/fortran)
   ADD_EXECUTABLE(landau landau.F90)
   TARGET_LINK_LIBRARIES(landau ${SELALIB_LIBRARIES} 
                                ${FFTPACK_LIBRARIES})

Makefile
========

SLL_ROOT is the path equal to the variable CMAKE_INSTALL_PREFIX::

   SLL_ROOT=/usr/local
   F90 = gfortran
   OPT = -O3
   
   F90FLAGS= -I${SLL_ROOT}/include -J${SLL_ROOT}/include/fortran
   LDFLAGS=-L${SLL_ROOT}/lib 
   LIBS= -lselalib -ldfftpack

   PROG= landau_1d

   all: $(PROG)

   $(PROG): $(OBJS)
       $(LD) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
 
   clean:
       rm -f $(PROG) $(OBJS) *.mod

   .SUFFIXES: $(SUFFIXES) .F90

   .F90.o:
       $(F90) $(F90FLAGS) -c $<

   .mod.o:
       $(F90) $(F90FLAGS) -c $*.F90

SConstruct
==========

If selalib is installed in /usr/local::

   import os

   SLL_ROOT='/usr/local'

   env = Environment( ENV=os.environ,
                      LIBS=['selalib','dfftpack'],
                      F90='ifort',
                      F90FLAGS = ['-O3'],
	              F90PATH = [SLL_ROOT+'include',SLL_ROOT+'include/fortran'],
                      LINK='ifort',
	              LIBPATH = [SLL_ROOT+'/lib'])

   env.Program('landau', ['landau.F90'])

