Installation
------------
       
Change the variable CMAKE_BUILD_TYPE to "Release".

Set the variable CMAKE_INSTALL_PREFIX to the path where you want to install selalib.

Set the variable BUILD_PACKAGE to "ON" to reduce the building time.

Fortran modules, header files and library will be installed to this path.
Just type:
~~~
   mkdir build; cd build
   cmake -DCMAKE_BUILD_TYPE="Release" \
         -DCMAKE_INSTALL_PREFIX=<install_dir> \
         -DBUILD_PACKAGE=1 ../
   make 
   make install
~~~

Use Selalib
-----------

To use the sequential version of selalib just add:
~~~
   #include "selalib.h"
~~~
and link your program with flag *-lselalib*. Or to use the parallel version 
of selalib just add:
~~~
   #include "selalib-mpi.h"
~~~
and link your program with flag *-lselalib-mpi*. 

Be careful, you could have to link with HDF5, FFTW, LAPACK, etc ...

Il you want to use some macros like SLL_ALLOCATE or SLL_ASSERT with gfortran, 
just add the compilation flag *-ffree-line-length-none*.

Find this example in directory : package/examples

CMakeLists.txt
--------------

Find config files for selalib and fftpack in directory called "cmake":
```cmake
   PROJECT(Landau)
   ENABLE_LANGUAGE(Fortran)
   CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
   SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})
   FIND_PACKAGE(SeLaLib)
   FIND_PACKAGE(FFTPACK)
   MESSAGE(STATUS "SELALIB_INCLUDES:${SELALIB_INCLUDES}")
   MESSAGE(STATUS "SELALIB_LIBRARIES:${SELALIB_LIBRARIES}")
   INCLUDE_DIRECTORIES(${SELALIB_INCLUDES})
   ADD_EXECUTABLE(landau landau.F90)
   TARGET_LINK_LIBRARIES(landau ${SELALIB_LIBRARIES} 
                                ${FFTPACK_LIBRARIES})
```

Makefile
--------

SLL_ROOT is the path equal to the variable CMAKE_INSTALL_PREFIX:
```
   SLL_ROOT=/usr/local
   F90 = gfortran
   OPT = -O3
   
   F90FLAGS= -I${SLL_ROOT}/include/selalib -ffree-line-length-none
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
```
SConstruct
----------

If selalib is installed in /usr/local:
```python
   import os

   SLL_ROOT='/usr/local'

   env = Environment( ENV=os.environ,
                      LIBS=['selalib','dfftpack'],
                      F90='ifort',
                      F90FLAGS = ['-O3'],
                      F90PATH = [SLL_ROOT+'/usr/include/selalib'],
                      LINK='ifort',
                      LIBPATH = [SLL_ROOT+'/usr/lib'])

   env.Program('landau', ['landau.F90'])
```