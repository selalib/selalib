Quick Start Guide
*****************

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

