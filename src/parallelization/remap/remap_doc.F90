! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory
! in file Doxyfile.in (line 691) if it is excluded.
! Type 'make doc' in build directory.
! To check the results, open :
! selalib/doc/build/html/doxygen/html/defgroup remaps.html
! The following lines will be read by doxygen to generate documentation:

!> @defgroup remap sll_remap
!> @brief
!> provides capabilities for global data reconfigurations in a parallel
!> machine.
!> @details
!> Suppose that a given dataset is represented by a multi-dimensional
!> array which is distributed on multiple processors. The specifics of this
!> distribution (which portion of the global array is contained in which
!> processor) are contained in an object called a 'layout'. To reconfigure
!> data, we need initial and final layouts, a remap 'plan' which uses the
!> given layouts for its initialization and finally an application of the
!> plan on the data described by the layouts. The data reconfiguration is an
!> out-of-place operation, so the module client is responsible for the
!> allocation of the appropriate arrays. Remap operates on multi-dimensional
!> arrays of several of the basic Fortran types.
!>
!> @author
!> Selalib team
!>
!>
!> <b> How to use it </b>
!> - Include : \code use sll_remapper \endcode
!> - Link with   <code>-lsll_remap</code>
!> - Add <code> use sll_remap </code>
