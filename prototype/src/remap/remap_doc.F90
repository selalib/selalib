! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/prototype/documentation/build/html/doxygen/html/namespaces.html 
! The following lines will be read by doxygen to generate documentation:


!> @namespace sll_remap 
!> @brief 
!> Remap provides capabilities for global data reconfigurations in a parallel 
!> machine. 
!> @details
!> All that is needed is a specification of the initial and final
!> configurations, the creation of a remap plan that uses these specifications
!> as input and an application of the plan. For details about the interface,
!> follow the link sll_remapper .
!> @details Suppose that a given dataset is represented by a multi-dimensional
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
!> <b> Headers file available </b>
!>  - sll_remap.h
!>
!> <b> Modules available </b>
!>  List fortran module available
!>  - sll_remap
!>
!> <b> How to use it </b>
!> - Header file : \code #include 'sll_remap.h' \endcode
!> - Link with   <code>-lsll_%s</code>
!> - Add <code> use sll_remap </code>
!>
!> <b> Examples </b>
!> -Add some fortran lines to explain how ti use the library
!> \code
!> call initialize(my_type, arg_1, arg_2)
!> call solve(my_type, your_result)
!> \endcode
!>
