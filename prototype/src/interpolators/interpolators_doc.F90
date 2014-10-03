! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/prototype/documentation/build/html/doxygen/html/namespaces.html 
! The following lines will be read by doxygen to generate documentation:


!> @namespace sll_interpolators 
!! This library provides implementations for the abstract interfaces
!! Define spline interpolation of values in data define on original grid at
!! points coordinates
!> @author Selalib team 
!> @details
!> Abstract classes name are:
!!  - sll_interpolate_1d_base
!!  - sll_interpolate_2d_base
!>
!> <b> Headers file available </b>
!>  - sll_interpolators.h
!>
!> <b> Modules available </b>
!>  List fortran module available
!>  - sll_interpolators
!>
!> <b> How to use it </b>
!> - Header file : \code #include 'sll_interpolators.h' \endcode
!> - Link with   <code>-lsll_interpolators</code>
!> - Add <code> use sll_interpolators </code>
!>
!> <b> Examples </b>
!> -Add some fortran lines to explain how ti use the library
!> \code
!> call initialize(my_type, arg_1, arg_2)
!> call solve(my_type, your_result)
!> \endcode
!>
