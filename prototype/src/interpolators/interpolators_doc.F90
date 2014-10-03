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
!! @brief
!! Classes for numerical interpolation.
!> @author Selalib team 
!! @details
!! Methods :
!! - compute_interpolants
!! - interpolate_value
!! - interpolate_derivative_eta1
!! - interpolate_array
!! - interpolate_array_disp
!! - reconstruct_array
!! - interpolate_array_values
!! - interpolate_pointer_values
!! - interpolate_array_derivatives
!! - interpolate_pointer_derivatives
!! - set_coefficients
!! - get_coefficients
!!
!> Classes available are:
!!  - sll_interpolate_1d_base
!!  - sll_interpolate_2d_base
!!  - sll_cubic_spline_interpolator_1d
!!  - sll_cubic_spline_interpolator_1d_nonuniform
!!  - sll_cubic_spline_interpolator_2d
!!  - sll_arbitrary_degree_spline_interpolator_1d
!!  - sll_arbitrary_degree_spline_interpolator_2d
!!  - sll_periodic_interpolator_1d
!!  - sll_lagrange_interpolator_1d
!>
!> <b> Headers file available </b>
!>  - sll_interpolators.h
!>
!> <b> Modules available </b>
!>  List fortran module available
!!  - sll_module_interpolators_1d_base
!!  - sll_module_interpolators_2d_base
!!  - sll_module_cubic_spline_interpolator_1d
!!  - sll_module_cubic_spline_interpolator_1d_nonuniform
!!  - sll_module_cubic_spline_interpolator_2d
!!  - sll_module_arbitrary_degree_spline_interpolator_1d
!!  - sll_module_arbitrary_degree_spline_interpolator_2d
!!  - sll_module_periodic_interpolator_1d
!!  - sll_module_lagrange_interpolator_1d
!>
!> <b> How to use it </b>
!> - Header file : \code #include 'sll_interpolators.h' \endcode
!> - Link with   <code>-lsll_interpolators</code>
!> - Add <code> use sll_interpolators </code>
!>
!> <b> Examples </b>
!> @todo Add some fortran lines to explain how to use the library
!>
!> @test cubic_spline_interpolator_1d 
!> @test cubic_spline_interpolator_2d 
!> @test arbitrary_degree_spline_interpolator_1d 
!> @test arbitrary_degree_spline_interpolator_2d 
!> @test periodic_interpolator_1d
