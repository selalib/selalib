! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/doc/build/html/doxygen/html/defgroup interpolatorss.html 
! The following lines will be read by doxygen to generate documentation:


!> @defgroup interpolators sll_interpolators 
!! @brief
!! Classes for numerical interpolation.
!! @authors Edwin-Chacon Golcher, Aurore Back and Pierre Navaro.
!! @details
!! PLEASE ADD DOCUMENTATION 
!!
!! Methods:
!! - compute_interpolants(...)
!! - interpolate_from_interpolant_value(...)
!! - interpolate_from_interpolant_derivative_eta1(...)
!! - interpolate_array(...)
!! - interpolate_array_disp(...)
!! - reconstruct_array(...)
!! - interpolate_from_interpolant_array(...)
!! - interpolate_pointer_values(...)
!! - interpolate_array_derivatives(...)
!! - interpolate_pointer_derivatives(...)
!! - set_coefficients(...)
!! - get_coefficients(...)
!!
!> Classes available are:
!!  - sll_m_interpolators_1d_base::sll_c_interpolator_1d
!!  - sll_m_interpolators_2d_base::sll_c_interpolator_2d
!!  - sll_m_cubic_spline_interpolator_1d::sll_cubic_spline_interpolator_1d
!!  - sll_m_cubic_spline_interpolator_1d_nonuniform::sll_cubic_spline_interpolator_1d_nonuniform
!!  - sll_m_cubic_spline_interpolator_2d::sll_cubic_spline_interpolator_2d
!!  - sll_m_arbitrary_degree_spline_interpolator_1d::sll_arbitrary_degree_spline_interpolator_1d
!!  - sll_m_arbitrary_degree_spline_interpolator_2d::sll_arbitrary_degree_spline_interpolator_2d
!!  - sll_periodic_interpolator_1d
!!  - sll_lagrange_interpolator_1d
!>
!> <b> How to use it </b>
!> - Include \code use <name_of_used_modules> \endcode in your file.
!> - Link with   <code>-lsll_interpolators</code>
!>
