!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> @ingroup interpolators
!> \brief abstract data type for 2d interpolation
!> @details
!> 
!> @todo
!> Function names should be reviewed and improved. What is the best way to
!> express that a derivative is in a particular direction? Why eta???
module sll_m_interpolators_2d_base
#include "sll_working_precision.h" 
use sll_m_boundary_condition_descriptors
implicit none
  
  !*************************************************************************
  !
  !                          2D Interpolators
  !
  !*************************************************************************
  
  
  !> Base class/basic interface for 2D interpolators
  type, abstract :: sll_c_interpolator_2d

   contains
     
     
     procedure(interpolator_two_arg_msg),  &
          deferred, pass(interpolator) :: interpolate_from_interpolant_value

     !>  Compute the value of the derivative of the interpolant with respect to eta1 of at a given abscissa \a x from the precomputed coefficients of the interpolator (stored in the interpolation object).
     procedure(interpolator_two_arg_msg),  &
          deferred, pass(interpolator) :: interpolate_from_interpolant_derivative_eta1

     !>  Compute the value of the derivative of the interpolant with respect to eta2 of at a given abscissa \a x from the precomputed coefficients of the interpolator (stored in the interpolation object).
     procedure(interpolator_two_arg_msg),  &
          deferred, pass(interpolator) :: interpolate_from_interpolant_derivative_eta2
     
     !> Compute the value of the interpolant at several given abscissae given as an array from given function values. Does not use a precomputed interpolant.
     procedure(interpolate_2d_array),      &
          pass, deferred :: interpolate_array
     
     !> Compute the value of the interpolant at all grid points shifted by the given displacement from function values. Does not use a precomputed interpolant.
     procedure(interpolate_2d_array_disp), &
          pass, deferred :: interpolate_array_disp
     

     !> Set value of coefficients of the interpolant to the given values.
     procedure(interpolator_2d_set_coeffs), &
          pass, deferred :: set_coefficients

     
     !> Check if interpolant was computed.
     procedure(interpolator_2d_logical_query), &
          pass, deferred :: coefficients_are_set
     
     !> Compute coefficients of the interpolants and stores it in the interpolation object.
     procedure(compute_coeffs_2d),&
          pass, deferred ::  compute_interpolants

     !> Extract the value of the precomputed coefficients of the interpolation from the interpolator object.
     procedure(get_coeffs_2d), &
          pass,deferred :: get_coefficients

     !> Delete the interpolator
     procedure(delete_interpolator_2d), & 
          pass, deferred :: delete
 
    ! generic, public :: delete => del !operator(delete) => del!
     
  end type sll_c_interpolator_2d
  


  !> Signature of interpolate_from_interpolant_value
  abstract interface
     function interpolator_two_arg_msg( interpolator, eta1, eta2 ) result(val)

       use sll_m_working_precision
       import sll_c_interpolator_2d
       sll_real64                                  :: val  !< interpolated value
       class(sll_c_interpolator_2d), intent(in) :: interpolator !< interpolator object
       sll_real64, intent(in)                      :: eta1 !< first coordinate of point where to interpolate
       sll_real64, intent(in)                      :: eta2 !< second coordinate of point where to interpolate

     end function interpolator_two_arg_msg

  end interface


  !> Compute interpolated values of n*m points
  abstract interface

     subroutine interpolate_2d_array(this,             &
          num_points1,      &
          num_points2,      &
          data_in,          &
          eta1,             &
          eta2,             &
          data_out)

       use sll_m_working_precision
       import sll_c_interpolator_2d
       class(sll_c_interpolator_2d), intent(in)    :: this !< interpolator object
       sll_int32,                       intent(in)    :: num_points1 !< number of points along first dimension
       sll_int32,                       intent(in)    :: num_points2 !< number of points along second dimension
       sll_real64,                      intent(in)    :: data_in(:,:) !< function values
       sll_real64,                      intent(in)    :: eta1(:,:)  !< values of the first coordinate for interpolation points
       sll_real64,                      intent(in)    :: eta2 (:,:) !< values of the second coordinate for interpolation points
       sll_real64,                      intent(out)   :: data_out(num_points1,num_points2) !< interpolated values

     end subroutine interpolate_2d_array

  end interface

  !> Signature of interpolate_array_disp
  abstract interface

     subroutine interpolate_2d_array_disp(this,        &
          num_points1, &
          num_points2, &
          data_in,     &
          alpha1,      &
          alpha2,      &
          data_out)

       use sll_m_working_precision
       import sll_c_interpolator_2d
       class(sll_c_interpolator_2d), intent(in)    :: this !< interpolator object
       sll_int32,                       intent(in)    :: num_points1  !< values of the first coordinate for interpolation points
       sll_int32,                       intent(in)    :: num_points2 !< values of the second coordinate for interpolation points
       sll_real64,                      intent(in)    :: data_in(:,:) !< function values
       sll_real64,                      intent(in)    :: alpha1(:,:) !< displacements along first dimension
       sll_real64,                      intent(in)    :: alpha2(:,:)  !< displacement along second dimesion
       sll_real64,                      intent(out)   :: data_out(num_points1,num_points2) !< interpolated values

     end subroutine interpolate_2d_array_disp

  end interface

  !> Signature of set_coefficients (Set the splines coefficients)
  abstract interface
     subroutine interpolator_2d_set_coeffs( &
          interpolator,&
          coeffs_1d,&
          coeffs_2d,&
          coeff2d_size1,&
          coeff2d_size2,&
          knots1,&
          size_knots1,&
          knots2,&
          size_knots2)
       use sll_m_working_precision
       import sll_c_interpolator_2d
       class(sll_c_interpolator_2d), intent(inout) :: interpolator !< interpolator object
       ! We allow the coefficients to be passed as 1d or 2d arrays. This allows
       ! for more flexibility for the children classes.
       sll_real64, dimension(:), intent(in), optional   :: coeffs_1d !< coefficients spezified in 1D array
       sll_real64, dimension(:,:), intent(in), optional :: coeffs_2d !< coefficients spezified as 2D array
       ! size coeffs 2D 
       sll_int32, intent(in), optional :: coeff2d_size1 !< size of 2D coeffs, first dimension
       sll_int32, intent(in), optional :: coeff2d_size2 !< size of 2D coeffs, second dimension
       sll_real64, dimension(:), intent(in), optional   :: knots1 !< knots first dimension
       sll_real64, dimension(:), intent(in), optional   :: knots2 !< knots second dimension
       sll_int32, intent(in), optional :: size_knots1 !< no. of knots first dimension
       sll_int32, intent(in), optional :: size_knots2 !< no. of knots second dimension
     end subroutine interpolator_2d_set_coeffs
  end interface

  !> Signature of coefficients_are_set (Check interpolator is computed)
  abstract interface
     function interpolator_2d_logical_query( interpolator ) result(res)
       import sll_c_interpolator_2d
       class(sll_c_interpolator_2d), intent(in) :: interpolator !< interpolator object
       logical :: res !< logical to specify if coefficients are set yes/no
     end function interpolator_2d_logical_query
  end interface

  !> Signature of compute_interpolants (Compute splines coefficients)
  abstract interface
     subroutine compute_coeffs_2d(interpolator, &
          data_array, &
          eta1_coords, &
          size_eta1_coords, &
          eta2_coords, &
          size_eta2_coords )
       use sll_m_working_precision
       import sll_c_interpolator_2d
       class(sll_c_interpolator_2d), intent(inout)  :: interpolator !< interpolator object
       sll_real64, dimension(:,:), intent(in)          :: data_array !< function values
       sll_real64, dimension(:), intent(in),optional   :: eta1_coords !< first coordinates of the grid points
       sll_real64, dimension(:), intent(in),optional   :: eta2_coords !< seoncd coordinates of the grid points
       sll_int32, intent(in), optional                 :: size_eta1_coords !< size of eta1_coords
       sll_int32, intent(in),optional                  :: size_eta2_coords !< size of eta2_coords
     end subroutine compute_coeffs_2d
  end interface
  
  !> Signature of get_coefficients (Get splines coefficients)
  abstract interface 
     function get_coeffs_2d(interpolator)
       use sll_m_working_precision
       import sll_c_interpolator_2d
       class(sll_c_interpolator_2d), intent(in) :: interpolator !< intepolator object
       sll_real64, dimension(:,:), pointer         :: get_coeffs_2d  !< value of the coefficients  
     end function get_coeffs_2d
  end interface

  !> Signature of delete (Deallocate the interpolator object)
  abstract interface 
     subroutine delete_interpolator_2d(interpolator)
       import sll_c_interpolator_2d
       class(sll_c_interpolator_2d), intent(inout) :: interpolator !< interpolator object
     end subroutine delete_interpolator_2d
  end interface


end module sll_m_interpolators_2d_base
