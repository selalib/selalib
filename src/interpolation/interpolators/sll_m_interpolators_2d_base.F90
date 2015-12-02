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
  type, abstract :: sll_interpolator_2d_base

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
     
  end type sll_interpolator_2d_base
  

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  !> Signature for interpolating function
  abstract interface
     function interpolator_two_arg_msg( interpolator, eta1, eta2 ) result(val)

       use sll_m_working_precision
       import sll_interpolator_2d_base
       sll_real64                                  :: val
       class(sll_interpolator_2d_base), intent(in) :: interpolator
       sll_real64, intent(in)                      :: eta1
       sll_real64, intent(in)                      :: eta2

     end function interpolator_two_arg_msg

  end interface

  !> Compute interpolated values of n*m points
  abstract interface

     subroutine interpolator_2d_array_msg( interpolator, data_array )

       use sll_m_working_precision
       import sll_interpolator_2d_base
       class(sll_interpolator_2d_base), intent(inout) :: interpolator
       sll_real64, dimension(:,:), intent(in)         :: data_array

     end subroutine interpolator_2d_array_msg

  end interface

  !> Compute interpolated values of n*m points
  abstract interface

     function interpolate_2d_array(this,             &
                                   num_points1,      &
                                   num_points2,      &
                                   data_in,          &
                                   eta1,             &
                                   eta2) result(res)

       use sll_m_working_precision
       import sll_interpolator_2d_base
       class(sll_interpolator_2d_base), intent(in)    :: this
       sll_int32, intent(in)                          :: num_points1 
       sll_int32, intent(in)                          :: num_points2 
       sll_real64, dimension(:,:), intent(in)         :: data_in
       sll_real64, dimension(:,:), intent(in)         :: eta1
       sll_real64, dimension(:,:), intent(in)         :: eta2  
       sll_real64, dimension(num_points1,num_points2) :: res

     end function interpolate_2d_array

  end interface

  !> Signature for interpolating function
  abstract interface

     function interpolate_2d_array_disp(this,        &
                                        num_points1, &
                                        num_points2, &
                                        data_in,     &
                                        alpha1,      &
                                        alpha2) result(res)

       use sll_m_working_precision
       import sll_interpolator_2d_base
       class(sll_interpolator_2d_base), intent(in)    :: this
       sll_int32, intent(in)                          :: num_points1  
       sll_int32, intent(in)                          :: num_points2 
       sll_real64, dimension(:,:), intent(in)         :: data_in
       sll_real64, dimension(:,:), intent(in)         :: alpha1
       sll_real64, dimension(:,:), intent(in)         :: alpha2  
       sll_real64, dimension(num_points1,num_points2) :: res

     end function interpolate_2d_array_disp

  end interface

  !> Set the splines coefficients
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
       import sll_interpolator_2d_base
       class(sll_interpolator_2d_base), intent(inout) :: interpolator
       ! We allow the coefficients to be passed as 1d or 2d arrays. This allows
       ! for more flexibility for the children classes.
       sll_real64, dimension(:), intent(in), optional   :: coeffs_1d
       sll_real64, dimension(:,:), intent(in), optional :: coeffs_2d
       ! size coeffs 2D 
       sll_int32, intent(in), optional :: coeff2d_size1
       sll_int32, intent(in), optional :: coeff2d_size2
       sll_real64, dimension(:), intent(in), optional   :: knots1
       sll_real64, dimension(:), intent(in), optional   :: knots2
       sll_int32, intent(in), optional :: size_knots1
       sll_int32, intent(in), optional :: size_knots2
     end subroutine interpolator_2d_set_coeffs
  end interface

  !> Check interpolator is computed
  abstract interface
     function interpolator_2d_logical_query( interpolator ) result(res)
       import sll_interpolator_2d_base
       class(sll_interpolator_2d_base), intent(in) :: interpolator
       logical :: res
     end function interpolator_2d_logical_query
  end interface

  !> Compute splines coefficients
  abstract interface
     subroutine compute_coeffs_2d(interpolator, &
          data_array, &
          eta1_coords, &
          size_eta1_coords, &
          eta2_coords, &
          size_eta2_coords )
       use sll_m_working_precision
       import sll_interpolator_2d_base
       class(sll_interpolator_2d_base), intent(inout)  :: interpolator
       sll_real64, dimension(:,:), intent(in)          :: data_array
       sll_real64, dimension(:), intent(in),optional   :: eta1_coords
       sll_real64, dimension(:), intent(in),optional   :: eta2_coords
       sll_int32, intent(in), optional                 :: size_eta1_coords
       sll_int32, intent(in),optional                  :: size_eta2_coords
     end subroutine compute_coeffs_2d
  end interface
  
  !> Get splines coefficients
  abstract interface 
     function get_coeffs_2d(interpolator)
       use sll_m_working_precision
       import sll_interpolator_2d_base
       class(sll_interpolator_2d_base), intent(in) :: interpolator
       sll_real64, dimension(:,:), pointer         :: get_coeffs_2d     
     end function get_coeffs_2d
  end interface

  !> Deallocate the interpolator object
  abstract interface 
     subroutine delete_interpolator_2d(interpolator)
       import sll_interpolator_2d_base
       class(sll_interpolator_2d_base), intent(inout) :: interpolator
     end subroutine delete_interpolator_2d
  end interface

#endif

end module sll_m_interpolators_2d_base
