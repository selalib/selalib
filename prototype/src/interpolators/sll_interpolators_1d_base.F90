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

!> \brief abstract data type for 1D interpolation and reconstruction
!>

module sll_module_interpolators_1d_base
#include "sll_working_precision.h"
  use sll_boundary_condition_descriptors
  implicit none

  type, abstract :: sll_interpolator_1d_base
   contains
     procedure(interpolator_1d_array_msg), deferred, pass(interpolator) :: &
          compute_interpolants
     procedure(interpolator_one_arg_msg), deferred, pass(interpolator) :: &
          interpolate_value
     procedure(interpolator_one_arg_msg), deferred, pass(interpolator) :: &
          interpolate_derivative_eta1
     procedure(interpolate_1d_array), pass, deferred :: interpolate_array
     procedure(interpolate_1d_array_at_displacement), pass, deferred :: interpolate_array_disp
     procedure(reconstruct_1d_array), pass, deferred :: reconstruct_array
     ! The following two are equivalent, and differ only by the type of
     ! the input and output data, one acts on 1d arrays, the other on 1d
     ! pointers. This is done for flexibility purposes.
     procedure(interpolator_1d_array_sub), deferred, pass(interpolator) :: &
          interpolate_array_values
     procedure(interpolator_1d_ptr_sub), deferred, pass(interpolator) :: &
          interpolate_pointer_values
     ! The following two are equivalent, and differ only by the type of
     ! the input and output data, oninterpolator_1d_array_sube acts on 1d arrays, the other on 1d
     ! pointers. This is done for flexibility purposes.
     procedure(interpolator_1d_array_sub), deferred, pass(interpolator) :: &
          interpolate_array_derivatives
     procedure(interpolator_1d_ptr_sub), deferred, pass(interpolator) :: &
          interpolate_pointer_derivatives
     ! Momentarily comment these out until we are actually going to have
     ! implementations for them.
     ! procedure(interpolate_1d_array), pass, deferred :: interpolate_array
     ! procedure(reconstruct_1d_array), pass, deferred :: reconstruct_array
      procedure(get_coeffs_1d), &
           pass,deferred :: get_coefficients
      procedure(interpolator_1d_set_coeffs), &
           pass, deferred :: set_coefficients
   end type sll_interpolator_1d_base

  sll_int32, parameter :: INTERP_PERIODIC_BC  = 0
  sll_int32, parameter :: INTERP_DIRICHLET_BC = 1
  sll_int32, parameter :: INTERP_NEUMANN_BC   = 2

 ! Signature of the interpolating function
  abstract interface
     function interpolator_one_arg_msg( interpolator, eta1 ) result(val)
       use sll_working_precision
       import :: sll_interpolator_1d_base
       sll_real64                                     :: val
       class(sll_interpolator_1d_base), intent(inout) :: interpolator
       sll_real64, intent(in)                         :: eta1
     end function interpolator_one_arg_msg
  end interface

  abstract interface
     subroutine interpolator_1d_array_msg( &
          interpolator, data_array,&
          eta_coords, &
          size_eta_coords)
       use sll_working_precision
       import :: sll_interpolator_1d_base
       class(sll_interpolator_1d_base), intent(inout) :: interpolator
       sll_real64, dimension(:), intent(in) :: data_array
       sll_real64, dimension(:), intent(in),optional  :: eta_coords
       sll_int32, intent(in),optional                 :: size_eta_coords
     end subroutine interpolator_1d_array_msg
  end interface

  abstract interface
     subroutine interpolator_1d_array_sub( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output_array )

       use sll_working_precision
       import :: sll_interpolator_1d_base
       class(sll_interpolator_1d_base), intent(in) :: interpolator
       sll_int32, intent(in)                :: num_pts
       sll_real64, dimension(:), intent(in) :: vals_to_interpolate
       sll_real64, dimension(:), intent(out):: output_array
     end subroutine interpolator_1d_array_sub
  end interface

  abstract interface
     subroutine interpolator_1d_ptr_sub( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output )

       use sll_working_precision
       import :: sll_interpolator_1d_base
       class(sll_interpolator_1d_base), intent(in) :: interpolator
       sll_int32, intent(in)                :: num_pts
       sll_real64, dimension(:), pointer :: vals_to_interpolate
       sll_real64, dimension(:), pointer :: output
     end subroutine interpolator_1d_ptr_sub
  end interface

  abstract interface
     function interpolate_1d_array(this, num_points, data, coordinates) &
          result(res)

       use sll_working_precision
       import sll_interpolator_1d_base
       class(sll_interpolator_1d_base), intent(in)     :: this
       sll_int32, intent(in)  :: num_points    ! size of output array
       sll_real64, dimension(:), intent(in) :: data  ! data to be interpolated
       ! points where output is desired
       sll_real64, dimension(:), intent(in) :: coordinates
       sll_real64, dimension(num_points)    :: res
     end function interpolate_1d_array
  end interface

  ! it is a bad practice to return large arrays like this. Must modify.
  abstract interface
     function interpolate_1d_array_at_displacement( &
       this, &
       num_points, &
       data, &
       alpha) result(res)

       use sll_working_precision
       import sll_interpolator_1d_base
       class(sll_interpolator_1d_base), intent(in)     :: this
       sll_int32, intent(in)  :: num_points    ! size of output array
       sll_real64, dimension(:), intent(in) :: data  ! data to be interpolated
       ! points where output is desired
       sll_real64, intent(in) :: alpha
       sll_real64, dimension(num_points)    :: res
     end function interpolate_1d_array_at_displacement
  end interface

  abstract interface
     function reconstruct_1d_array(this, num_points, data) result(res)
       use sll_working_precision
       import sll_interpolator_1d_base
       class(sll_interpolator_1d_base), intent(in)     :: this
       sll_int32, intent(in)     :: num_points    ! size of output array
       sll_real64, dimension(:), intent(in)  :: data ! data to be interpolated
       sll_real64, dimension(num_points)     :: res
     end function reconstruct_1d_array
  end interface

  abstract interface
     subroutine interpolator_1d_set_coeffs( interpolator, coeffs )
       use sll_working_precision
       import sll_interpolator_1d_base
       class(sll_interpolator_1d_base), intent(inout) :: interpolator
       ! We allow the coefficients to be passed as 1d or 2d arrays. This allows
       ! for more flexibility for the children classes.
       sll_real64, dimension(:), intent(in), optional   :: coeffs
     end subroutine interpolator_1d_set_coeffs
  end interface


   abstract interface
     function get_coeffs_1d(interpolator)
       use sll_working_precision
       import sll_interpolator_1d_base
       class(sll_interpolator_1d_base), intent(in) :: interpolator
       sll_real64, dimension(:), pointer         :: get_coeffs_1d
     end function get_coeffs_1d
  end interface

end module sll_module_interpolators_1d_base
