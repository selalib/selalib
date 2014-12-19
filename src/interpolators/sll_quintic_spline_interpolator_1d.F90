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

#define sll_interpolator class(sll_quintic_spline_interpolator_1d)

!> @ingroup interpolators
!> @brief 
!! Interpolator 1d using quintic splines on regular mesh
!! @details
!! the following provides an implementation for the abstract interface
!! sll_interpolator_1d and define spline interpolation of values in 
!! data define on original grid at points coordinates
!!
module sll_module_quintic_spline_interpolator_1d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"

use sll_module_interpolators_1d_base
use sll_quintic_splines

implicit none
private

!> Quintic spline interpolator 1d
type, extends(sll_interpolator_1d_base), public :: sll_quintic_spline_interpolator_1d

  sll_real64, dimension(:),   pointer  :: x        !< points position
  sll_int32                            :: n        !< number of points
  sll_int32                            :: ind1     !< number of points
  sll_int32                            :: indn     !< number of points
  sll_real64, dimension(:,:), pointer  :: cf       !< values and derivatives
  sll_real64, dimension(:),   pointer  :: h        !< work array
  sll_int32                            :: bc_min   !< boundary condition
  sll_int32                            :: bc_max   !< boundary condition
  sll_real64                           :: xmin
  sll_real64                           :: xmax
  sll_real64                           :: slope_min
  sll_real64                           :: slope_max

contains

  !> PLEASE ADD DOCUMENTATION
  procedure, pass(interpolator) :: initialize => initialize_interpolator
  !> PLEASE ADD DOCUMENTATION
  procedure :: compute_interpolants => compute_interpolants_1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_value => interpolate_value_1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_derivative_eta1 => interpolate_deriv1_1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_array_values => interpolate_values_1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_pointer_values => interpolate_pointer_values_1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_array_derivatives => interpolate_derivatives_1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_pointer_derivatives => interpolate_pointer_derivatives_1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_array => spline_interpolate1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: interpolate_array_disp => spline_interpolate1d_disp
  !> PLEASE ADD DOCUMENTATION
  procedure :: reconstruct_array => reconstruct_array_1d 
  !> PLEASE ADD DOCUMENTATION
  procedure :: set_coefficients => set_coefficients_1d
  !> PLEASE ADD DOCUMENTATION
  procedure :: get_coefficients => get_coefficients_1d

end type sll_quintic_spline_interpolator_1d

!> Pointer to quintic spline interpolator implementation 1D
type :: sll_quintic_spline_interpolator_1d_ptr
  type(sll_quintic_spline_interpolator_1d), pointer :: interp
end type sll_quintic_spline_interpolator_1d_ptr

!> Deallocate the interpolator object
interface sll_delete
   module procedure delete_1d
end interface sll_delete

public new_quintic_spline_interpolator_1d
public sll_delete

contains  ! ****************************************************************


function new_quintic_spline_interpolator_1d( num_points, &
                                             xmin,       &
                                             xmax,       &
                                             bc_min,     &
                                             bc_max,     &
                                             slope_min,  &
                                             slope_max ) result(res)

  sll_interpolator,  pointer           :: res

  sll_int32,  intent(in)               :: num_points
  sll_real64, intent(in)               :: xmin
  sll_real64, intent(in)               :: xmax
  sll_int32,  intent(in)               :: bc_min
  sll_int32,  intent(in)               :: bc_max
  sll_real64, intent(in), optional     :: slope_min
  sll_real64, intent(in), optional     :: slope_max

  sll_int32                            :: ierr

  call initialize_1d_interpolator( res,        &
                                   num_points, &
                                   xmin,       &
                                   xmax,       &
                                   bc_min,     &
                                   bc_max,     &
                                   slope_min,  &
                                   slope_max )

end function new_quintic_spline_interpolator_1d

!---------------------------------------------------------------------------

subroutine initialize_interpolator( interpolator, &
                                    num_points,   &
                                    xmin,         &
                                    xmax,         &
                                    bc_min,       &
                                    bc_max,       &
                                    slope_min,    &
                                    slope_max )

  sll_interpolator,  intent(inout) :: interpolator

  sll_int32,  intent(in)           :: num_points
  sll_real64, intent(in)           :: xmin
  sll_real64, intent(in)           :: xmax
  sll_int32,  intent(in)           :: bc_min
  sll_int32,  intent(in)           :: bc_max
  sll_real64, intent(in), optional :: slope_min
  sll_real64, intent(in), optional :: slope_max
  sll_int32                        :: ierr
  sll_int32                        :: i
  sll_real64                       :: delta

  interpolator%n      = num_points
  interpolator%xmin   = xmin
  interpolator%xmax   = xmax
  interpolator%bc_min = bc_min
  interpolator%bc_max = bc_max

  if (present(slope_min) .and. present(slope_max)) then
    interpolator%slope_min = slope_min
    interpolator%slope_max = slope_max
  else
    interpolator%slope_min = 0.0_f64
    interpolator%slope_max = 0.0_f64
  end if

  SLL_ALLOCATE(interpolator%x(num_points),ierr)
  delta = (xmax - xmin) / (num_points - 1)
  do i = 1, num_points
     interpolator%x(i) = xmin + (i-1)*delta
  end do

  interpolator%ind1 = -1
  interpolator%indn = -1

  SLL_ALLOCATE(interpolator%h(6*num_points-3),ierr)

  SLL_ALLOCATE(interpolator%cf(3,num_points),ierr)

end subroutine

!---------------------------------------------------------------------------

function spline_interpolate1d(this,        &
                              num_points,  &
                              data,        &
                              coordinates) result(f_interp)

sll_interpolator,  intent(in)        :: this
sll_int32,  intent(in)               :: num_points
sll_real64, dimension(:), intent(in) :: coordinates
sll_real64, dimension(:), intent(in) :: data
sll_real64, dimension(num_points)    :: f_interp

sll_real64, allocatable :: c(:,:)
sll_real64, allocatable :: h(:)
sll_int32               :: ierr
sll_real64              :: f(3)
sll_int32               :: i

SLL_ALLOCATE(c(3,this%n),ierr)
SLL_ALLOCATE(h(6*this%n-3),ierr)

call inspl5(this%n, this%x, this%ind1, this%indn, c, h)

do i = 1, num_points
  call splin5(this%n, this%x, c, coordinates(i), f(1:3))
  f_interp(i) = f(1)
end do


end function

!---------------------------------------------------------------------------

function spline_interpolate1d_disp(this,        &
                                   num_points,  &
                                   data,        &
                                   alpha) result(f_interp)

sll_interpolator,                 intent(in) :: this

sll_real64, intent(in)                       :: alpha
sll_real64, dimension(:),         intent(in) :: data
sll_int32,                        intent(in) :: num_points
sll_real64, dimension(num_points)            :: f_interp
sll_real64, dimension(num_points)            :: coordinates

sll_real64                           :: delta
sll_real64                           :: length
sll_int32                            :: i
sll_real64                           :: x
sll_real64                           :: dx

dx = (this%xmax-this%xmin)/(num_points-1)
do i = 1, num_points
  coordinates(i) = this%xmin + (i-1)*dx
end do

do i = 1, num_points
  x = coordinates(i) - alpha
  SLL_ASSERT(x >= this%xmin)
  SLL_ASSERT(x <= this%xmax)
end do

SLL_ERROR("not implemented")

end function

!---------------------------------------------------------------------------

subroutine compute_interpolants_1d( interpolator, &
                                    data_array,   &
                                    eta_coords,   &
                                    size_eta_coords)

sll_interpolator, intent(inout)                :: interpolator
sll_real64, dimension(:), intent(in)           :: data_array
sll_real64, dimension(:), intent(in), optional :: eta_coords
sll_int32, intent(in),optional                 :: size_eta_coords

if(present(eta_coords) .or. present(size_eta_coords))then

  SLL_ERROR('Not implemented')

else

  call inspl5(interpolator%n,    &
              interpolator%x,    &
              interpolator%ind1, &
              interpolator%indn, &
              interpolator%cf,   &
              interpolator%h)

endif

end subroutine compute_interpolants_1d

!---------------------------------------------------------------------------

subroutine interpolate_values_1d( interpolator,        &
                                  num_pts,             &
                                  vals_to_interpolate, &
                                  output_array )

  sll_interpolator,  intent(in)  :: interpolator
  sll_int32,         intent(in)  :: num_pts
  sll_real64,        intent(in)  :: vals_to_interpolate(:)
  sll_real64,        intent(out) :: output_array(:)

  sll_int32                      :: i

  do i = 1, num_pts
    output_array(i) = interpolate_value_1d(interpolator,vals_to_interpolate(i))
  end do

end subroutine interpolate_values_1d

!---------------------------------------------------------------------------

subroutine interpolate_pointer_values_1d( interpolator,        &
                                          num_pts,             &
                                          vals_to_interpolate, &
                                          output )

  sll_interpolator,    intent(in)    :: interpolator
  sll_int32,           intent(in)    :: num_pts
  sll_real64, pointer                :: vals_to_interpolate(:)
  sll_real64, pointer                :: output(:)


end subroutine interpolate_pointer_values_1d

!---------------------------------------------------------------------------

subroutine interpolate_derivatives_1d( interpolator,        &
                                       num_pts,             &
                                       vals_to_interpolate, &
                                       output_array )

  sll_interpolator,         intent(in)  :: interpolator
  sll_int32,                intent(in)  :: num_pts
  sll_real64, dimension(:), intent(in)  :: vals_to_interpolate
  sll_real64, dimension(:), intent(out) :: output_array


end subroutine interpolate_derivatives_1d

!---------------------------------------------------------------------------

subroutine interpolate_pointer_derivatives_1d( interpolator,        &
                                               num_pts,             &
                                               vals_to_interpolate, &
                                               output )

  sll_interpolator,  intent(in)       :: interpolator
  sll_int32,  intent(in)              :: num_pts
  sll_real64, dimension(:), pointer   :: vals_to_interpolate
  sll_real64, dimension(:), pointer   :: output


end subroutine interpolate_pointer_derivatives_1d

!---------------------------------------------------------------------------

function interpolate_value_1d( interpolator, eta1 ) result(val)

  sll_interpolator, intent(in) :: interpolator
  sll_real64,       intent(in) :: eta1
  sll_real64                   :: val

  sll_real64                   :: f(3)

  call splin5(interpolator%n,  &
              interpolator%x,  &
              interpolator%cf, &
              eta1,            &
              f(1:3))

  val = f(1)

end function

!---------------------------------------------------------------------------

function interpolate_deriv1_1d( interpolator, eta1 ) result(val)
  sll_interpolator, intent(in) :: interpolator
  sll_real64             :: val
  sll_real64, intent(in) :: eta1

  
end function interpolate_deriv1_1d


!---------------------------------------------------------------------------

function reconstruct_array_1d(this, num_points, data) result(res)

  sll_interpolator, intent(in) :: this
  sll_int32, intent(in)                :: num_points
  sll_real64, dimension(:), intent(in) :: data
  sll_real64, dimension(num_points)    :: res

  res(:) = 0.0_f64

  SLL_WARNING('#warning reconstruct_array is dummy')

end function reconstruct_array_1d

!---------------------------------------------------------------------------

subroutine delete_1d( this )

  sll_interpolator :: this


end subroutine delete_1d

!---------------------------------------------------------------------------

subroutine set_coefficients_1d( interpolator, coeffs )

  sll_interpolator,  intent(inout) :: interpolator
  sll_real64, dimension(:), intent(in), optional :: coeffs

  if (present(coeffs)) then
    print*, interpolator%n
  end if

  SLL_ERROR('set_coefficients_1d() not implemented yet.')

end subroutine set_coefficients_1d

!---------------------------------------------------------------------------

function get_coefficients_1d(interpolator)

  sll_interpolator,  intent(in)      :: interpolator
  sll_real64, dimension(:), pointer  :: get_coefficients_1d

  get_coefficients_1d = interpolator%x
  SLL_ERROR('get_coefficients_1d() not been implemented yet.')

end function get_coefficients_1d

end module sll_module_quintic_spline_interpolator_1d
