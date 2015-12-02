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
module sll_m_quintic_spline_interpolator_1d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_errors.h"

use sll_m_interpolators_1d_base
use sll_m_quintic_splines

implicit none
private

!> Quintic spline interpolator 1d
type, extends(sll_c_interpolator_1d), public :: sll_quintic_spline_interpolator_1d

  sll_real64, dimension(:),   pointer  :: x          !< points position
  sll_int32                            :: n          !< number of points
  sll_int32                            :: ind1       !< number of points
  sll_int32                            :: indn       !< number of points
  sll_real64, dimension(:,:), pointer  :: cf         !< values and derivatives
  sll_real64, dimension(:),   pointer  :: h          !< work array
  sll_int32                            :: bc_min     !< boundary condition
  sll_int32                            :: bc_max     !< boundary condition
  sll_real64                           :: x_min      !< left boundary
  sll_real64                           :: x_max      !< right boundary
  sll_real64                           :: value_min  !< left boundary value
  sll_real64                           :: value_max  !< right boundary value
  sll_real64                           :: slope_min  !< left boundary derivative
  sll_real64                           :: slope_max  !< right boundary derivative
  sll_real64                           :: dx         !< size step

contains

  procedure :: initialize                       !< Initialize
  procedure :: compute_interpolants             !< Compute splines
  procedure :: interpolate_from_interpolant_value              !< Interpolate single value
  procedure :: interpolate_from_interpolant_derivative_eta1      !< Compute derivative
  procedure :: interpolate_from_interpolant_array         !< Interpolate array values
  procedure :: interpolate_from_interpolant_derivatives_eta1    !< Return derivatives
  procedure :: interpolate_array                !< Interpolate an array
  procedure :: interpolate_array_disp           !< Return an array after displacement
  procedure :: set_coefficients                 !< Not implemented
  procedure :: get_coefficients                 !< Not implemented

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
public set_values_at_boundary

contains  ! ****************************************************************

!> Create a new quintic splines interpolator
function new_quintic_spline_interpolator_1d( num_points, &
                                             x_min,      &
                                             x_max,      &
                                             bc_min,     &
                                             bc_max ) result(res)

  sll_interpolator,  pointer :: res

  sll_int32,  intent(in)     :: num_points
  sll_real64, intent(in)     :: x_min
  sll_real64, intent(in)     :: x_max
  sll_int32,  intent(in)     :: bc_min
  sll_int32,  intent(in)     :: bc_max

  call initialize( res,        &
                   num_points, &
                   x_min,      &
                   x_max,      &
                   bc_min,     &
                   bc_max)

end function new_quintic_spline_interpolator_1d

!---------------------------------------------------------------------------

subroutine initialize( interpolator, &
                       num_points,   &
                       x_min,        &
                       x_max,        &
                       bc_min,       &
                       bc_max)

  sll_interpolator,  intent(inout) :: interpolator
  sll_int32,         intent(in   ) :: num_points
  sll_real64,        intent(in   ) :: x_min
  sll_real64,        intent(in   ) :: x_max
  sll_int32,         intent(in   ) :: bc_min
  sll_int32,         intent(in   ) :: bc_max

  character(len=*),  parameter     :: this_sub_name = 'initialize'
  sll_int32                        :: ierr
  sll_int32                        :: i
  sll_real64                       :: dx

  if (bc_min == SLL_PERIODIC .or. bc_max == SLL_PERIODIC) then
    SLL_ERROR( this_sub_name, 'These boundary conditions are not implemented' )
  end if

  interpolator%n      = num_points
  interpolator%x_min  = x_min
  interpolator%x_max  = x_max
  interpolator%bc_min = bc_min
  interpolator%bc_max = bc_max

  SLL_ALLOCATE(interpolator%x(num_points),ierr)
  dx = (x_max - x_min) / (num_points - 1)
  do i = 1, num_points
     interpolator%x(i) = x_min + (i-1)*dx
  end do
  interpolator%dx = dx

  interpolator%ind1 = +1
  interpolator%indn = +1

  SLL_ALLOCATE(interpolator%h(6*num_points-3),ierr)
  SLL_ALLOCATE(interpolator%cf(3,num_points),ierr)

end subroutine

!> Initialization of the boundary for interpolator arbitrary degree splines 1d.
!> @param[in]  value_min  contains the value in the min
!> @param[in]  value_max contains the value in the max
!> @param[in]  slope_min  contains the value in the min for derivative
!> @param[in]  slope_max contains the value in the max for derivative
!> @param[out] interpolator the type sll_arbitrary_degree_spline_interpolator_1d
subroutine set_values_at_boundary( interpolator, &
                                   value_min,    &
                                   value_max,    &
                                   slope_min,    &
                                   slope_max)

  sll_interpolator, intent(inout)  :: interpolator
  sll_real64, intent(in), optional :: value_min
  sll_real64, intent(in), optional :: value_max
  sll_real64, intent(in), optional :: slope_min
  sll_real64, intent(in), optional :: slope_max
  
  if (present(value_min)) interpolator%value_min = value_min
  if (present(value_max)) interpolator%value_max = value_max
  if (present(slope_min)) interpolator%slope_min = slope_min
  if (present(slope_max)) interpolator%slope_max = slope_max
  
end subroutine set_values_at_boundary

!--------------------------------------------------------------------------


subroutine compute_interpolants( interpolator, &
                                 data_array,   &
                                 eta_coords,   &
                                 size_eta_coords)

  sll_interpolator, intent(inout)        :: interpolator
  sll_real64,       intent(in)           :: data_array(:)
  sll_real64,       intent(in), optional :: eta_coords(:)
  sll_int32,        intent(in), optional :: size_eta_coords
  sll_int32                              :: n
  
  n = interpolator%n
  interpolator%cf(1,:) = data_array(:)
  interpolator%cf(1,1) = interpolator%value_min
  interpolator%cf(1,n) = interpolator%value_max
  interpolator%cf(2,1) = interpolator%slope_min
  interpolator%cf(2,n) = interpolator%slope_max 

  if(present(eta_coords) .and. present(size_eta_coords))then
  
    call inspl5(size_eta_coords,  &
                eta_coords,       &
                1,                &
                1,                &
                interpolator%cf,  &
                interpolator%h)

     interpolator%n = size_eta_coords
     interpolator%x = eta_coords
 
  else
  
    call inspl5(interpolator%n,    &
                interpolator%x,    &
                interpolator%ind1, &
                interpolator%indn, &
                interpolator%cf,   &
                interpolator%h)

  end if

end subroutine compute_interpolants

!---------------------------------------------------------------------------

function interpolate_from_interpolant_value( interpolator, eta1 ) result(val)

  sll_interpolator, intent(in) :: interpolator
  sll_real64,       intent(in) :: eta1
  sll_real64                   :: val
  sll_real64                   :: res

  if (eta1 < interpolator%x_min) then
    res = interpolator%x_min
  else if (eta1 > interpolator%x_max) then
    res = interpolator%x_max
  else 
    res = eta1
  end if

  call splin5(interpolator%n,  &
              interpolator%x,  &
              interpolator%cf, &
              res,             &
              0,               &
              val)

end function

!---------------------------------------------------------------------------

subroutine interpolate_array(this,        &
     num_pts,  &
     data,        &
     coordinates, &
     output_array) 

  sll_interpolator,  intent(in)        :: this
  sll_int32,  intent(in)               :: num_pts
  sll_real64, dimension(num_pts), intent(in) :: coordinates
  sll_real64, dimension(:), intent(in) :: data
  sll_real64, dimension(num_pts), intent(out) :: output_array
  
  sll_real64, allocatable :: c(:,:)
  sll_real64, allocatable :: h(:)
  sll_int32               :: ierr
  sll_int32               :: i
  sll_int32               :: n
  
  SLL_CLEAR_ALLOCATE(c(1:3,1:this%n),ierr)
  SLL_CLEAR_ALLOCATE(h(1:6*this%n-3),ierr)
  
  SLL_ASSERT(size(data) == this%n)
  c(1,:) = data
  
  n = num_pts
  c(1,1) = this%value_min
  c(1,n) = this%value_max
  c(2,1) = this%slope_min
  c(2,n) = this%slope_max 
  
  call inspl5(this%n,this%x,this%ind1,this%indn,c,h)
  do i = 1, num_pts
    call splin5(this%n,this%x,c,coordinates(i),0,output_array(i))
  end do

end subroutine interpolate_array


!---------------------------------------------------------------------------

function interpolate_from_interpolant_derivative_eta1( interpolator, eta1 ) result(val)

  sll_interpolator, intent(in) :: interpolator
  sll_real64             :: val
  sll_real64, intent(in) :: eta1

  call splin5(interpolator%n,  &
              interpolator%x,  &
              interpolator%cf, &
              eta1,            &
              1,               &
              val)

end function interpolate_from_interpolant_derivative_eta1


!---------------------------------------------------------------------------

subroutine interpolate_from_interpolant_array( interpolator,        &
     num_pts,             &
     vals_to_interpolate, &
     output_array )

  sll_interpolator,  intent(in)  :: interpolator
  sll_int32,         intent(in)  :: num_pts
  sll_real64,        intent(in)  :: vals_to_interpolate(num_pts)
  sll_real64,        intent(out) :: output_array(num_pts)

  sll_int32                      :: i

  do i = 1, num_pts
    output_array(i) = interpolate_from_interpolant_value(interpolator,vals_to_interpolate(i))
  end do
  
end subroutine interpolate_from_interpolant_array


!---------------------------------------------------------------------------

subroutine interpolate_from_interpolant_derivatives_eta1( interpolator,        &
                                          num_pts,             &
                                          vals_to_interpolate, &
                                          output_array )

  sll_interpolator,         intent(in)  :: interpolator
  sll_int32,                intent(in)  :: num_pts
  sll_real64, dimension(:), intent(in)  :: vals_to_interpolate
  sll_real64, dimension(:), intent(out) :: output_array

  sll_real64, dimension(:,:), allocatable :: c
  sll_real64, dimension(:),   allocatable :: h
  sll_int32                               :: ierr
  sll_int32                               :: i
  sll_int32                               :: n

  SLL_ALLOCATE(c(3,num_pts),ierr)
  SLL_ALLOCATE(h(6*num_pts-3),ierr)
  
  c(1,:) = vals_to_interpolate

  n = num_pts
  c(1,1) = interpolator%value_min
  c(1,n) = interpolator%value_max
  c(2,1) = interpolator%slope_min
  c(2,n) = interpolator%slope_max 

  call inspl5(num_pts, interpolator%x, interpolator%ind1, interpolator%indn, c, h)

  do i = 1, num_pts
    call splin5(num_pts, interpolator%x, c, interpolator%x(i), 1, output_array(i))
  end do

end subroutine interpolate_from_interpolant_derivatives_eta1

!---------------------------------------------------------------------------

subroutine delete_1d( this )

  sll_interpolator :: this

  deallocate(this%x)
  deallocate(this%cf)
  deallocate(this%h)

end subroutine delete_1d

!---------------------------------------------------------------------------

subroutine set_coefficients( interpolator, coeffs )

  sll_interpolator, intent(inout)           :: interpolator
  sll_real64      , intent(in   ), optional :: coeffs(:)

  character(len=*), parameter :: this_sub_name = 'set_coefficients'

  if (present(coeffs)) then
    print*, interpolator%n
  end if

  SLL_ERROR( this_sub_name, 'Not implemented yet.' )

end subroutine set_coefficients

!---------------------------------------------------------------------------

function get_coefficients( interpolator ) result( coeffs )

  sll_interpolator, intent(in) :: interpolator
  sll_real64, pointer          :: coeffs(:)

  character(len=*), parameter  :: this_fun_name = 'get_coefficients'

  coeffs = interpolator%x
  SLL_ERROR( this_fun_name, 'Not implemented yet.' )

end function get_coefficients

!---------------------------------------------------------------------------

subroutine interpolate_array_disp( this,        &
                                 num_pts,  &
                                 data,        &
                                 alpha, &
                                 output_array)

  sll_interpolator, intent(in)    :: this
  sll_int32,        intent(in)    :: num_pts
  sll_real64,       intent(inout) :: data(num_pts)
  sll_real64,       intent(in)    :: alpha
  sll_real64,       intent(inout) :: output_array(num_pts)
  
  character(len=*), parameter  :: this_fun_name = 'interpolate_array_disp'

  print*, this%n
  print*, size(data), alpha, size(output_array), num_pts
  
  SLL_ERROR( this_fun_name, 'Not implemented.' )

end subroutine interpolate_array_disp



end module sll_m_quintic_spline_interpolator_1d
