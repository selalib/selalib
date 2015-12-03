!**************************************************************
!  Copyright INRIA
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
!> @brief
!> Class interpolator and methods for bspline interpolator
!> @details
!> This interpolator works for regular spaced mesh points.
module sll_m_bspline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_errors.h"

use sll_m_bsplines
use sll_m_interpolators_1d_base

implicit none
private

!> Class for arbitrary degree spline 1d interpolator
type, public, extends(sll_c_interpolator_1d) :: sll_bspline_interpolator_1d

  type(sll_bspline_1d), pointer :: bspline    !< bspline data
  sll_int32                     :: num_pts
  sll_int32                     :: spl_deg
  sll_real64                    :: eta_min
  sll_real64                    :: eta_max
  sll_int32                     :: bc_type
  sll_real64                    :: value_l          = 0.0_f64
  logical                       :: compute_value_l  = .false.
  sll_real64                    :: value_r         = 0.0_f64
  logical                       :: compute_value_r = .false.
  sll_real64                    :: slope_l          = 0.0_f64
  logical                       :: compute_slope_l  = .false.
  sll_real64                    :: slope_r         = 0.0_f64
  logical                       :: compute_slope_r = .false.

contains

  !> Initialize the interpolator
  procedure :: initialize=>initialize_bs1d_interpolator
  !> Set spline coefficients
  procedure :: set_coefficients => set_coefficients_bs1d
  !> Compute interpolants
  procedure :: compute_interpolants => compute_interpolants_bs1d
  !> Interpolate single value
  procedure :: interpolate_from_interpolant_value => interpolate_value_bs1d
  !> Interpolate an array (subroutine) 
  procedure :: interpolate_from_interpolant_array => interpolate_values_bs1d
  !> Compute derivatives
  procedure :: interpolate_from_interpolant_derivative_eta1 => interpolate_derivative_bs1d
  !> Compute derivatives array
  procedure :: interpolate_array_derivatives => interpolate_derivatives_bs1d
  !> Interpolate an array (function)
  procedure :: interpolate_array => interpolate_array_bs1d
  !> Interpolate an array after displacement
  procedure :: interpolate_array_disp => interpolate_1d_array_disp_bs1d
  !> Interpolate an array after displacement
  procedure :: interpolate_array_disp_inplace => interpolate_1d_array_disp_inplace_bs1d
  !> Get splines coefficients
  procedure :: get_coefficients => get_coefficients_bs1d
  !> Destory the derived type and free memory
  procedure :: delete => delete_bs1d_interpolator

end type sll_bspline_interpolator_1d

!> Deallocate
interface sll_delete
   module procedure delete_bs1d_interpolator
end interface sll_delete

public sll_delete 
public new_bspline_interpolator_1d
public set_values_at_boundary1d
public initialize_bs1d_interpolator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Delete interpolator arbitrary degree splines.
!> @details
!> The parameters are
!> @param interpolator the type sll_bspline_interpolator_1d
subroutine delete_bs1d_interpolator( interpolator )

class(sll_bspline_interpolator_1d), intent(inout) :: interpolator

call delete_bspline_1d(interpolator%bspline)

end subroutine delete_bs1d_interpolator

!> @brief Initialization of a pointer interpolator arbitrary degree splines 1d.
!> @details To have the interpolator arbitrary degree splines 1d such as a pointer
!>
!> @param[in] num_pts the number of points
!> @param[in] eta_min the minimun
!> @param[in] eta_max the maximun
!> @param[in] bc_l  the boundary condition at left
!> @param[in] bc_r the boundary condition at right
!> @param[in] spl_deg the degree of B-spline
!> @return the type interpolator arbitrary degree splines 1d

function new_bspline_interpolator_1d( &
  num_pts,                            &
  eta_min,                            &
  eta_max,                            &
  spl_deg,                            &
  bc_type,                            &
  slope_l,                            &
  slope_r) result(interpolator)

class(sll_bspline_interpolator_1d),pointer :: interpolator

sll_int32,  intent(in) :: num_pts
sll_real64, intent(in) :: eta_min
sll_real64, intent(in) :: eta_max
sll_int32,  intent(in) :: bc_type
sll_int32,  intent(in) :: spl_deg
sll_real64, optional   :: slope_l
sll_real64, optional   :: slope_r

sll_int32              :: ierr

SLL_ALLOCATE(interpolator,ierr)


if ( present(slope_l) .and. present(slope_r) ) then

  call initialize_bs1d_interpolator( interpolator,  &
                                     num_pts,  &
                                     eta_min,  &
                                     eta_max,  &
                                     spl_deg,  &
                                     bc_type,  &
                                     slope_l,  &
                                     slope_r)
else

  call initialize_bs1d_interpolator( interpolator,  &
                                     num_pts,  &
                                     eta_min,  &
                                     eta_max,  &
                                     bc_type,  &
                                     spl_deg  )

end if
                                           

end function new_bspline_interpolator_1d

!> @brief Initialization of interpolator arbitrary degree splines 1d.
!> @details To have the interpolator arbitrary degree splines 1d
!>
!> @param[in]  num_pts  the number of points
!> @param[in]  eta_min  the minimun
!> @param[in]  eta_max  the maximun
!> @param[in]  bc_type  the boundary condition (periodic or not)
!> @param[in]  spl_deg  the degree of B-spline
!> @param[out] interpolator the type sll_bspline_interpolator_1d

subroutine initialize_bs1d_interpolator( interpolator,  &
                                         num_pts,  &
                                         eta_min,  &
                                         eta_max,  &
                                         spl_deg,  &
                                         bc_type,  &
                                         slope_l,  &
                                         slope_r)

class(sll_bspline_interpolator_1d), intent(inout) :: interpolator

sll_int32,       intent(in) :: num_pts
sll_real64,      intent(in) :: eta_min
sll_real64,      intent(in) :: eta_max
sll_int32,       intent(in) :: spl_deg
sll_int32,       intent(in) :: bc_type
sll_real64,      optional   :: slope_l
sll_real64,      optional   :: slope_r

interpolator%num_pts = num_pts
interpolator%spl_deg = spl_deg
interpolator%eta_min = eta_min
interpolator%eta_max = eta_max

if (present(slope_l) .and. present(slope_r)) then

  interpolator%bspline => new_bspline_1d( num_pts, &
                                          spl_deg, &
                                          eta_min, &
                                          eta_max, &
                                          bc_type, &
                                          slope_l, &
                                          slope_r)
else

  interpolator%bspline => new_bspline_1d( num_pts, &
                                          spl_deg, &
                                          eta_min, &
                                          eta_max, &
                                          bc_type)
end if

end subroutine initialize_bs1d_interpolator


!> Set values at the boundaries for the interpolator 1d.
!> The parameters are
!> @param[in]  value_l contains the value in the left
!> @param[in]  value_r contains the value in the right
!> @param[in]  slope_l contains the value in the left for derivative
!> @param[in]  slope_r contains the value in the right for derivative
!> @param[out] interpolator the type sll_bspline_interpolator_1d
subroutine set_values_at_boundary1d( interpolator, &
                                     value_l, &
                                     value_r, &
                                     slope_l, &
                                     slope_r)

class(sll_bspline_interpolator_1d), intent(inout) :: interpolator

sll_real64, intent(in), optional :: value_l
sll_real64, intent(in), optional :: value_r
sll_real64, intent(in), optional :: slope_l
sll_real64, intent(in), optional :: slope_r

if (present(value_l)) then
  interpolator%value_l = value_l
  interpolator%compute_value_l = .false.
end if

if (present(value_r)) then
  interpolator%value_r = value_r
  interpolator%compute_value_r = .false.
end if

if (present(slope_l)) then
  interpolator%slope_l = slope_l
  interpolator%compute_slope_l = .false.
end if

if (present(slope_r)) then
  interpolator%slope_r = slope_r
  interpolator%compute_slope_r = .false.
end if

end subroutine set_values_at_boundary1d

!> @brief computing the coefficients spline with a given
!>  data_array 1D cooresponding at the values of a function
!> @details 
!>  on eta_coords of size size_eta_coords
!>  if the eta_coords and eta_coords is not given
!>  we consider that the values of the function is on the points in the mesh_1d
!>
!> The parameters are
!> @param[in]  data_array the 1d arrays corresponding at the values of a function
!> @param[in]  eta_coords the 1d arrays
!> @param[in]  size_eta_coords the size of eta_coords
!> @param[out] interpolator the type sll_bspline_interpolator_1d
subroutine compute_interpolants_bs1d( interpolator,    &
                                      data_array,      &
                                      eta_coords,      &
                                      size_eta_coords)

class(sll_bspline_interpolator_1d), &
            intent(inout)           :: interpolator
sll_real64, intent(in   )           :: data_array(:)
sll_real64, intent(in   ), optional :: eta_coords(:)
sll_int32,  intent(in   ), optional :: size_eta_coords

character(len=*), parameter :: this_sub_name = 'compute_interpolants_bs1d'

if(present(eta_coords) .or. present(size_eta_coords)) then
   SLL_ERROR( this_sub_name, 'This case is not yet implemented' )
end if

call compute_bspline_1d ( interpolator%bspline, data_array )

end subroutine compute_interpolants_bs1d

    
!> @brief Interpolation on the points eta using
!> the arbitrary degree splines interpolator 1d
!> @details computing the values with the interpolator 
!> arbitrary degree splines 1d
!> on the points eta of arbitrary degree splines 1d
!> @param[in] interpolator the type sll_bspline_interpolator_1d
!> @param[in] eta1 the point
!> @return val the values on the points eta
function interpolate_value_bs1d( interpolator, eta1) result(val)

class(sll_bspline_interpolator_1d), intent(in)  :: interpolator

sll_real64, intent(in)          :: eta1
sll_real64                      :: val
sll_real64                      :: res

res = eta1

if (interpolator%bc_type == SLL_PERIODIC) then ! periodic

  if( res < interpolator%eta_min ) then
     res = res+interpolator%bspline%length
  else if( res >  interpolator%eta_max ) then
     res = res-interpolator%bspline%length
  end if

end if

val = interpolate_value_1d( interpolator%bspline, res)

end function interpolate_value_bs1d


!> @brief initializing the coefficients of splines.
!> @details  initializing the coefficients of splines
!>  fot the arbitrary degree splines interpolator 1d
!> The parameters are
!> @param interpolator the type sll_bspline_interpolator_1d
!> @param[in] coeffs the 1d arrays corresponding of the splines coefficients
!> @param[out] interpolator the type sll_bspline_interpolator_1d

subroutine set_coefficients_bs1d( interpolator, coeffs)

class(sll_bspline_interpolator_1d), intent(inout)  :: interpolator
sll_real64, dimension(:), optional, intent(in)     :: coeffs

print*, 'num_pts =', interpolator%num_pts
if (present(coeffs)) print*, size(coeffs)
stop ' set_coefficients_bs1d not implemented '

end subroutine set_coefficients_bs1d

  
!> @brief First derivative interpolation on the point eta
!> @details computing the values of the first derivative
!> with the interpolator arbitrary degree splines 1d
!> on the points eta of arbitrary degree splines 1d
!>
!> The parameters are
!> @param interpolator the type sll_bspline_interpolator_1d
!> @param[in] eta1 the point
!> @return val the values on the point eta of the first derivative

function interpolate_derivative_bs1d( interpolator, eta1 ) result(val)

class(sll_bspline_interpolator_1d), intent(in)  :: interpolator

sll_real64, intent(in)           :: eta1
sll_real64                       :: val
sll_real64                       :: res

res = eta1

if (interpolator%bspline%bc_type == SLL_PERIODIC ) then 

  if( res < interpolator%eta_min ) then
    res = res+interpolator%bspline%length
  else if( res >  interpolator%eta_max ) then
    res = res-interpolator%bspline%length
  end if

end if

SLL_ASSERT( res >= interpolator%eta_min )
SLL_ASSERT( res <= interpolator%eta_max )

val = interpolate_derivative_1d( interpolator%bspline, res )

end function interpolate_derivative_bs1d

subroutine interpolate_array_bs1d( this,         &
                                 num_pts,   &
                                 data,   &
                                 coordinates, &
                                 output_array)

class(sll_bspline_interpolator_1d), intent(in) :: this
sll_int32,  intent(in)               :: num_pts
sll_real64, dimension(num_pts), intent(in) :: coordinates
sll_real64, dimension(:), intent(in) :: data
sll_real64, dimension(num_pts), intent(out)    :: output_array

call compute_bspline_1d( this%bspline, data)
call interpolate_array_values_1d( this%bspline, num_pts, coordinates, output_array)

end subroutine interpolate_array_bs1d

subroutine interpolate_1d_array_disp_bs1d( this,       &
                                         num_pts, &
                                         data,       &
                                         alpha, &
                                         output_array)

class(sll_bspline_interpolator_1d), intent(in) :: this
sll_int32,                          intent(in) :: num_pts
sll_real64, dimension(:),           intent(in) :: data
sll_real64, intent(in)                         :: alpha
sll_real64, dimension(num_pts),intent(out)              :: output_array

output_array = -1000000._f64*alpha*data*this%spl_deg
stop 'interpolate_1d_array_disp_bs1d: not implemented.'

end subroutine interpolate_1d_array_disp_bs1d


subroutine interpolate_1d_array_disp_inplace_bs1d( this,       &
                                         num_pts, &
                                         data,       &
                                         alpha)

class(sll_bspline_interpolator_1d), intent(in) :: this
sll_int32,                          intent(in) :: num_pts
sll_real64, dimension(num_pts),           intent(inout) :: data
sll_real64, intent(in)                         :: alpha

stop 'interpolate_1d_array_disp_bs1d: not implemented.'

end subroutine interpolate_1d_array_disp_inplace_bs1d


function get_coefficients_bs1d(interpolator)

class(sll_bspline_interpolator_1d), intent(in)    :: interpolator
sll_real64, dimension(:), pointer            :: get_coefficients_bs1d

get_coefficients_bs1d => interpolator%bspline%bcoef

end function get_coefficients_bs1d

subroutine interpolate_values_bs1d( interpolator,        &
                                    num_pts,             &
                                    vals_to_interpolate, &
                                    output_array )

class(sll_bspline_interpolator_1d), intent(in)  :: interpolator
sll_int32,                          intent(in)  :: num_pts
sll_real64, dimension(num_pts),           intent(in)  :: vals_to_interpolate
sll_real64, dimension(num_pts),           intent(out) :: output_array

call interpolate_array_values_1d(interpolator%bspline, &
                              num_pts,              &
                              vals_to_interpolate,  &
                              output_array)

end subroutine interpolate_values_bs1d


subroutine interpolate_derivatives_bs1d( interpolator,        &
                                         num_pts,             &
                                         vals_to_interpolate, &
                                         output_array )

class(sll_bspline_interpolator_1d), intent(in)  :: interpolator
sll_int32,                          intent(in)  :: num_pts
sll_real64, dimension(:),           intent(in)  :: vals_to_interpolate
sll_real64, dimension(:),           intent(out) :: output_array

call interpolate_array_derivatives_1d(interpolator%bspline, &
                                   num_pts,              &
                                   vals_to_interpolate,  &
                                   output_array)

end subroutine interpolate_derivatives_bs1d


end module sll_m_bspline_interpolator_1d
