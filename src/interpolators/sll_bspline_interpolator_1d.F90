!**************************************************************
!  Copyright INRIA
!  Authors :
!     Aurore
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
module sll_module_bspline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_utilities.h"

use sll_bsplines
use sll_module_interpolators_1d_base

implicit none
private

!> Class for arbitrary degree spline 1d interpolator
type, public, extends(sll_interpolator_1d_base) :: sll_bspline_interpolator_1d

  type(sll_bspline_1d) :: bspline             !< bspline data
  sll_int32            :: num_pts             !< nodes number
  sll_real64           :: eta_min             !< left boundary
  sll_real64           :: eta_max             !< right boundary
  sll_int32            :: bc_left             !< left boundary type
  sll_int32            :: bc_right            !< right boundary type
  sll_int32            :: bc_selector         !< boundary combination
  sll_int32            :: spline_degree       !< spline degree
  sll_real64, pointer  :: eta                 !< node positions
  sll_int32            :: size_coeffs         !< coeffs array dimension
  sll_real64           :: slope_left          !< left boundary derivative
  sll_real64           :: slope_right         !< right boundary derivative
  sll_real64           :: value_left          !< left boundary value
  sll_real64           :: value_right         !< right boundary value
  logical              :: compute_slope_left  !< true
  logical              :: compute_slope_right !< true
  logical              :: compute_value_left  !< true
  logical              :: compute_value_right !< true

contains

  !> Initialize the interpolator
  procedure :: initialize=>initialize_bs1d_interpolator
  !> Set spline coefficients
  procedure :: set_coefficients => set_coefficients_bs1d
  !> Compute interpolants
  procedure :: compute_interpolants => compute_interpolants_bs1d
  !> Interpolate single value
  procedure :: interpolate_value => interpolate_value_bs1d
  !> Interpolate an array (subroutine) 
  procedure :: interpolate_array_values => interpolate_values_bs1d
  !> Interpolate a pointer to array 
  procedure :: interpolate_pointer_values => interpolate_pointer_values_bs1d
  !> Compute derivatives
  procedure :: interpolate_derivative_eta1 => interpolate_derivative_bs1d
  !> Compute derivatives array
  procedure :: interpolate_array_derivatives => interpolate_derivatives_bs1d
  !> Compute derivatives array pointer
  procedure :: interpolate_pointer_derivatives =>interpolate_pointer_derivatives_bs1d
  !> Interpolate an array (function)
  procedure :: interpolate_array => interpolate_array_bs1d
  !> Interpolate an array after displacement
  procedure :: interpolate_array_disp => interpolate_1d_array_disp_bs1d
  !> Get splines coefficients
  procedure :: get_coefficients => get_coefficients_bs1d
  !> Not implemented
  procedure :: reconstruct_array
  !> Destory the derived type and free memory
  procedure :: delete => delete_b1d_interpolator

end type sll_bspline_interpolator_1d

!> Deallocate
interface sll_delete
   module procedure delete_b1d_interpolator
end interface sll_delete

public sll_delete 
public new_b1d_interpolator
public set_values_at_boundary1d
public initialize_bs1d_interpolator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> @brief Delete interpolator arbitrary degree splines.
!> @details
!> The parameters are
!> @param interpolator the type sll_bspline_interpolator_1d
subroutine delete_b1d_interpolator( interpolator )

class(sll_bspline_interpolator_1d), intent(inout) :: interpolator
sll_int32 :: ierr

end subroutine delete_b1d_interpolator

!> @brief Initialization of a pointer interpolator arbitrary degree splines 1d.
!> @details To have the interpolator arbitrary degree splines 1d such as a pointer
!>
!> @param[in] num_pts the number of points
!> @param[in] eta_min the minimun
!> @param[in] eta_max the maximun
!> @param[in] bc_left  the boundary condition at left
!> @param[in] bc_right the boundary condition at right
!> @param[in] spline_degree the degree of B-spline
!> @return the type interpolator arbitrary degree splines 1d

function new_b1d_interpolator( num_pts,       &
                                               eta_min,       &
                                               eta_max,       &
                                               bc_left,       &
                                               bc_right,      &
                                               spline_degree) &
result(interpolator)

class(sll_bspline_interpolator_1d),pointer :: interpolator
sll_int32,  intent(in) :: num_pts
sll_real64, intent(in) :: eta_min
sll_real64, intent(in) :: eta_max
sll_int32,  intent(in) :: bc_left
sll_int32,  intent(in) :: bc_right
sll_int32,  intent(in) :: spline_degree

sll_int32              :: ierr

SLL_ALLOCATE(interpolator,ierr)

call initialize_bs1d_interpolator( interpolator, &
                                   num_pts,      &
                                   eta_min,      &
                                   eta_max,      &
                                   bc_left,      &
                                   bc_right,     &
                                   spline_degree)

end function new_b1d_interpolator

!> @brief Initialization of interpolator arbitrary degree splines 1d.
!> @details To have the interpolator arbitrary degree splines 1d
!>
!> @param[in]  num_pts the number of points
!> @param[in]  eta_min the minimun
!> @param[in]  eta_max the maximun
!> @param[in]  bc_left  the boundary condition at left
!> @param[in]  bc_right the boundary condition at right
!> @param[in]  spline_degree the degree of B-spline
!> @param[out] interpolator the type sll_bspline_interpolator_1d

subroutine initialize_bs1d_interpolator( interpolator, &
                                         num_pts,      &
                                         eta_min,      &
                                         eta_max,      &
                                         bc_left,      &
                                         bc_right,     &
                                         spline_degree)

class(sll_bspline_interpolator_1d), intent(inout) :: interpolator

sll_int32,       intent(in) :: num_pts
sll_real64,      intent(in) :: eta_min
sll_real64,      intent(in) :: eta_max
sll_int32,       intent(in) :: bc_left
sll_int32,       intent(in) :: bc_right
sll_int32,       intent(in) :: spline_degree

sll_int32                   :: ierr
sll_int32                   :: tmp
sll_int32                   :: bc_selector
sll_int32                   :: i, k
sll_real64                  :: delta_eta
character(len=*), parameter :: this_sub_name = 'initialize_bs1d_interpolator'

! do some argument checking...
if(((bc_left == SLL_PERIODIC).and.(bc_right.ne. SLL_PERIODIC))) then
   print *, 'initialize_b1d_interpolator, ERROR: ', &
        'if one boundary condition is specified as periodic, then ', &
        'both must be. Error in first direction.'
end if

bc_selector = 0

if( bc_left  == SLL_DIRICHLET ) bc_selector = bc_selector + 1
if( bc_left  == SLL_NEUMANN   ) bc_selector = bc_selector + 2
if( bc_left  == SLL_HERMITE   ) bc_selector = bc_selector + 4
if( bc_right == SLL_DIRICHLET ) bc_selector = bc_selector + 8
if( bc_right == SLL_NEUMANN   ) bc_selector = bc_selector + 16
if( bc_right == SLL_HERMITE   ) bc_selector = bc_selector + 32

interpolator%spline_degree       = spline_degree
interpolator%eta_min             = eta_min
interpolator%eta_max             = eta_max
interpolator%bc_left             = bc_left
interpolator%bc_right            = bc_right
interpolator%bc_selector         = bc_selector
interpolator%num_pts             = num_pts
interpolator%compute_slope_left  = .true.
interpolator%compute_slope_right = .true.
interpolator%compute_value_left  = .true.
interpolator%compute_value_right = .true.

k = spline_degree+1

select case (bc_selector)

case (0) ! 1. periodic

  call initialize_bspline_1d(interpolator%bspline, &
                             num_pts,              &
                             k,                    &
                             eta_min,              &
                             eta_max,              &
                             SLL_PERIODIC          )

case (9) ! 2. dirichlet-left, dirichlet-right

  interpolator%value_left  = 0.0_f64
  interpolator%value_right = 0.0_f64

  call initialize_bspline_1d(interpolator%bspline, &
                             num_pts,              &
                             k,                    &
                             eta_min,              &
                             eta_max,              &
                             SLL_DIRICHLET         )

case default

  call initialize_bspline_1d(interpolator%bspline, &
                             num_pts,              &
                             k,                    &
                             eta_min,              &
                             eta_max,              &
                             SLL_HERMITE           )

  interpolator%slope_right = 0.0_f64
  interpolator%slope_left = 0.0_f64
    
end select

end subroutine initialize_bs1d_interpolator


!> Initialization of the boundary for interpolator arbitrary degree splines 1d.
!> The parameters are
!> @param[in]  value_left  contains the value in the left
!> @param[in]  value_right contains the value in the right
!> @param[in]  slope_left  contains the value in the left for derivative
!> @param[in]  slope_right contains the value in the right for derivative
!> @param[out] interpolator the type sll_bspline_interpolator_1d
subroutine set_values_at_boundary1d( interpolator, &
                                     value_left,   &
                                     value_right,  &
                                     slope_left,   &
                                     slope_right)

class(sll_bspline_interpolator_1d), intent(inout) :: interpolator

sll_real64, intent(in), optional :: value_left
sll_real64, intent(in), optional :: value_right
sll_real64, intent(in), optional :: slope_left
sll_real64, intent(in), optional :: slope_right

sll_int32 :: bc_left
sll_int32 :: bc_right
sll_int64 :: bc_selector

bc_left = interpolator%bc_left
bc_right= interpolator%bc_right
bc_selector = interpolator%bc_selector

if (present(value_left)) then
  interpolator%value_left = value_left
  interpolator%compute_value_left = .FALSE.
end if

if (present(value_right)) then
  interpolator%value_right = value_right
  interpolator%compute_value_right = .FALSE.
end if

if (present(slope_left)) then
  interpolator%slope_left = slope_left
  interpolator%compute_slope_left = .FALSE.
end if

if (present(slope_right)) then
  interpolator%slope_right = slope_right
  interpolator%compute_slope_right = .FALSE.
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
sll_int32       , parameter :: sz_deriv = 2
sll_int32                   :: point_locate_eta_derivative( sz_deriv )
sll_real64                  :: data_array_derivative      ( sz_deriv )
sll_int32                   :: sz
sll_int32                   :: k
sll_real64                  :: period
sll_int32                   :: order
sll_int32                   :: ierr

k = interpolator%bspline%k

if(present(eta_coords) .or. present(size_eta_coords)) then
   SLL_ERROR( this_sub_name, 'This case is not yet implemented' )
end if

sz = interpolator%num_pts

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
sll_int32                       :: size_coeffs
sll_real64                      :: res

size_coeffs = interpolator%size_coeffs

res = eta1

if (interpolator%bc_selector == 0) then ! periodic

  if( res < interpolator%eta_min ) then
     res = res+interpolator%eta_max-interpolator%eta_min
  else if( res >  interpolator%eta_max ) then
     res = res+interpolator%eta_min-interpolator%eta_max
  end if

end if

val = interpolate_value( interpolator%bspline, res)

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
sll_real64, dimension(:), intent(in), optional :: coeffs

sll_int32  :: sp_deg
sll_int32  :: num_cells
sll_int32  :: tmp
sll_int32  :: i
sll_real64 :: eta_min
sll_real64 :: eta_max
sll_real64 :: delta
sll_int32  :: nb_spline_eta
sll_real64 :: eta

sp_deg    = interpolator%spline_degree
num_cells = interpolator%num_pts - 1
eta_min   = interpolator%eta_min
eta_max   = interpolator%eta_max
delta     = (eta_max - eta_min)/num_cells

tmp = (sp_deg + 1)/2

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

if (interpolator%bc_selector == 0 ) then !periodic

  if( res < interpolator%eta_min ) then
    res = res+interpolator%eta_max-interpolator%eta_min
  else if( res >  interpolator%eta_max ) then
    res = res+interpolator%eta_min-interpolator%eta_max
  end if

end if

SLL_ASSERT( res >= interpolator%eta_min )
SLL_ASSERT( res <= interpolator%eta_max )

val = interpolate_derivative( interpolator%bspline, res )

end function interpolate_derivative_bs1d

function interpolate_array_bs1d( this,         &
                                 num_points,   &
                                 data,   &
                                 coordinates) result(res)

class(sll_bspline_interpolator_1d), intent(in) :: this

sll_int32,  intent(in)               :: num_points
sll_real64, dimension(:), intent(in) :: coordinates
sll_real64, dimension(:), intent(in) :: data
sll_real64, dimension(num_points)    :: res
sll_int32                            :: i

call compute_bspline_1d( this%bspline, data)
call interpolate_array_values( this%bspline, num_points, coordinates, res)

end function interpolate_array_bs1d

function interpolate_1d_array_disp_bs1d( &
     this,        &
     num_points, &
     data,     &
     alpha) result(res)

class(sll_bspline_interpolator_1d), intent(in)    :: this
sll_int32, intent(in)                          :: num_points
sll_real64, dimension(:), intent(in)         :: data
sll_real64, intent(in)         :: alpha
sll_real64, dimension(num_points) :: res

print *, 'interpolate_1d_array_disp_bs1d: not implemented.'
res = -1000000._f64*alpha*data*this%spline_degree

end function interpolate_1d_array_disp_bs1d

function get_coefficients_bs1d(interpolator)

class(sll_bspline_interpolator_1d), intent(in)    :: interpolator
sll_real64, dimension(:), pointer            :: get_coefficients_bs1d

get_coefficients_bs1d => interpolator%bspline%bcoef

end function get_coefficients_bs1d

subroutine interpolate_values_bs1d( interpolator,        &
                                    num_pts,             &
                                    vals_to_interpolate, &
                                    output_array )

class(sll_bspline_interpolator_1d),  intent(in) :: interpolator
sll_int32,  intent(in)                 :: num_pts
sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
sll_real64, dimension(:), intent(out)  :: output_array

call interpolate_array_values(interpolator%bspline, &
                              num_pts,              &
                              vals_to_interpolate,  &
                              output_array)

end subroutine interpolate_values_bs1d

subroutine interpolate_pointer_values_bs1d( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output )

class(sll_bspline_interpolator_1d),  intent(in) :: interpolator
sll_int32,  intent(in)            :: num_pts
sll_real64, dimension(:), pointer :: vals_to_interpolate
sll_real64, dimension(:), pointer :: output
sll_int32 :: idx

SLL_ASSERT(num_pts==size(vals_to_interpolate))
do idx=1,num_pts
      output(idx)=interpolate_value_bs1d( &
                            interpolator, &
                            vals_to_interpolate(idx))
enddo

end subroutine interpolate_pointer_values_bs1d

subroutine interpolate_derivatives_bs1d( &
     interpolator, &
     num_pts, &
     vals_to_interpolate, &
     output_array )

  class(sll_bspline_interpolator_1d),  intent(in) :: interpolator
  sll_int32,  intent(in)                 :: num_pts
  sll_real64, dimension(:), intent(in) :: vals_to_interpolate
  sll_real64, dimension(:), intent(out) :: output_array
  sll_int32 :: idx

  SLL_ASSERT(num_pts==size(vals_to_interpolate))
  do idx=1,num_pts
        output_array(idx)=interpolate_derivative_bs1d( &
                              interpolator, &
                              vals_to_interpolate(idx))
  enddo

end subroutine interpolate_derivatives_bs1d

subroutine interpolate_pointer_derivatives_bs1d( &
     interpolator, &
     num_pts, &
     vals_to_interpolate, &
     output )

class(sll_bspline_interpolator_1d),  intent(in) :: interpolator
sll_int32,  intent(in)              :: num_pts
sll_real64, dimension(:), pointer   :: vals_to_interpolate
sll_real64, dimension(:), pointer   :: output
sll_int32 :: idx

SLL_ASSERT(num_pts==size(vals_to_interpolate))
do idx=1,num_pts
      output(idx)=interpolate_derivative_bs1d( &
                            interpolator, &
                            vals_to_interpolate(idx))
enddo

end subroutine interpolate_pointer_derivatives_bs1d

function reconstruct_array(this, num_points, data) result(res)

  class(sll_bspline_interpolator_1d),  intent(in) :: this
  sll_int32, intent(in)                :: num_points! size of output array
  sll_real64, dimension(:), intent(in) :: data   ! data to be interpolated
  sll_real64, dimension(num_points)    :: res
  res(:) = -1000000.0_f64*data*this%spline_degree
  print*, 'reconstruct_array 1d not implemented yet'

end function reconstruct_array

end module sll_module_bspline_interpolator_1d
