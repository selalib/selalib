module sll_cubic_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#ifndef STDF95
use sll_module_interpolators_1d_base
#endif
use sll_splines
  implicit none

#ifdef STDF95
  type                                    ::  cubic_spline_1d_interpolator
#else  
  type, extends(sll_interpolator_1d_base) ::  cubic_spline_1d_interpolator
#endif
     sll_real64, dimension(:), allocatable :: interpolation_points 
     sll_int32                     :: num_points ! size
     sll_int32                     :: bc_type
     type(sll_spline_1D), pointer  :: spline
#ifdef STDF95
#else
   contains
     procedure, pass(interpolator) :: initialize => initialize_cs1d_interpolator
     procedure :: compute_interpolants => compute_interpolants_cs1d
     procedure :: interpolate_value => interpolate_value_cs1d
     procedure :: interpolate_derivative_eta1 => interpolate_deriv1_cs1d
     procedure, pass:: interpolate_array => spline_interpolate1d
     procedure, pass:: interpolate_array_disp => spline_interpolate1d_disp
     procedure, pass:: reconstruct_array
     !generic :: initialize => initialize_cs1d_interpolator
#endif
  end type cubic_spline_1d_interpolator


  
contains  ! ****************************************************************


  ! the following provides an implementation for the abstract interface 
  !interpolate1d
  !> Define spline interpolation of values in data define on original grid at 
  !> points coordinates
  ! Issues with the following function:
  ! - entities referenced through "this" are modified, violating the declared
  !   intent.
  ! - it is probably better to convert this into a subroutine, since data_out
  !   will be allocated on the stack (too big an array will crash the program),
  !   and some copy operation might be involved when "catching" the results.

#ifdef STDF95
  function cubic_spline_interpolate_array(this, num_points, data, coordinates) &
       result(data_out)
    type(cubic_spline_1d_interpolator),  intent(in)       :: this
#else
  function spline_interpolate1d(this, num_points, data, coordinates) &
       result(data_out)
    class(cubic_spline_1d_interpolator),  intent(in)       :: this
#endif
    !class(sll_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    ! local variables
    sll_int32 :: ierr
    ! compute the interpolating spline coefficients
    call compute_spline_1D( data, this%bc_type, this%spline )
    call interpolate_array_values( coordinates, data_out, num_points, &
         this%spline )
  end function 

#ifdef STDF95
  function cubic_spline_interpolate_array_at_displacement(this, num_points, &
       data, coordinates) &
       result(data_out)
    type(cubic_spline_1d_interpolator),  intent(in)       :: this
#else
  function spline_interpolate1d_disp(this, num_points, data, alpha) &
       result(data_out)
    class(cubic_spline_1d_interpolator),  intent(in)       :: this
#endif
    !class(sll_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64,  intent(in)   :: alpha
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    ! local variables
    sll_real64, dimension(num_points)      :: coordinates
    sll_real64 :: length, delta
    sll_int32 :: i
    sll_int32 :: ierr
    ! compute the interpolating spline coefficients
    call compute_spline_1D( data, this%bc_type, this%spline )
    ! compute array of coordinates where interpolation is performed from displacement
    length = this%interpolation_points(num_points) - this%interpolation_points(1)
    delta = this%interpolation_points(2) - this%interpolation_points(1)
!    if (this%bc_type == PERIODIC_SPLINE) then
       coordinates(1) = this%interpolation_points(1) + modulo(this%interpolation_points(1) - alpha, length)
       print*, coordinates(1), this%interpolation_points(1), this%interpolation_points(num_points)
       do i = 2, num_points      
          coordinates(i) = modulo(coordinates(i-1) + delta, length)
          print*, coordinates(i)
       end do
!!$    else
!!$       if (alpha < 0 ) then 
!!$          coordinates(1) = this%interpolation_points(1)
!!$          do i = 2, num_points
!!$             coordinates(i) = max(coordinates(i-1) + delta, this%interpolation_points(1))
!!$          end do
!!$       else
!!$          coordinates(num_points) = this%interpolation_points(num_points)
!!$          do i = num_points-1, 1, -1
!!$             coordinates(i) = min(coordinates(i+1) - delta, &
!!$                  this%interpolation_points(num_points))
!!$          end do
!!$       endif
!!$    end if
    call interpolate_array_values( coordinates, data_out, num_points, &
         this%spline )
  end function
  
#ifdef STDF95
  subroutine cubic_spline_compute_interpolants( interpolator, data_array )
    type(cubic_spline_1d_interpolator), intent(inout)  :: interpolator
#else
  subroutine compute_interpolants_cs1d( interpolator, data_array )
    class(cubic_spline_1d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64, dimension(:), intent(in)               :: data_array
    call compute_spline_1D_bis( data_array, interpolator%spline )
  end subroutine

#ifdef STDF95
  function cubic_spline_interpolate_value( interpolator, eta1 ) result(val)
    type(cubic_spline_1d_interpolator), intent(inout) :: interpolator
#else
  function interpolate_value_cs1d( interpolator, eta1 ) result(val)
    class(cubic_spline_1d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    val = interpolate_value( eta1, interpolator%spline )
  end function
  
#ifdef STDF95
  function cubic_spline_interpolate_derivative_eta1( interpolator, eta1 ) result(val)
    type(cubic_spline_1d_interpolator), intent(inout)  :: interpolator
#else
  function interpolate_deriv1_cs1d( interpolator, eta1 ) result(val)
    class(cubic_spline_1d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    val = interpolate_derivative(eta1,interpolator%spline)
  end function

  !> create new spline object
#ifdef STDF95
  subroutine cubic_spline_initialize( &
#else
  subroutine initialize_cs1d_interpolator( &
#endif
    interpolator, &
    num_points, &
    xmin, &
    xmax, &
    bc_type, &
    slope_left, &
    slope_right )

#ifdef STDF95
    type(cubic_spline_1d_interpolator),  intent(inout)  :: interpolator 
#else
    class(cubic_spline_1d_interpolator),  intent(inout) :: interpolator 
#endif
    sll_int32,  intent(in)               :: num_points
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type
    sll_real64, intent(in), optional     :: slope_left
    sll_real64, intent(in), optional     :: slope_right
    sll_int32                            :: ierr
    sll_int32  :: i  
    sll_real64 :: delta
    
    interpolator%num_points = num_points
    SLL_ALLOCATE(interpolator%interpolation_points(num_points),ierr)
    interpolator%interpolation_points(1) = xmin
    delta = (xmax - xmin) / (num_points - 1)
    do i = 2, num_points
       interpolator%interpolation_points(i) = interpolator%interpolation_points(i-1) + delta
    end do
    interpolator%bc_type = bc_type
    if (present(slope_left).and.present(slope_right)) then
       interpolator%spline => new_spline_1D(num_points, xmin, xmax, bc_type, slope_left, slope_right )
    else
       interpolator%spline => new_spline_1D(num_points, xmin, xmax, bc_type)
    end if
  end subroutine

  function reconstruct_array(this, num_points, data) result(res)
    ! dummy procedure
#ifdef STDF95
    type(cubic_spline_1d_interpolator), intent(in)      :: this
#else
    class(cubic_spline_1d_interpolator), intent(in)     :: this
#endif
       sll_int32, intent(in)                :: num_points! size of output array
       sll_real64, dimension(:), intent(in) :: data   ! data to be interpolated 
       sll_real64, dimension(num_points)    :: res
       res(:) = 0.0_f64
  end function reconstruct_array
  
end module sll_cubic_spline_interpolator_1d
