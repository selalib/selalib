!> \brief data type for 1D interpolation
!>
!> 


module sll_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_splines
use sll_weno
  implicit none

  type :: interpolator_1d
     type(sll_spline_1D), pointer :: spl
     type(sll_weno_1d), pointer   :: weno
     procedure(realarray_function_interpolator1d_int_realarray_realarray), pointer, nopass:: interpolate_1d
  end type interpolator_1d

  abstract interface
     function realarray_function_interpolator1d_int_realarray_realarray(interp, num_points, coordinates, data) result(res)
       use sll_working_precision
       import interpolator_1d
       type(interpolator_1d), pointer, intent(in)     :: interp
       sll_int32, intent(in)                         :: num_points    ! number of interpolation points
       sll_real64, dimension(:), intent(in)          :: coordinates   ! points where output is desired
       sll_real64, dimension(:), intent(in)          :: data          ! data to be interpolated 
       sll_real64, dimension(num_points)             :: res
     end function realarray_function_interpolator1d_int_realarray_realarray
  end interface
contains
  function new_interpolator_1d(what_interpolator, num_points, xmin, xmax, bc_type) result(interp)
    type(interpolator_1d), pointer :: interp
    character(len=*) :: what_interpolator
    sll_int32 :: num_points
    sll_real64 :: xmin, xmax
    sll_int32 :: bc_type
    sll_int32 :: ierr

    SLL_ALLOCATE(interp,ierr)
    select case(what_interpolator)
    case ('spline')
       SLL_ALLOCATE(interp%spl,ierr)
       interp%spl => new_spline_1d( num_points, xmin, xmax, bc_type )
       interp%interpolate_1d => interpolate_1d_spline
    case ('weno')
       SLL_ALLOCATE(interp%weno,ierr)
       interp%weno => new_WENO_1D(num_points, xmin, xmax)
       interp%interpolate_1d => interpolate_1d_weno
    case default
       print*, 'interpolator ', what_interpolator, ' not implemented'
       stop
    end select

  end function new_interpolator_1d

  ! Define spline interpolation of values in data define on original grid at points coordinates
  function interpolate_1d_spline(interp, num_points, coordinates, data) result(data_out)
    type(interpolator_1d), pointer, intent(in)       :: interp
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    
    ! compute the interpolating spline coefficients
    call compute_spline_1D( data, interp%spl%bc_type, interp%spl )
    call interpolate_array_values( coordinates, data_out, num_points, interp%spl )
    
  end function interpolate_1d_spline

  function interpolate_1d_weno(interp, num_points, coordinates, data) result(data_out)
    type(interpolator_1d), pointer, intent(in)       :: interp
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    
    data_out = interpolate_WENO_1D(interp%weno, num_points, coordinates, data)

  end function interpolate_1d_weno

end module sll_interpolator_1d
