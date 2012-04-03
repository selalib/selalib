!> \brief abstract data type for 1D interpolation and reconstruction
!> 

module sll_interpolator_1d
#include "sll_working_precision.h"
  implicit none

  type, abstract :: interpolator_1d_base
   contains
     procedure(interpolate_1d_array), pass, deferred :: interpolate_array
     procedure(reconstruct_1d_array), pass, deferred :: reconstruct_array
  end type interpolator_1d_base

  abstract interface
     function interpolate_1d_array(this, num_points, data, coordinates) result(res)
       use sll_working_precision
       import interpolator_1d_base
       class(interpolator_1d_base), intent(in)     :: this
       sll_int32, intent(in)                 :: num_points    ! size of output array
       sll_real64, dimension(:), intent(in)  :: data          ! data to be interpolated 
       sll_real64, dimension(:), intent(in)  :: coordinates   ! points where output is desired
       sll_real64, dimension(num_points)     :: res
     end function interpolate_1d_array
  end interface
  abstract interface
     function reconstruct_1d_array(this, num_points, data) result(res)
       use sll_working_precision
       import interpolator_1d_base
       class(interpolator_1d_base), intent(in)     :: this
       sll_int32, intent(in)                 :: num_points    ! size of output array
       sll_real64, dimension(:), intent(in)  :: data          ! data to be interpolated 
       sll_real64, dimension(num_points)     :: res
     end function reconstruct_1d_array
  end interface


 end module sll_interpolator_1d
