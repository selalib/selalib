!> \brief Abstract data type for 1D interpolation
!>
!> 


module sll_interpolator1d_interface
#include "sll_working_precision.h"

  implicit none
  
  type, abstract :: interpolator1d
   contains
     procedure(realarray_function_interpolator1d_int_realarray), deferred, pass:: interpolate1d
  end type interpolator1d
  
  abstract interface
     function realarray_function_interpolator1d_int_realarray(this, num_points, data) result(res)
       use sll_working_precision
       class(interpolator1d), intent(in)     :: this
       sll_int32, intent(in)                 :: num_points    ! number of interpolation points
       sll_real64, dimension(:), intent(in)  :: data          ! data to be interpolated 
       sll_real64, dimension(num_points)     :: res
     end function realarray_function_interpolator1d_int_realarray
  end interface


end module sll_interpolator1d_interface
