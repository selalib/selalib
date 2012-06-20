!> \brief abstract data type for 1D interpolation and reconstruction
!> 

module sll_module_interpolators_1d_base
#include "sll_working_precision.h"
  implicit none

  type, abstract :: sll_interpolator_1d_base
   contains
     ! must delete deferred keyword for the moment because sll_WENO
     ! (sll_WENO_node_interpolation.F90) is an extends
     ! type of interpolator_1d_base and i don't want touch this module.
     procedure(interpolator_1d_array_msg), deferred, pass(interpolator) :: &
          compute_interpolants
     procedure(interpolator_one_arg_msg), deferred, pass(interpolator) :: &
          interpolate_value
     procedure(interpolator_one_arg_msg), deferred, pass(interpolator) :: &
          interpolate_derivative_eta1
     
     procedure(interpolate_1d_array), pass, deferred :: interpolate_array
     procedure(reconstruct_1d_array), pass, deferred :: reconstruct_array
  end type sll_interpolator_1d_base

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
     subroutine interpolator_1d_array_msg( interpolator, data_array )
       use sll_working_precision
       import :: sll_interpolator_1d_base
       class(sll_interpolator_1d_base), intent(inout) :: interpolator
       sll_real64, dimension(:), intent(in) :: data_array
     end subroutine interpolator_1d_array_msg
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


end module sll_module_interpolators_1d_base
