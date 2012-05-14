module sll_interpolators_base
#include "sll_working_precision.h" 
  implicit none
  
  !*************************************************************************
  !
  !                          2D Interpolators
  !
  !*************************************************************************
  
  ! Base class/basic interface for 2D interpolators
  
  ! TO BE RESOLVED:
  ! Function names should be reviewed and improved. What is the best way to
  ! express that a derivative is in a particular direction? Why eta???
  type, abstract :: interpolator_2d_base
   contains
     procedure(interpolator_array_msg), deferred, pass(interpolator) :: &
          compute_interpolants
     procedure(interpolator_two_arg_msg), deferred, pass(interpolator) :: &
          interpolate_value
     procedure(interpolator_two_arg_msg), deferred, pass(interpolator) :: &
          interpolate_derivative_eta1
     procedure(interpolator_two_arg_msg), deferred, pass(interpolator) :: &
          interpolate_derivative_eta2
  end type interpolator_2d_base
  
  ! Signature of the interpolating function  
  abstract interface
     function interpolator_two_arg_msg( interpolator, eta1, eta2 ) result(val)
       use sll_working_precision
       import :: interpolator_2d_base
       sll_real64                  :: val
       class(interpolator_2d_base), intent(in) :: interpolator
       sll_real64, intent(in)      :: eta1
       sll_real64, intent(in)      :: eta2
     end function interpolator_two_arg_msg
  end interface

  abstract interface
     subroutine interpolator_array_msg( interpolator, data_array )
       use sll_working_precision
       import :: interpolator_2d_base
       class(interpolator_2d_base), intent(inout) :: interpolator
       sll_real64, dimension(:,:), intent(in) :: data_array
     end subroutine interpolator_array_msg
  end interface
  
  
end module sll_interpolators_base
