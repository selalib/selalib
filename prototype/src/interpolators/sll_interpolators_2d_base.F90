module sll_module_interpolators_2d_base
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
  type, abstract :: sll_interpolator_2d_base
   contains
     procedure(interpolator_2d_array_msg), deferred, pass(interpolator) :: &
          compute_interpolants
     procedure(interpolator_two_arg_msg), deferred, pass(interpolator) :: &
          interpolate_value
     procedure(interpolator_two_arg_msg), deferred, pass(interpolator) :: &
          interpolate_derivative_eta1
     procedure(interpolator_two_arg_msg), deferred, pass(interpolator) :: &
          interpolate_derivative_eta2
     procedure(interpolate_2d_array), pass, deferred :: interpolate_array
     procedure(interpolate_2d_array_disp), pass, deferred :: interpolate_array_disp
  end type sll_interpolator_2d_base
  
  ! Signature of the interpolating function  
  abstract interface
     function interpolator_two_arg_msg( interpolator, eta1, eta2 ) result(val)
       use sll_working_precision
       import :: sll_interpolator_2d_base
       sll_real64                  :: val
       class(sll_interpolator_2d_base), intent(in) :: interpolator
       sll_real64, intent(in)      :: eta1
       sll_real64, intent(in)      :: eta2
     end function interpolator_two_arg_msg
  end interface

  abstract interface
     subroutine interpolator_2d_array_msg( interpolator, data_array )
       use sll_working_precision
       import :: sll_interpolator_2d_base
       class(sll_interpolator_2d_base), intent(inout) :: interpolator
       sll_real64, dimension(:,:), intent(in) :: data_array
     end subroutine interpolator_2d_array_msg
  end interface

  abstract interface
     function interpolate_2d_array(this,num_points1,num_points2,data_in,eta1,eta2) &
          result(res)
       use sll_working_precision
       import sll_interpolator_2d_base
       class(sll_interpolator_2d_base), intent(in)     :: this
       sll_int32, intent(in)  :: num_points1 
       sll_int32, intent(in)  :: num_points2 
       sll_real64, dimension(:,:), intent(in) :: data_in
       sll_real64, dimension(:,:), intent(in) :: eta1, eta2  
       sll_real64, dimension(num_points1,num_points2)    :: res
     end function interpolate_2d_array
  end interface

  abstract interface
     function interpolate_2d_array_disp(this,num_points1,num_points2,data_in,alpha1,alpha2) &
          result(res)
       use sll_working_precision
       import sll_interpolator_2d_base
       class(sll_interpolator_2d_base), intent(in)     :: this
       sll_int32, intent(in)  :: num_points1  
       sll_int32, intent(in)  :: num_points2 
       sll_real64, dimension(:,:), intent(in) :: data_in
       sll_real64, dimension(:,:), intent(in) :: alpha1, alpha2  
       sll_real64, dimension(num_points1,num_points2)    :: res
     end function interpolate_2d_array_disp
  end interface

end module sll_module_interpolators_2d_base
