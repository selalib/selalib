module test_function_module
#include "sll_working_precision.h"
implicit none

contains

  function test_func(x)
    intrinsic :: cos
    sll_real64 :: test_func
    sll_real64, intent(in) :: x
    test_func = x*x*cos(x)
  end function test_func

  function one(x)
    sll_real64 :: one
    sll_real64, intent(in) :: x
    one = 1.0_f64
  end function one

  function test_func_2d(x,y)
    intrinsic  :: cos
    intrinsic  :: sin
    sll_real64 :: test_func_2d
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    test_func_2d = x*y*cos(x)*sin(y)
  end function test_func_2d

  function one_2d(x, y)
    sll_real64 :: one_2d
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    one_2d = 1.0_f64
  end function one_2d


end module test_function_module
