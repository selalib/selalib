module test_function_module
#include "sll_working_precision.h"
  implicit none

contains

  function test_func(x)
    intrinsic :: dcos
    sll_real64 :: test_func
    sll_real64, intent(in) :: x
    test_func = x*x*dcos(x)
  end function test_func

  function one(x)
    sll_real64 :: one
    sll_real64, intent(in) :: x
    one = 1.0_f64
  end function one

end module test_function_module
