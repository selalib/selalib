module geometry_functions
#include "sll_working_precision.h"

contains
  function default_x1 ( eta1, eta2 )
    sll_real64  :: default_x1
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    default_x1 = eta1
  end function default_x1

  function default_x2 ( eta1, eta2 )
    sll_real64  :: default_x2
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    default_x2 = eta2
  end function default_x2

  function default_jac11 ( eta1, eta2 )
    sll_real64  :: default_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    default_jac11 = 1.0_f64
  end function default_jac11

    function default_jac12 ( eta1, eta2 )
    sll_real64  :: default_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    default_jac12 = 0.0_f64
  end function default_jac12

  function default_jac21 ( eta1, eta2 )
    sll_real64  :: default_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    default_jac21 = 0.0_f64
  end function default_jac21

  function default_jac22 ( eta1, eta2 )
    sll_real64  :: default_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    default_jac22 = 1.0_f64
  end function default_jac22
end module geometry_functions
