module user_geometry_functions
#include "sll_working_precision.h"
use sll_constants
  implicit none
contains
#define MAKE_AFFINE_MAP(f_name, eta, a, b) \
  function f_name(eta1,eta2); \
    sll_real64, intent(in) :: eta1, eta2; \
    sll_real64 :: f_name; \
    f_name = ((b)-(a))*eta + (a); \
  end function f_name

#define MAKE_LINEAR_MAP(f_name, eta, a) \
  function f_name(eta1,eta2); \
    sll_real64, intent(in) :: eta1, eta2; \
    sll_real64 :: f_name; \
    f_name = (a)*eta; \
  end function f_name

#define XMIN  0.0_f64
#define XMAX  (4*sll_pi)
#define VMIN  (-6.0_f64)
#define VMAX  6.0_f64  

  MAKE_AFFINE_MAP(x1_cartesian, eta1, XMIN, XMAX)
  MAKE_AFFINE_MAP(x2_cartesian, eta2, VMIN, VMAX)
  MAKE_LINEAR_MAP(jac1_cartesian, eta1, XMAX-XMIN)
  MAKE_LINEAR_MAP(jac2_cartesian, eta2, VMAX-VMIN)


end module user_geometry_functions
