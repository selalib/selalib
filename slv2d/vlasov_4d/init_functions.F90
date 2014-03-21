module init_functions
#include "sll_working_precision.h"
#include "sll_constants.h"

implicit none

integer, parameter :: LANDAU_X_CASE = 1
integer, parameter :: LANDAU_Y_CASE = 2
integer, parameter :: LANDAU_COS_PROD_CASE = 3
integer, parameter :: LANDAU_COS_SUM_CASE  = 4
integer, parameter :: TSI_CASE  = 5

integer, parameter :: VA_VALIS          = 0
integer, parameter :: VA_OLD_FUNCTION   = 1
integer, parameter :: VA_VLASOV_POISSON = 2

integer, parameter :: METH_BSL_CUBIC_SPLINES = 0
integer, parameter :: METH_CSL_CUBIC_SPLINES = 1
integer, parameter :: METH_CSL_LAG3 = 2
integer, parameter :: METH_CSL_PPM1 = 3
integer, parameter :: METH_CSL_PPM2 = 4


contains

  function landau_1d(eps,kx, x, v2)
    sll_real64 :: landau_1d
    sll_real64, intent(in) :: x, kx, v2, eps

    landau_1d = (1._f64+eps*cos(kx*x))/(2*sll_pi)*exp(-0.5_f64*v2)

  end function landau_1d

  function landau_cos_prod(eps,kx, ky, x, y, v2)
    sll_real64 :: landau_cos_prod
    sll_real64, intent(in) :: x, kx
    sll_real64, intent(in) :: y, ky
    sll_real64, intent(in) :: eps, v2

    landau_cos_prod = (1._f64+eps*cos(kx*x)*cos(ky*y))/(2*sll_pi)*exp(-0.5_f64*v2)

  end function landau_cos_prod

  function landau_cos_sum(eps,kx, ky, x, y, v2)
    sll_real64 :: landau_cos_sum
    sll_real64, intent(in) :: x, kx
    sll_real64, intent(in) :: y, ky
    sll_real64, intent(in) :: eps, v2

     landau_cos_sum = (1._f64+eps*cos(kx*(x+y)))/(2*sll_pi)*exp(-0.5_f64*v2)

  end function landau_cos_sum

  function tsi(eps, kx, x, vx, v2)
    sll_real64 :: tsi
    sll_real64, intent(in) :: x, kx, vx
    sll_real64, intent(in) :: eps, v2

    tsi = (1._f64+eps*cos(kx*x))/(2*sll_pi)*exp(-0.5_f64*v2)*vx*vx

  end function tsi

end module init_functions

