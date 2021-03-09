program test_wedge_product_generaldim
#include "sll_working_precision.h"

   use sll_m_wedge_product_generaldim

   implicit none

   sll_real64, dimension(2, 1) :: r, s
   sll_real64, dimension(2, 1) :: t

   sll_real64, dimension(3, 1) :: m, n
   sll_real64, dimension(3, 1) :: o

   r(:, 1) = [sll_real64 :: 1, 2]
   s(:, 1) = [sll_real64 :: 3, 4]
   t(:, :) = sll_f_cross_product_2d(r, s)

   if (r(1, 1) /= 1.0_f64 .and. r(2, 1) /= 2.0_f64) stop 'FAILED'
   m(:, 1) = [sll_real64 :: 1, 2, 3]
   n(:, 1) = [sll_real64 :: 4, 5, 6]
   o(:, :) = sll_f_cross_product_3d(m, n)

   if (o(1, 1) /= -3.0_f64 .and. o(2, 1) /= 6.0_f64 .and. o(3, 1) /= -3.0_f64) stop 'FAILED'

   print *, 'PASSED'

end program test_wedge_product_generaldim
