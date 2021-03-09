program test_primitives
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
   use sll_m_constants, only: sll_p_pi
   use sll_m_primitives

   implicit none

   sll_real64, dimension(:), allocatable :: f
   sll_real64, dimension(:), allocatable :: x
   sll_int32                             :: i
   sll_int32                             :: n
   sll_real64                            :: xm

   n = 128

   allocate (f(n + 1))
   allocate (x(n + 1))

   do i = 1, n + 1
      x(i) = real(i - 1, f64)*(2.0_f64*sll_p_pi)/real(n, f64)
      f(i) = cos(x(i))
   end do

   call sll_s_function_to_primitive(f, x, n, xm)

   if (abs(xm) > 1d-7) stop 'FAILED'

   if (sum((sin(x) - f)) + sll_p_pi > 1d-7) stop 'FAILED'

   call sll_s_primitive_to_function(f, x, n, xm)

   if (sum(abs(cos(x) - f)) > 1d-7) stop 'FAILED'

   print *, 'PASSED'

end program test_primitives
