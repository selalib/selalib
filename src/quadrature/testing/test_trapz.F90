program test_trapz

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_constants, only: sll_p_pi
   use sll_m_trapz_integration, only: sll_o_trapz_integrate_1d

   use test_function_module

   implicit none

   sll_int32  :: i, j
   sll_real64 :: s, f
   sll_real64, dimension(10) :: x

   write (*, '(5x, 2a16 )') 'trapz', 'Exact value'
   do i = 2, 10
      do j = 1, i
         x(j) = real(j - 1, f64)*0.5_f64*sll_p_pi/real(i - 1, f64)
      end do
      write (*, '(a, i2, a, 5f16.12)') 'n = ', i, ': ', &
         sll_o_trapz_integrate_1d(test_func, x, i), &
         0.4674011002723395_f64
   end do

   write (*, *) " ------- "

   s = 0.0_f64
   do i = 2, 10
      do j = 1, i
         x(j) = real(j - 1, f64)*1.0_f64/real(i - 1, f64)
      end do
      f = sll_o_trapz_integrate_1d(one, x, i)
      write (*, '(a, i2, a, 4f16.12)') 'n = ', i, ': ', f
      s = s + f
   end do
   print *, 'Exact value: '
   write (*, '(f22.15)') 1.00000

   if (abs(s - 9.0_f64) > 1d-7) then
      print *, 'FAILED'
   else
      print *, 'PASSED'
   end if

end program test_trapz
