program test_fekete

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_box_splines, only: sll_s_write_connectivity

   use sll_m_fekete_integration, only: &
      sll_s_fekete_order_num, &
      sll_f_fekete_points_and_weights

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   sll_int32  :: i, j, n
   sll_int32  :: ierr
   sll_real64 :: s

   character(len=18) :: string

   sll_real64, dimension(2, 3) :: pxy1
   sll_real64, dimension(2, 3) :: pxy2
   sll_real64, dimension(:, :), allocatable :: xyw
   sll_int32  :: rule

   write (*, "(/,a)") "*********************************** "
   write (*, "(a)") "       FEKETE QUAD TEST       "
   write (*, "(a)") "*********************************** "

!Definition of first triangle
   pxy1(:, 1) = (/0._f64, 0._f64/)
   pxy1(:, 2) = (/1._f64, 0._f64/)
   pxy1(:, 3) = (/0._f64, 1._f64/)

!Definition of first triangle
   pxy2(:, 1) = (/1._f64, 0._f64/)
   pxy2(:, 2) = (/1._f64, 1._f64/)
   pxy2(:, 3) = (/0._f64, 1._f64/)

   rule = 2
   call sll_s_fekete_order_num(rule, n)
   SLL_CLEAR_ALLOCATE(xyw(1:3, 1:n), ierr)

   write (*, "(a)") " Computing Fekete points and weights on reference triangle "
   write (*, "(/,a)") "           x                   y                    w"
   xyw = sll_f_fekete_points_and_weights(pxy1, rule)

   write (string, '( "(",I2,"f20.15)" )') n
   do j = 1, n
      write (*, string) (xyw(i, j), i=1, 3)
   end do

   s = sum(xyw(3, :))
   print *, "sum weights = ", s

   if (abs(s - 1.0_f64) > 1d-7) then
      print *, 'FAILED'
   else
      print *, 'PASSED'
   end if

end program test_fekete
