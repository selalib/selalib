program test_lagrange
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic

   use sll_m_lagrange_interpolation, only: &
      sll_s_compact_derivative_weight, &
      sll_s_compute_stencil_plus, &
      sll_f_lagrange_interpolate, &
      sll_o_weight_product_x1

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define NP 20
#define XMIN  0.0_f64
#define DELTA 1.0_f64

   sll_int32 :: i
   sll_real64, dimension(1:NP) :: xi
   sll_real64, dimension(1:NP) :: yi
   sll_real64                  :: res
   sll_real64                  :: tmp
!  sll_real64 :: p012, p123
   !-----> for test of derivatives
   sll_real64 :: w(-20:20)
   !sll_int32 :: ierr
   sll_int32 :: p
   sll_int32 :: r, s

   do i = 1, NP
      xi(i) = XMIN + real(i - 1, f64)*DELTA
!     yi(i) = test_f(xi(i))
      yi(i) = square(xi(i))
   end do

   ! We could add more points in between the nodes...
   do i = 1, NP
      tmp = xi(i)
      res = sll_f_lagrange_interpolate(tmp, 3, xi, yi)
      print *, '#interpolated value = ', res, '. Correct value = ', yi(i)
   end do

   !----------->
   ! test of derivatives
   do p = 0, 5
      call sll_s_compute_stencil_plus(p, r, s)
      call sll_s_compact_derivative_weight(w(r:s), r, s)
      print *, r, s
      print *, w(r:s)
      call sll_o_weight_product_x1(yi, xi, NP - 1, w(r:s), r, s, sll_p_periodic)
      do i = 1, NP
         print *, '#derivative(i) = ', xi(i), 2._f64*(XMIN + real(i - 1, f64)*DELTA)
      end do
   end do

contains

   function test_f(x)
      sll_real64 :: test_f
      sll_real64, intent(in) :: x
      test_f = 1.0_f64 + 2.0_f64*x*x + x*x*x
   end function test_f

   function square(x)
      sll_real64 :: square
      sll_real64, intent(in) :: x
      square = x*x
   end function square
end program test_lagrange
