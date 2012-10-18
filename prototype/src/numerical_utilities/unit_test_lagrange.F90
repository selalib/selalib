program test_lagrange
#include "sll_working_precision.h"
  use lagrange_interpolation
  implicit none

#define NP 4
#define XMIN  0.0_f64
#define DELTA 1.0_f64

  sll_int32 :: i
  sll_real64, dimension(1:NP) :: xi
  sll_real64, dimension(1:NP) :: yi
  sll_real64                  :: res
  sll_real64                  :: tmp
  sll_real64 :: p012, p123

  do i=1,NP
     xi(i) = XMIN + real(i-1,f64)*DELTA
!     yi(i) = test_f(xi(i))
     yi(i) = square(xi(i))
  end do

  ! We could add more points in between the nodes...
  do i=1, NP
     tmp = xi(i)
     res = lagrange_interpolate( tmp, 3, xi, yi )
     print *, 'interpolated value = ', res, '. Correct value = ', yi(i)
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
