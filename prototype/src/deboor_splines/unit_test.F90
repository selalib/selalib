program test_deboor_splines
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_constants.h"
use sll_module_deboor_splines_1d

implicit none

sll_int32  :: n
sll_int32  :: n_der
sll_int32  :: deg_spline
sll_int32  :: ordre
sll_int32  :: ierr
sll_real64 :: d
sll_real64 :: res
sll_int32  :: i
sll_int32  :: j
sll_real64 :: tstart
sll_real64 :: tend

sll_real64, dimension(:), allocatable :: knots
sll_real64, dimension(:), allocatable :: f
sll_real64, dimension(:), allocatable :: df
sll_real64, dimension(:), allocatable :: x,x_test
sll_int32,  dimension(:), allocatable :: x_der
sll_real64, dimension(:), allocatable :: bcoef

deg_spline = 3
n          = 64
n_der      = 2
ordre      = deg_spline+1

SLL_ALLOCATE(knots(n+ordre+n_der),ierr)
SLL_ALLOCATE(f(n),ierr)
SLL_ALLOCATE(df(n_der),ierr)
SLL_ALLOCATE(x(n),ierr)
SLL_ALLOCATE(x_test(n),ierr)
SLL_ALLOCATE(x_der(n_der),ierr)
SLL_ALLOCATE(bcoef(n+n_der),ierr)

call cpu_time(tstart)

do j = 1, 1

  !Data positions and values
  do i = 1 , n
    x(i) = (i-1.)/(n-1)
    f(i) = cos(2*sll_pi*x(i))
  end do

  d = sum(abs(f(2:n)-f(1:n-1)))
  x_test(1) =  x(1)
  x_test(n) =  x(n)
  
  do i = 2,n
    x_test(i) = x_test(i-1) +  abs(f(i)-f(i-1))/d
  end do

  knots(1:ordre) = x(1)

  do i = ordre+1, n+deg_spline-1
    knots(i) = x(i-deg_spline)
  end do
  
  do i = n+deg_spline,n+ordre+n_der
    knots(i) = x(n)
  end do

  x_der(1) = 1
  df(1)    = -sin(2*sll_pi*x(1))*2*sll_pi
  x_der(2) = n
  df(2)    = -sin(2*sll_pi*x(n))*2*sll_pi

  call splint_der( x,      &
                   f,      &
                   x_der,  &
                   df,     &
                   knots,  &
                   n,      &
                   n_der,  &
                   ordre,  &
                   bcoef)

  do i = 1 , n
    res=bvalue( knots, bcoef, n+n_der, ordre, x(i), 0)
    write(10,*) x(i), x_test(i), f(i), res
    res=bvalue( knots, bcoef, n+n_der, ordre, x(i), 1)
    write(11,*) x(i), x_test(i), -sin(2*sll_pi*x(i))*2*sll_pi, res
  end do

end do

call cpu_time(tend)

print*, tend - tstart

DEALLOCATE(knots)
DEALLOCATE(f)
DEALLOCATE(df)
DEALLOCATE(x)
DEALLOCATE(x_der)
DEALLOCATE(bcoef)

print*, 'PASSED'
end program test_deboor_splines
