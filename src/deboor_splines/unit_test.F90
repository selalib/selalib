program unit_test
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_constants.h"
use sll_module_deboor_splines_1d

implicit none

sll_real64, dimension(:), pointer :: knots
sll_real64, dimension(:), allocatable :: value_func
sll_real64, dimension(:), allocatable :: value_der
sll_real64, dimension(:), allocatable :: x,points_test
sll_int32,  dimension(:), allocatable :: points_der
sll_int32                         :: n
sll_int32                         :: m
sll_int32                         :: k
sll_real64, dimension(:), allocatable :: matrices
sll_real64, dimension(:), pointer :: bcoef
sll_int32                         :: iflag
sll_int32                         :: ierr
sll_real64                        :: d
sll_real64                        :: res
sll_int32                         :: i

k = 3
n = 10
m = 2

SLL_ALLOCATE(knots(n+(k+1)+m),ierr)
SLL_ALLOCATE(value_func(n),ierr)
SLL_ALLOCATE(value_der(m),ierr)
SLL_ALLOCATE(x(n),ierr)
SLL_ALLOCATE(points_test(n),ierr)
SLL_ALLOCATE(points_der(m),ierr)
SLL_ALLOCATE(matrices((2*(k+1)-1)*(n+m)),ierr)
SLL_ALLOCATE(bcoef(n+m),ierr)

!PN : Set point values and abscissae
do i = 1 , n
 x(i) = 0.0_f64 + (i-1)*1./(n-1)
 value_func(i) = cos(2*sll_pi*x(i))
end do

d = sum( abs( value_func(2:n)-value_func(1:n-1)))
points_test(1) =  x(1)
points_test(n) =  x(n)

do i = 2,n
  points_test(i) = points_test(i-1) +  abs( value_func(i)-value_func(i-1))/d
end do

!PN : Initialize knots
knots(1:(k+1)) = x(1)
do i = (k+1)+1, n + k-1
  knots(i) = x(i-k)
end do
do i = n +k,n + (k+1) + m
  knots(i) = x(n)
end do

points_der(1) = 1
value_der(1)  = -sin(2*sll_pi*x(1))*2*sll_pi
points_der(2) = n
value_der(2)  = -sin(2*sll_pi*x(n))*2*sll_pi

call splint_der(&
     x,value_func,&
     points_der,value_der,&
     knots,n,m,(k+1),&
     matrices, bcoef, iflag )

!PN : interpolate values
print*, " Values "
do i = 1 , n
  res=bvalue( knots, bcoef, n+m, (k+1), x(i), 0)
  print"(3(a8,f17.12))", 'approx',res,'exact', value_func(i), 'error', res-value_func(i)
end do
!PN : interpolate values of first derivative
!print*, " Values of first derivative "
!do i = 1,n
!  res=bvalue( knots, bcoef, n+m, (k+1), x(i), 1)
!  print"(a8,f17.12,a8,f17.12)", 'approx',res,'exact', value_der(i)
!end do
   
print*, 'PASSED'
end program unit_test
