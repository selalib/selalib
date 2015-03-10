program lagrange_test
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use sll_lagrange_interpolation
use sll_constants

implicit none

sll_int32  :: i
sll_int32  :: ierr
sll_int32  :: d
sll_int32  :: n
sll_real64 :: diff
sll_real64 :: alpha
sll_real64 :: xmin
sll_real64 :: xmax
sll_real64 :: dx

sll_real64, dimension(:), allocatable :: xi
sll_real64, dimension(:), allocatable :: yi
sll_real64, dimension(:), allocatable :: xp
sll_real64, dimension(:), allocatable :: yp

type(sll_lagrange_interpolation_1D), pointer :: l_i

d     = 2
n     = 100
alpha = 0.2_f64

SLL_ALLOCATE(xi(1:n), ierr)
SLL_ALLOCATE(yi(1:n), ierr)
SLL_ALLOCATE(xp(1:n), ierr)
SLL_ALLOCATE(yp(1:n), ierr)

!data initialization
xmin = 0.0_f64
xmax = 1.0_f64
dx   = (xmax - xmin) / (n-1)

do i=1,n
 xi(i) = (i-1)*dx
 yi(i) = f(xi(i))
 xp(i) = xi(i)+alpha*dx
 yp(i) = f(xp(i))
end do

diff = 0.0_f64

l_i => new_lagrange_interpolation_1D(n,                 &
                                     xmin,              &
                                     xmax,              &
                                     PERIODIC_LAGRANGE, &
                                     d)

call compute_lagrange_interpolation_1D(alpha,l_i)

call interpolate_array_values(yi,l_i)

diff = maxval(abs(yp-l_i%data_out))

if(diff<0.0001) then
 print *, ""
 print *, "Lagrange interpolation unit test: PASSED"
 print *, "error = ",diff
else
 print *, ""
 print *, "Lagrange interpolation unit test: FAILED"
 print *, "error = ",diff
end if

deallocate(xi)
deallocate(yi)
deallocate(xp)
deallocate(yp)

contains 

function f(x)

  sll_real64 :: f
  sll_real64 :: x
  
  f=cos(2*sll_pi*x)

end function f

end program

