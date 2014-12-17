program test_quintic_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_constants.h"

use sll_quintic_splines

sll_int32, parameter :: n = 100    ! number of interpolation points
sll_real64           :: x(n)       ! vector of abscissae
sll_real64           :: f(3,n)     ! interpolated values and its derivatives
sll_real64           :: xx(n)      ! interpolated values abscissae
sll_real64           :: cf(1:3,n)  ! ordinates, first derivatibes, second derivatives
sll_int32            :: ind1, indn ! boundary conditions switches at x(1) and x(n)
                                   ! for both switches one has the correspondance
                                   ! = 1 type 1
                                   ! = 0 type 2
                                   ! =-1 type 3
sll_real64           :: h(6*n-3)   ! auxiliary vector

sll_int32            :: i
sll_real64           :: dx
sll_real64           :: x_min = -2
sll_real64           :: x_max = +2
sll_real64           :: err(3)

print*, 'Test quintic splines low level function'

dx = (x_max-x_min)/(n-1)
do i = 1, n
  x(i) = x_min + (i-1)*dx
  cf(1,i) =  x(i)**3-2*x(i)**2-x(i)+5
end do

call inspl5(n,x,ind1,indn,cf,h)

call random_number(xx)
xx = x_min + xx * (x_max-x_min)

err(1) = 0.0_f64
err(2) = 0.0_f64
err(3) = 0.0_f64
open(33, file='quintic.dat')
do i = 1, n

  call splin5(n,x,cf,xx(i),f(:,i))

  err(1) = err(1) + ( g(xx(i)) - f(1,i))**2 * dx
  err(2) = err(2) + ( dg(xx(i)) - f(2,i))**2 * dx
  err(3) = err(3) + ( ddg(xx(i)) - f(3,i))**2 * dx
  write(33,*) xx(i), f(1:3,i)

end do
close(33)

print"(' error on interpolated value ',f25.20)", sqrt(err(1))
print"(' error on first derivative   ',f25.20)", sqrt(err(2))
print"(' error on second derivative  ',f25.20)", sqrt(err(3))

print*, 'PASSED'

contains

function g(x)
sll_real64 :: x
sll_real64 :: g
  g =  x**3-2*x**2-x+5
end function g
function dg(x)
sll_real64 :: x
sll_real64 :: dg
  dg =  3*x**2-4*x-1
end function dg
function ddg(x)
sll_real64 :: x
sll_real64 :: ddg
  ddg =  6*x-4
end function ddg

end program test_quintic_splines
