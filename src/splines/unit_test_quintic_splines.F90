program test_quintic_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_constants.h"

use sll_quintic_splines

sll_int32, parameter :: n = 64     ! number of interpolation points
sll_real64           :: x(n)       ! vector of abscissae
sll_real64           :: f(3)       ! interpolated values and its derivatives
sll_real64           :: xx         ! interpolated values abscissae
sll_real64           :: cf(1:3,n)  ! ordinates, first derivatibes, second derivatives
sll_int32            :: ind1, indn ! boundary conditions switches at x(1) and x(n)
                                   ! for both switches one has the correspondance
                                   ! = 1 type 1
                                   ! = 0 type 2
                                   ! =-1 type 3
sll_real64           :: h(6*n-3)   ! auxiliary vector

sll_int32            :: i
sll_real64           :: dx
sll_real64           :: x_min = 0.0_f64
sll_real64           :: x_max = 1.0_f64
sll_real64           :: err(3)

print*, 'Test quintic splines low level function'

dx = (x_max-x_min)/(n-1)
do i = 1, n
  x(i)    = x_min + (i-1)*dx
  cf(1,i) = g(x(i))
  cf(2,i) = dg(x(i))
  cf(3,i) = ddg(x(i))
end do

ind1 = -1
indn = -1
call inspl5(n,x,ind1,indn,cf,h)


open(33, file='quintic_0.dat')
open(34, file='quintic_1.dat')
open(35, file='quintic_2.dat')
err(:) = 0.0_f64
do i = 1, 1000

  xx = x_min + (i-1)/999.*(x_max-x_min)
  call splin5(n,x,cf,xx,f(1:3))

  err(1) = err(1) + ( g(xx)   - f(1))**2 
  err(2) = err(2) + ( dg(xx)  - f(2))**2 
  err(3) = err(3) + ( ddg(xx) - f(3))**2 

  write(33,*) xx, f(1), g  (xx)
  write(34,*) xx, f(2), dg (xx)
  write(35,*) xx, f(3), ddg(xx)

end do
close(33)
close(34)
close(35)

print"(' error on interpolated value ',f25.20)", sqrt(err(1))
print"(' error on first derivative   ',f25.20)", sqrt(err(2))
print"(' error on second derivative  ',f25.20)", sqrt(err(3))

print*, 'PASSED'

contains

function g(x)
sll_real64 :: x
sll_real64 :: g
  g =  0.25_f64*x*x*x*x
end function g
function dg(x)
sll_real64 :: x
sll_real64 :: dg
  dg = x*x*x
end function dg
function ddg(x)
sll_real64 :: x
sll_real64 :: ddg
  ddg = 3.0_f64*x*x
end function ddg

end program test_quintic_splines
