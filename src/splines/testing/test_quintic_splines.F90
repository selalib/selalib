program test_quintic_splines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_quintic_splines, only: &
    inspl5, &
    splin5

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
sll_real64           :: err(4)

print*, 'Test quintic splines low level function'

dx = (x_max-x_min)/(n-1)

do i = 1, n
  x(i)    = x_min + (i-1)*dx
  cf(1,i) = g(x(i))
  cf(2,i) = 0.0_f64 
  cf(3,i) = 0.0_f64 
end do

ind1 = 1  ! value function and first derivative are given
indn = 1
cf(2,1) = dg(x_min)
cf(2,n) = dg(x_max)
cf(3,1) = ddg(x_min)
cf(3,n) = ddg(x_max)

call inspl5(n,x,ind1,indn,cf,h)
  
err(:) = 0.0_f64
do i = 1, 1000
  
  xx = x_min + (i-1)/999.*(x_max-x_min)
  call splin5(n,x,cf,xx,0,f(1))
  call splin5(n,x,cf,xx,1,f(2))
  call splin5(n,x,cf,xx,2,f(3))
  
  err(1) = err(1) + ( g(xx)   - f(1))**2 
  err(2) = err(2) + ( dg(xx)  - f(2))**2 
  err(3) = err(3) + ( ddg(xx) - f(3))**2 
  write(48,*) xx, f(1)-g(xx)
  
end do

print"(' error on interpolated value ',f25.20)", sqrt(err(1))
print"(' error on first derivative   ',f25.20)", sqrt(err(2))
print"(' error on second derivative  ',f25.20)", sqrt(err(3))

print*, 'PASSED'

contains

function g(x)

  sll_real64 :: x
  sll_real64 :: g

  g =  sin(x)

end function g

function dg(x)

  sll_real64 :: x
  sll_real64 :: dg

  dg = cos(x)

end function dg

function ddg(x)

  sll_real64 :: x
  sll_real64 :: ddg

  ddg = -sin(x)

end function ddg

end program test_quintic_splines
