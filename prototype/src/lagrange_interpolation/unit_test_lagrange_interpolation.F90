program test_lagrange
#include "sll_working_precision.h"
 use sll_lagrange_interpolation
implicit none
sll_int32 :: i,j
sll_real64 :: res,diff
sll_real64,dimension(1:4) ::xi,yi,coord,coef

coord(1)=1.0_f64
coord(2)=2.5_f64
coord(3)=3.2_f64
coord(4)=4.5_f64
diff=0.0_f64
do i=1,4
 xi(i)=i
 yi(i)=f(xi(i))
end do
do i=1,4
 res=lagrange_interpolation(coord(i),xi,yi,3)
 print*,"interpolated value = ", res, " , Correct value = ",f(coord(i))
 diff=diff+(f(coord(i))-res)
end do

if(diff<1e-5) then
 print*,"Lagrange interpolation is OK"
else
 print*,"Error of lagrange interpolation"
end if

diff=0.0_f64
coef=lagrange_interpolation(xi,yi,3)
do i=1,4
 res=calcul_newton(coord(i),xi,coef,3)
 print*,"interpolated value = ", res, " , Correct value = ",f(coord(i))
 diff=diff+(f(coord(i))-res)
end do

if(diff<1e-5) then
 print*,"Newton Lagrange interpolation is OK"
else
 print*,"Error of Newton lagrange interpolation"
end if

contains 

function f(x)
sll_real64 :: x,f
f=3*x*x+2*x+1.0_f64
end function

function testf(x)
sll_real64 :: x,testf
testf=1+x*x
end function

end program
