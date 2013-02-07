program test_lagrange_1d
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

!!!test Lagrange interpolation 1d
do i=1,4
 xi(i)=i
 yi(i)=f(xi(i))
end do
do i=1,4
 res=lagrange_interpolation(coord(i),xi,yi,3)
 print*,"interpolated value = ", res, " , Correct value = ",f(coord(i))
 diff=diff+abs(f(coord(i))-res)
end do

if(diff<1e-5) then
 print*,"Lagrange interpolation 1d is OK"
else
 print*,"Error : lagrange interpolation 1d"
end if

!!!test Newton Lagrange interpolation
diff=0.0_f64
coef=lagrange_interpolation(xi,yi,3)
do i=1,4
 res=compute_newton_interpolation(coord(i),xi,coef,3)
 print*,"interpolated value = ", res, " , Correct value = ",f(coord(i))
 diff=diff+abs(f(coord(i))-res)
end do

if(diff<1e-5) then
 print*,"Newton Lagrange interpolation is OK"
else
 print*,"Error : Newton lagrange interpolation"
end if

contains 

function f(x)
sll_real64 :: x,f
f=3*x*x+2*x+1.0_f64
end function

end program
