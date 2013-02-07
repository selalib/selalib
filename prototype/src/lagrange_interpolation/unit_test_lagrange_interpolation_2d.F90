program test_lagrange_2d
#include "sll_working_precision.h"
 use sll_lagrange_interpolation
implicit none
sll_int32 :: i,j
sll_real64 :: res,diff,x,y
sll_real64,dimension(1:4) :: xi
sll_real64,dimension(1:5) :: yi
sll_real64,dimension(1:4,1:5) :: zi
sll_real64,dimension(1:9) ::coord

do i=1,4
 xi(i)=-1.0_f64+(i-1)*0.5_f64
end do
do j=1,5
 yi(j)=2.0_f64+(j-1)*0.5_f64
end do
do i=1,4
 do j=1,5
 zi(i,j)=f(xi(i),yi(j))
 end do
end do

coord(1)=1.0_f64
coord(2)=3.0_f64
coord(3)=4.0_f64
coord(4)=2.0_f64
coord(5)=8.0_f64
coord(6)=7.0_f64
coord(7)=6.0_f64
coord(8)=4.0_f64
coord(9)=6.0_f64
do i=1,4
 do j=1,5
  x=coord(i)
  y=coord(4+j)
  res=lagrange_interpolation(x,y,xi,yi,zi,3,4)
  print*,"interpolated value = ", res, " , Correct value = ",f(x,y)
  diff=diff+abs(f(x,y)-res)
 end do
end do

if(diff<1e-5) then
 print*,"Lagrange interpolation 2d is OK"
else
 print*,"Error : Lagrange interpolation 2d "
end if


contains

function f(x,y)
sll_real64 ::f,x,y
f=x*x+2.0_f64*x*y+y+3.0_f64

end function

end program
