program test_lagrange_3d
#include "sll_working_precision.h"
 use sll_lagrange_interpolation
implicit none
sll_int32 :: i,j,k
sll_real64 :: res,diff,x,y,z
sll_real64,dimension(1:4) :: xi,zi
sll_real64,dimension(1:5) :: yi
sll_real64,dimension(1:4,1:5,1:4) :: fi
sll_real64,dimension(1:13) ::coord

do i=1,4
 xi(i)=-1.0_f64+(i-1)*0.5_f64
 zi(i)=-3.0_f64+(i-1)*2.0_f64
end do
do j=1,5
 yi(j)=2.0_f64+(j-1)*0.5_f64
end do
do i=1,4
 do j=1,5
  do k=1,4
   fi(i,j,k)=f(xi(i),yi(j),zi(k))
  end do
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
coord(10)=2.0_f64
coord(11)=3.0_f64
coord(12)=5.0_f64
coord(13)=1.0_f64
do i=1,4
 do j=1,5
  do k=1,4
   x=coord(i)
   y=coord(4+j)
   z=coord(9+k)
   res=lagrange_interpolation(x,y,z,xi,yi,zi,fi,3,4,3)
   print*,"interpolated value = ", res, " , Correct value = ",f(x,y,z)
   diff=diff+abs(f(x,y,z)-res)
  end do
 end do
end do

if(diff<1e-5) then
 print*,"Lagrange interpolation 3d is OK"
else
 print*,"Error : Lagrange interpolation 3d "
end if


contains

function f(x,y,z)
sll_real64 ::f,x,y,z
f=x*x+2.0_f64*x*y+y+3.0_f64+z*x+z*z*z

end function

end program

