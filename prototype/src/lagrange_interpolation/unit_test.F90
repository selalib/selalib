program test_lagrange
#include "sll_working_precision.h"
 use sll_lagrange_interpolation
implicit none
sll_int32 :: i,j
sll_real64 :: res,diff
sll_real64,dimension(1:4) ::xi,fi,coord,wj


coord(1)=1.3_f64
coord(2)=2.5_f64
coord(3)=3.2_f64
coord(4)=4.5_f64
do i=1,4
 xi(i)=i
 fi(i)=f(xi(i))
end do 

!indirect
diff=0.0_f64
wj=lagrange_interpolation_wj(xi,3)
do i=1,4
 res=compute_lagrange_interpolation(coord(i),xi,fi,wj,3)
 print*,"interpolated value = ", res, " , Correct value = ",f(coord(i))
 diff=diff+abs(f(coord(i))-res)
end do

if(diff<1e-5) then
 print*,"Lagrange interpolation indirect is OK"
else
 print*,"Error : lagrange interpolation indirect"
end if

!direct
diff=0.0_f64
do i=1,4
 res=lagrange_interpolation_direct(coord(i),xi,fi,3)
 print*,"interpolated value = ", res, " , Correct value = ",f(coord(i))
 diff=diff+abs(f(coord(i))-res)
end do

if(diff<1e-5) then
 print*,"Lagrange interpolation direct is OK"
else
 print*,"Error : lagrange interpolation direct"
end if


contains 

function f(x)
sll_real64 :: x,f
f=3*x*x+2*x+1.0_f64
end function

end program

