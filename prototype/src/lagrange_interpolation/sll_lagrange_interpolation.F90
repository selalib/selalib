module sll_lagrange_interpolation
#include "sll_working_precision.h"
implicit none

contains

!! compute parameter wj, for lagrange interpolation
function lagrange_interpolation_wj(xi,degree) result(wj)
sll_int32 ::i,j
sll_int32,intent(in) :: degree
sll_real64,dimension(1:degree+1),intent(in) :: xi
sll_real64,dimension(1:degree+1) :: wj

wj=1.0_f64
do j=1,degree+1
 do i=1,degree+1
  if(i/=j) then
   wj(j)=wj(j)*(xi(j)-xi(i))
  end if
 end do
 wj(j)=1/wj(j)
end do

end function

function compute_lagrange_interpolation(x,xi,fi,wj,degree) result(res)
sll_int32 ::i
sll_int32,intent(in) :: degree
sll_real64,dimension(1:degree+1),intent(in) :: xi,wj,fi
sll_real64, intent(in) :: x
sll_real64 :: res,sum1,sum2

sum1=0.0_f64
sum2=0.0_f64
do i=1,degree+1
 sum1=sum1+fi(i)*wj(i)/(x-xi(i))
 sum2=sum2+wj(i)/(x-xi(i))
end do
res=sum1/sum2
end function

function lagrange_interpolation_direct(x,xi,fi,degree) result(res)
sll_int32 ::i,j
sll_int32, intent(in) :: degree
sll_real64, dimension(1:degree+1),intent(in) :: xi,fi
sll_real64,dimension(1:degree+1) ::wj
sll_real64,intent(in) :: x
sll_real64 :: res,sum1,sum2


sum1=0.0_f64
sum2=0.0_f64
wj=1.0_f64

do j=1,degree+1
 do i=1,degree+1
  if(i/=j) then
   wj(j)=wj(j)*(xi(j)-xi(i))
  end if
 end do
 wj(j)=1/wj(j)
end do


do i=1,degree+1
 sum1=sum1+fi(i)*wj(i)/(x-xi(i))
 sum2=sum2+wj(i)/(x-xi(i))
end do
res=sum1/sum2

end function

end module

