module sll_lagrange_interpolation
#include "sll_working_precision.h"
implicit none

contains

function lagrange_interpolation(x,xi,yi,degree)
sll_int32 ::i,j
sll_int32, intent(in) :: degree
sll_real64 :: lagrange_interpolation,li
sll_real64,intent(in) :: x
sll_real64,dimension(1:degree+1),intent(in) ::xi,yi

li=1.0_f64
lagrange_interpolation=0.0_f64
do i=1,degree+1
 do j=1,degree+1
  if(j/=i)then
   li=li*(x-xi(j))/(xi(i)-xi(j))
  end if
 end do
 lagrange_interpolation=lagrange_interpolation+yi(i)*li
 li=1.0_f64
end do
end function

! function neuville_lagrange_interpolation(x,xi,yi,degree,epsilon)
! sll_int32 ::i,j
! sll_int32, intent(in) :: degree
! sll_real64 :: lagrange_interpolation,li
! sll_real64,intent(in) :: x,epsilon
! sll_real64,dimension(1:degree+1),intent(in) ::xi,yi
! 
! end function

end module 
