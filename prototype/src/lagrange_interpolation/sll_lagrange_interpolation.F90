module sll_lagrange_interpolation
#include "sll_working_precision.h"
implicit none

  interface lagrange_interpolation
    module procedure lagrange_interpolation_1d, lagrange_interpolation_2d, &
                     neville_lagrange_interpolation, newton_lagrange_interpolation
  end interface

contains

function lagrange_interpolation_1d(x,xi,yi,degree) result (res)
sll_int32 ::i,j
sll_int32, intent(in) :: degree
sll_real64 :: li,res
sll_real64,intent(in) :: x
sll_real64,dimension(1:degree+1),intent(in) ::xi,yi

li=1.0_f64
res=0.0_f64
do i=1,degree+1
 do j=1,degree+1
  if(j/=i)then
   li=li*(x-xi(j))/(xi(i)-xi(j))
  end if
 end do
 res=res+yi(i)*li
 li=1.0_f64
end do
end function

function newton_lagrange_interpolation(xi,yi,degree) result(coef)
sll_int32 ::i,j
sll_int32, intent(in) :: degree
sll_real64,dimension(1:degree+1),intent(in) ::xi,yi
sll_real64,dimension(1:degree+1) :: coef

coef=yi
do i=1,degree+1
 do j=degree+1,i+1,-1
  coef(j)=(coef(j)-coef(j-1))/(xi(j)-xi(j-i))
 end do
end do

end function

function calcul_newton(x,xi,coef,degree) result(res)
sll_int32 :: i
sll_int32,intent(in) :: degree
sll_real64,intent(in) :: x
sll_real64 :: res
sll_real64,dimension(1:degree+1),intent(in) :: coef,xi
sll_real64,dimension(1:degree+1) :: base

base(1)=1.0_f64
res=coef(1)
do i=2,degree+1
 base(i)=base(i-1)*(x-xi(i-1))
 res=res+base(i)*coef(i)
end do

end function


! function newton_lagrange_interpolation(x,xi,yi,degree) result(res)
! sll_int32 ::i,j
! sll_int32, intent(in) :: degree
! sll_real64 :: res,li
! sll_real64,intent(in) :: x
! sll_real64,dimension(1:degree+1),intent(in) ::xi,yi
! sll_real64,dimension(1:degree+1) :: coef,base
! 
! coef=yi
! do i=1,degree+1
!  do j=degree+1,i+1,-1
!   coef(j)=(coef(j)-coef(j-1))/(xi(j)-xi(j-i))
!  end do
! end do
! 
! base(1)=1.0_f64
! do i=2,degree+1
!  base(i)=base(i-1)*(x-xi(i-1))
! end do
! 
! res=0
! do i=1,degree+1
!  res=res+base(i)*coef(i)
! end do
! print*,""
! print*,"coef",coef
! print*,"base",base
! 
! end function


function neville_lagrange_interpolation(x,xi,yi,degree,epsilon) result(res)
sll_int32 ::i,j
sll_int32, intent(in) :: degree
sll_real64 :: res,li
sll_real64,intent(in) :: x,epsilon
sll_real64,dimension(1:degree+1),intent(in) ::xi,yi

print*,"pas encore implementé"
end function


function lagrange_interpolation_2d(coord,xi,yi,zi,degree) result(res)
sll_int32 ::i,j,k
sll_int32, intent(in) :: degree
sll_real64 :: li,lj,res
sll_real64,dimension(2),intent(in) :: coord
sll_real64,dimension(1:degree+1),intent(in) ::xi,yi,zi

print*,"pas encore implementé"
end function

end module 
