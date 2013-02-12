module sll_lagrange_interpolation
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

implicit none

 type :: sll_lagrange_interpolation_1D
   sll_real64                           :: result !result=p(x) where p is the polynomial of interpolation
   sll_real64                           :: x
   sll_int32                            :: degree
   sll_int32                            :: bc_type
   sll_real64,dimension(:),pointer      :: wj 
   sll_real64,dimension(:),pointer      :: xi 
   sll_real64,dimension(:),pointer      :: fi ! fi(i)=function(xi(i))
 end type sll_lagrange_interpolation_1D

 integer, parameter :: PERIODIC_SPLINE = 0, HERMITE_SPLINE = 1

interface delete
  module procedure delete_lagrange_interpolation_1D
end interface

contains  !*****************************************************************************


function new_lagrange_interpolation_1D(xi,fi,degree)
type(sll_lagrange_interpolation_1D), pointer :: new_lagrange_interpolation_1D
sll_int32 ::i,j,ierr
sll_int32,intent(in) :: degree
sll_real64,dimension(1:degree+1),intent(in) :: xi,fi
sll_real64,dimension(1:degree+1) :: wj

SLL_ALLOCATE( new_lagrange_interpolation_1D, ierr )

wj=1.0_f64
do j=1,degree+1
 do i=1,degree+1
  if(i/=j) then
   wj(j)=wj(j)*(xi(j)-xi(i))
  end if
 end do
 wj(j)=1/wj(j)
end do

SLL_ALLOCATE(new_lagrange_interpolation_1D%xi(degree),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%fi(degree),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%wj(degree),ierr)
new_lagrange_interpolation_1D%degree=degree
new_lagrange_interpolation_1D%xi=xi
new_lagrange_interpolation_1D%fi=fi
new_lagrange_interpolation_1D%wj=wj
end function new_lagrange_interpolation_1D

subroutine compute_lagrange_interpolation_1D(x,lagrange_interpolation)
type(sll_lagrange_interpolation_1D), pointer :: lagrange_interpolation
sll_int32 ::i
sll_real64, intent(in) :: x
sll_real64 :: sum1,sum2,result

lagrange_interpolation%x=x
lagrange_interpolation%result=result

sum1=0.0_f64
sum2=0.0_f64
do i=1,lagrange_interpolation%degree+1
 sum1=sum1+lagrange_interpolation%fi(i)*lagrange_interpolation%wj(i)/(x-lagrange_interpolation%xi(i))
 sum2=sum2+lagrange_interpolation%wj(i)/(x-lagrange_interpolation%xi(i))
end do
lagrange_interpolation%result=sum1/sum2
end subroutine compute_lagrange_interpolation_1D

subroutine delete_lagrange_interpolation_1D( lagrange_interpolation )
type(sll_lagrange_interpolation_1D), pointer :: lagrange_interpolation
sll_int32                    :: ierr
  SLL_ASSERT( associated(lagrange_interpolation) )
  SLL_DEALLOCATE( lagrange_interpolation%xi, ierr )
  SLL_DEALLOCATE( lagrange_interpolation%fi, ierr )
  SLL_DEALLOCATE( lagrange_interpolation%wj, ierr )
  SLL_DEALLOCATE( lagrange_interpolation, ierr )
end subroutine delete_lagrange_interpolation_1D

end module

