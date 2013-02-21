module sll_lagrange_interpolation
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

implicit none

 type :: sll_lagrange_interpolation_1D
   sll_int32                            :: d
   sll_int32                            :: num_points
   sll_int32                            :: bc_type
   sll_real64                           :: alpha
   sll_real64, dimension(:,:), pointer  :: wj 
   sll_real64, dimension(:), pointer    :: xi 
   sll_real64, dimension(:), pointer    :: x
   sll_real64, dimension(:), pointer    :: fi ! fi(i)=function(xi(i))
   sll_real64, dimension(:), pointer    :: data_out !result=p(x) where p is the polynomial of interpolation
 end type sll_lagrange_interpolation_1D

 integer, parameter :: DEFAULT_LAGRANGE = 0, PERIODIC_LAGRANGE = 1, HERMITE_LAGRANGE = 2

interface delete
  module procedure delete_lagrange_interpolation_1D
end interface

contains  !*****************************************************************************


function new_lagrange_interpolation_1D(xi,fi,d,num_points,alpha,bc_type)
type(sll_lagrange_interpolation_1D), pointer :: new_lagrange_interpolation_1D
sll_int32 ::i,j,ierr
sll_int32,intent(in) :: d, num_points
sll_int32,intent(in),optional :: bc_type
sll_real64 :: alpha
sll_real64,dimension(1:num_points+1),intent(in) :: xi,fi

SLL_ALLOCATE( new_lagrange_interpolation_1D, ierr )

if(present(bc_type))then
 new_lagrange_interpolation_1D%bc_type=bc_type
else
 new_lagrange_interpolation_1D%bc_type=1
end if

SLL_ALLOCATE(new_lagrange_interpolation_1D%xi(num_points),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%fi(num_points),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%wj(d,num_points),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%x(num_points),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%data_out(num_points),ierr)
new_lagrange_interpolation_1D%d=d
new_lagrange_interpolation_1D%num_points=num_points
new_lagrange_interpolation_1D%xi=xi
new_lagrange_interpolation_1D%fi=fi
new_lagrange_interpolation_1D%alpha=alpha
end function new_lagrange_interpolation_1D


subroutine compute_lagrange_interpolation_1D(xi,lagrange)
type(sll_lagrange_interpolation_1D), pointer :: lagrange
sll_int32 :: i,j,k,bc_type
sll_real64,dimension(1:lagrange%num_points+1),intent(in) :: xi
sll_real64,dimension(1:lagrange%d+1,1:lagrange%num_points+1) :: wj

bc_type=lagrange%bc_type

select case(bc_type)


case (PERIODIC_LAGRANGE)
 print*,"pas encore de periodique"
case (HERMITE_LAGRANGE)
 wj=1.0_f64
 do k=1,lagrange%num_points+1
  do j=1,lagrange%d+1
   do i=1,lagrange%d+1
    if(i/=j) then
     wj(k,j)=wj(k,j)*(xi(j)-xi(i))
    end if
   end do
   wj(k,j)=1/wj(k,j)
  end do
 end do
case default
   print *, 'ERROR: compute_lagrange_interpolation_1D(): not recognized boundary condition'
   STOP
end select

lagrange%wj=wj
end subroutine compute_lagrange_interpolation_1D


subroutine interpolate_array_values(x,lagrange)
type(sll_lagrange_interpolation_1D), pointer :: lagrange
sll_int32 ::i,j,k,bc_type
sll_real64,dimension(1:lagrange%num_points), intent(in) :: x
sll_real64 :: sum1,sum2

bc_type=lagrange%bc_type
lagrange%x=x

select case(bc_type)
case (PERIODIC_LAGRANGE)
 print*,"pas encore implementer"
case (HERMITE_LAGRANGE)
 do k=1,lagrange%num_points+1
  do j=1,lagrange%d+1
   sum1=0.0_f64
   sum2=0.0_f64
   do i=1,lagrange%d+1
    sum1=sum1+lagrange%fi(i)*lagrange%wj(k,i)/(x(j)-lagrange%xi(i))
    sum2=sum2+lagrange%wj(k,i)/(x(j)-lagrange%xi(i))
   end do
   lagrange%data_out(j)=sum1/sum2
  end do
 end do
case default
   print *, 'ERROR: compute_lagrange_interpolation_1D(): not recognized boundary condition'
   STOP
end select
end subroutine interpolate_array_values 


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

