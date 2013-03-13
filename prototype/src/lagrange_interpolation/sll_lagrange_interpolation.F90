module sll_lagrange_interpolation
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

implicit none

 type :: sll_lagrange_interpolation_1D
   sll_int32                            :: d
   sll_int32                            :: num_points
   sll_int32                            :: bc_type
   sll_int32                            :: indice_decalage
   sll_real64                           :: alpha
   sll_real64                           :: xmin
   sll_real64                           :: xmax
   sll_real64, dimension(:,:), pointer  :: wj
   sll_real64, dimension(:), pointer    :: xiadd 
   sll_real64, dimension(:), pointer    :: fiadd 
   sll_real64, dimension(:), pointer    :: data_out !result=p(x) where p is the polynomial of interpolation
 end type sll_lagrange_interpolation_1D

 integer, parameter :: PERIODIC_LAGRANGE = 0, HERMITE_LAGRANGE = 1

interface delete
  module procedure delete_lagrange_interpolation_1D
end interface

contains  !*****************************************************************************


function new_lagrange_interpolation_1D(num_points,xmin,xmax,bc_type,d)
type(sll_lagrange_interpolation_1D), pointer :: new_lagrange_interpolation_1D
sll_int32 ::i,j,ierr
sll_int32,intent(in) :: d, num_points
sll_int32,intent(in),optional :: bc_type
sll_real64 :: xmin,xmax

SLL_ALLOCATE( new_lagrange_interpolation_1D, ierr )


SLL_ALLOCATE(new_lagrange_interpolation_1D%xiadd(num_points+2*d-1),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%fiadd(num_points+2*d-1),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%wj(2*d,num_points),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%data_out(num_points),ierr)
new_lagrange_interpolation_1D%d=d
new_lagrange_interpolation_1D%xmin=xmin
new_lagrange_interpolation_1D%xmax=xmax
new_lagrange_interpolation_1D%num_points=num_points
end function new_lagrange_interpolation_1D

!compute_lagrange_interpolation_1D écrit les wj, qui permettent de calculer le polynome
subroutine compute_lagrange_interpolation_1D(fi,alpha,lagrange)
type(sll_lagrange_interpolation_1D), pointer :: lagrange
sll_int32 :: i,j,k,indice_decalage
sll_real64 :: xij,xii,h
sll_real64, intent(in) :: alpha
sll_real64,dimension(1:lagrange%num_points) :: xi
sll_real64,dimension(1:lagrange%num_points),intent(in) :: fi
sll_real64,dimension(1:lagrange%d*2-1+lagrange%num_points) :: xiadd,fiadd
sll_real64,dimension(1:2*lagrange%d,1:lagrange%num_points) :: wj


h=(lagrange%xmax-lagrange%xmin)/(lagrange%num_points-1)
do i=0,lagrange%num_points-1
xi(i+1)=lagrange%xmin+i*h
end do
if(alpha<0)then
indice_decalage=alpha/h-1
else
indice_decalage=alpha/h
end if

select case(lagrange%bc_type)
case (PERIODIC_LAGRANGE)

  !remplissage de xiadd et fiadd en ajoutant les noeuds fictifs à xi
  do i=1,lagrange%num_points+2*lagrange%d-1
  xiadd(i)=xi(modulo(i-lagrange%d+indice_decalage,lagrange%num_points)+1)
  fiadd(i)=fi(modulo(i-lagrange%d+indice_decalage,lagrange%num_points)+1)
  end do

case (HERMITE_LAGRANGE)

  !remplissage de xiadd en ajoutant les noeuds fictifs à xi
  do i=1,lagrange%d-1
  xiadd(i)=xi(1)-(lagrange%d-i)*h+indice_decalage
  end do
  do i=lagrange%d,lagrange%d+lagrange%num_points-1
  xiadd(i)=xi(i-lagrange%d+1)+indice_decalage
  end do
  do i=lagrange%d+lagrange%num_points,lagrange%num_points+2*lagrange%d-1
  xiadd(i)=xi(lagrange%num_points)+(i-lagrange%d-lagrange%num_points+1)*h+indice_decalage
  end do
  !remplissage du tableau fiadd
  j=1
  do i=1,lagrange%num_points+2*lagrange%d-1
   if(xiadd(i).lt.xi(1))then
    fiadd(i)=fi(1)
   else if(xiadd(i).gt.xi(lagrange%num_points))then
    fiadd(i)=fi(lagrange%num_points)
   else
    fiadd(i)=fi(j)
    j=j+1
   end if
  end do

case default
  print *, 'ERROR: compute_lagrange_interpolation_1D(): not recognized boundary condition'
  STOP
end select

!remplissage du tableau wj
wj=1.0_f64
do k=1,lagrange%num_points
 do j=-lagrange%d+1,lagrange%d
  xij=xiadd(k+lagrange%d-1+j)
  do i=-lagrange%d+1,lagrange%d
   if(i/=j) then
    xii=xiadd(k+lagrange%d-1+i)
    wj(j+lagrange%d,k)=wj(j+lagrange%d,k)*(xij-xii)
   endif  
  end do
 wj(j+lagrange%d,k)=1.0_f64/wj(j+lagrange%d,k)
 end do
end do

lagrange%wj=wj
lagrange%xiadd=xiadd
lagrange%fiadd=fiadd
lagrange%alpha=alpha

end subroutine compute_lagrange_interpolation_1D


subroutine interpolate_array_values(lagrange)
type(sll_lagrange_interpolation_1D), pointer :: lagrange
sll_int32 ::j,k
sll_real64 :: sum1,sum2,xk,h
sll_real64, dimension(1:lagrange%num_points) :: x

h=(lagrange%xmax-lagrange%xmin)/(lagrange%num_points-1)

do k=1,lagrange%num_points
 sum1=0.0_f64
 sum2=0.0_f64
 xk=lagrange%xmin+(k-1)*h+lagrange%alpha
 do j=1,lagrange%d*2
  sum1=sum1+lagrange%fiadd(k-1+j)*lagrange%wj(j,k)/(xk-lagrange%xiadd(k-1+j))
  sum2=sum2+lagrange%wj(j,k)/(xk-lagrange%xiadd(k-1+j))
 end do
 lagrange%data_out(k)=sum1/sum2
end do

end subroutine interpolate_array_values 


subroutine delete_lagrange_interpolation_1D( lagrange_interpolation )
type(sll_lagrange_interpolation_1D), pointer :: lagrange_interpolation
sll_int32                    :: ierr
  SLL_ASSERT( associated(lagrange_interpolation) )
  SLL_DEALLOCATE( lagrange_interpolation%wj, ierr )
  SLL_DEALLOCATE( lagrange_interpolation, ierr )
end subroutine delete_lagrange_interpolation_1D

end module

