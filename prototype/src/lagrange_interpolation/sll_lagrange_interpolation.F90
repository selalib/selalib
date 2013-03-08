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
   sll_real64, dimension(:,:), pointer  :: wj
   sll_real64, dimension(:), pointer    :: xiadd 
   sll_real64, dimension(:), pointer    :: fiadd 
   sll_real64, dimension(:), pointer    :: xi 
   sll_real64, dimension(:), pointer    :: x
   sll_real64, dimension(:), pointer    :: fi ! fi(i)=function(xi(i))
   sll_real64, dimension(:), pointer    :: data_out !result=p(x) where p is the polynomial of interpolation
 end type sll_lagrange_interpolation_1D

 integer, parameter :: PERIODIC_LAGRANGE = 0, HERMITE_LAGRANGE = 1

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
 new_lagrange_interpolation_1D%bc_type=0
end if

SLL_ALLOCATE(new_lagrange_interpolation_1D%xi(num_points),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%xiadd(num_points+2*d-1),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%fiadd(num_points+2*d-1),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%fi(num_points),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%wj(2*d,num_points),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%x(num_points),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%data_out(num_points),ierr)
new_lagrange_interpolation_1D%d=d
new_lagrange_interpolation_1D%num_points=num_points
new_lagrange_interpolation_1D%xi=xi
new_lagrange_interpolation_1D%fi=fi
new_lagrange_interpolation_1D%alpha=alpha
end function new_lagrange_interpolation_1D

!compute donne un tableau avec les wj, qui permettent de calculer le polynome
subroutine compute_lagrange_interpolation_1D(xi,lagrange)
type(sll_lagrange_interpolation_1D), pointer :: lagrange
sll_int32 :: i,j,k,indice_decalage
sll_real64 :: xij,xii,h
sll_real64,dimension(1:lagrange%num_points),intent(in) :: xi
sll_real64,dimension(1:lagrange%d*2-1+lagrange%num_points) :: xiadd,fiadd
sll_real64,dimension(1:2*lagrange%d,1:lagrange%num_points) :: wj

!calcul de h
h=xi(2)-xi(1)
!calcul d'indice_decalage en fonction de alpha
if(lagrange%alpha<0)then
indice_decalage=lagrange%alpha/h-1
else
indice_decalage=lagrange%alpha/h
end if


select case(lagrange%bc_type)
case (PERIODIC_LAGRANGE)

  !remplissage de xiadd et fiadd en ajoutant les noeuds fictifs à xi
  do i=1,lagrange%num_points+2*lagrange%d-1
  xiadd(i)=xi(modulo(i-lagrange%d+indice_decalage,lagrange%num_points)+1)
  fiadd(i)=lagrange%fi(modulo(i-lagrange%d+indice_decalage,lagrange%num_points)+1)
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
    fiadd(i)=lagrange%fi(1)
   else if(xiadd(i).gt.xi(lagrange%num_points))then
    fiadd(i)=lagrange%fi(lagrange%num_points)
   else
    fiadd(i)=lagrange%fi(j)
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

end subroutine compute_lagrange_interpolation_1D


subroutine interpolate_array_values(x,lagrange)
type(sll_lagrange_interpolation_1D), pointer :: lagrange
sll_int32 ::j,k
sll_real64 :: sum1,sum2
sll_real64, dimension(1:lagrange%num_points), intent(in) :: x

do k=1,lagrange%num_points
 sum1=0.0_f64
 sum2=0.0_f64
 do j=1,lagrange%d*2
  sum1=sum1+lagrange%fiadd(k-1+j)*lagrange%wj(j,k)/(x(k)-lagrange%xiadd(k-1+j))
  sum2=sum2+lagrange%wj(j,k)/(x(k)-lagrange%xiadd(k-1+j))
 end do
 lagrange%data_out(k)=sum1/sum2
end do

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

