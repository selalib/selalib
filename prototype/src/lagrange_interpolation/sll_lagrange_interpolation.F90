module sll_lagrange_interpolation
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

implicit none

 type :: sll_lagrange_interpolation_1D
   sll_int32                            :: d !moitié du stencil
   sll_int32                            :: nb_cell 
   sll_int32                            :: bc_type
   sll_int32                            :: indice_decalage
   sll_real64                           :: alpha
   sll_real64                           :: xmin
   sll_real64                           :: xmax
   sll_real64, dimension(:), pointer  :: wj
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
sll_int32,intent(in) :: d, num_points,bc_type
sll_real64 :: xmin,xmax

SLL_ALLOCATE( new_lagrange_interpolation_1D, ierr )


SLL_ALLOCATE(new_lagrange_interpolation_1D%wj(2*d),ierr)
SLL_ALLOCATE(new_lagrange_interpolation_1D%data_out(num_points),ierr)
new_lagrange_interpolation_1D%d=d
new_lagrange_interpolation_1D%xmin=xmin
new_lagrange_interpolation_1D%xmax=xmax
new_lagrange_interpolation_1D%nb_cell=num_points-1
end function new_lagrange_interpolation_1D

!compute_lagrange_interpolation_1D écrit les wj, qui permettent de calculer le polynome
subroutine compute_lagrange_interpolation_1D(alpha,lagrange)
type(sll_lagrange_interpolation_1D), pointer :: lagrange
sll_int32 :: i,j,k,indice_decalage
sll_real64 :: h
sll_real64, intent(in) :: alpha
sll_int32,dimension(1:4*lagrange%d-2) :: tableautest
sll_real64,dimension(1:2*lagrange%d) :: wj


h=(lagrange%xmax-lagrange%xmin)/(lagrange%nb_cell)
if(alpha<0)then
indice_decalage=alpha/h-1
else
indice_decalage=alpha/h
end if


select case(lagrange%bc_type)
case (PERIODIC_LAGRANGE)

do i=i,2*lagrange%d-1
tableautest(i)=2*lagrange%d-1-(i-1)
tableautest(i+2*lagrange%d-1)=i
end do
!print*,"tableau test",tableautest

!remplissage de wj
wj(:)=1.0_f64
do i=1,lagrange%d
 do j=1,2*lagrange%d-1
  wj(i)=wj(i)*tableautest(i+j-1)
 end do
 wj(i)=((-1)**(lagrange%d+i))*(h**(lagrange%d-1))*wj(i)
end do
do i=1,lagrange%d
 wj(i+lagrange%d)=-wj(lagrange%d-i+1)
end do
wj=1.0_f64/wj

print*,"wj dans principale",wj

case (HERMITE_LAGRANGE)

print*,"not yet"

!  !remplissage de xiadd en ajoutant les noeuds fictifs à xi
!  do i=1,lagrange%d-1
!  xiadd(i)=xi(1)-(lagrange%d-i)*h+indice_decalage
!  end do
!  do i=lagrange%d,lagrange%d+lagrange%num_points-1
!  xiadd(i)=xi(i-lagrange%d+1)+indice_decalage
!  end do
!  do i=lagrange%d+lagrange%num_points,lagrange%num_points+2*lagrange%d-1
!  xiadd(i)=xi(lagrange%num_points)+(i-lagrange%d-lagrange%num_points+1)*h+indice_decalage
!  end do
!  !remplissage du tableau fiadd
!  j=1
!  do i=1,lagrange%num_points+2*lagrange%d-1
!   if(xiadd(i).lt.xi(1))then
!    fiadd(i)=fi(1)
!   else if(xiadd(i).gt.xi(lagrange%num_points))then
!    fiadd(i)=fi(lagrange%num_points)
!   else
!    fiadd(i)=fi(j)
!    j=j+1
!   end if
!  end do

case default
  print *, 'ERROR: compute_lagrange_interpolation_1D(): not recognized boundary condition'
  STOP
end select

lagrange%wj=wj
lagrange%alpha=alpha

end subroutine compute_lagrange_interpolation_1D


subroutine interpolate_array_values(fi,lagrange)
type(sll_lagrange_interpolation_1D), pointer :: lagrange
sll_int32 ::i,j,indice_decalage
sll_real64 :: sum1,sum2,beta,h
sll_real64,dimension(1:lagrange%nb_cell+1),intent(in) :: fi


select case(lagrange%bc_type)
case (PERIODIC_LAGRANGE)
h=(lagrange%xmax-lagrange%xmin)/(lagrange%nb_cell)
if(lagrange%alpha<0)then
indice_decalage=lagrange%alpha/h-1
beta=h+mod(lagrange%alpha,h)
else
indice_decalage=lagrange%alpha/h
beta=mod(lagrange%alpha,h)
end if

print*,"h",h,"alpha",lagrange%alpha,"beta",beta

do j=1,2*lagrange%d
 lagrange%wj(j)=lagrange%wj(j)/(beta+real(lagrange%d,f64)-real(j,f64))
end do

do i=1,lagrange%nb_cell+1
sum1=0.0_f64
sum2=0.0_f64
 do j=1,lagrange%d*2
  sum1=sum1+lagrange%wj(j)*fi(modulo(indice_decalage+(i-1)+(j-1)-(lagrange%d-1),lagrange%nb_cell)+1)
  sum2=sum2+lagrange%wj(j)
!print*,"dans 1 fi",fi(modulo(indice_decalage+(i-1)+(j-1)-(lagrange%d-1),lagrange%nb_cell)+1),"modulo",modulo(indice_decalage+(i-1)+(j-1)-(lagrange%d-1),lagrange%nb_cell)+1
!print*,"dans normal wj/bla",lagrange%wj(j),"sum1",sum1,"sum2",sum2
 end do
!print*,""
lagrange%data_out(i)=sum1/sum2
end do


! do i=1,lagrange%nb_cell+1
!  do j=1,lagrange%d*2
!   sum1=sum1+lagrange%wj(j)/(beta+lagrange%d-j)*fi(modulo(indice_decalage+i-1-(lagrange%d-1)+(j-1),lagrange%nb_cell+1)+1)
!   sum2=sum2+lagrange%wj(j)/(beta+lagrange%d-j)
!  end do
! lagrange%data_out(i)=sum1/sum2
! end do


case(HERMITE_LAGRANGE)
print*,"not yet"

! do k=1,lagrange%num_points
!  sum1=0.0_f64
!  sum2=0.0_f64
!  xk=lagrange%xmin+(k-1)*h+lagrange%alpha
!  do j=1,lagrange%d*2
!   sum1=sum1+lagrange%fiadd(k-1+j)*lagrange%wj(j,k)/(xk-lagrange%xiadd(k-1+j))
!   sum2=sum2+lagrange%wj(j,k)/(xk-lagrange%xiadd(k-1+j))
!  end do
!  lagrange%data_out(k)=sum1/sum2
! end do

case default
  print *, 'ERROR: compute_lagrange_interpolation_1D(): not recognized boundary condition'
  STOP
end select

end subroutine interpolate_array_values 
 
 
 subroutine delete_lagrange_interpolation_1D( lagrange_interpolation )
 type(sll_lagrange_interpolation_1D), pointer :: lagrange_interpolation
 sll_int32                    :: ierr
   SLL_ASSERT( associated(lagrange_interpolation) )
   SLL_DEALLOCATE( lagrange_interpolation%wj, ierr )
   SLL_DEALLOCATE( lagrange_interpolation, ierr )
 end subroutine delete_lagrange_interpolation_1D

end module

