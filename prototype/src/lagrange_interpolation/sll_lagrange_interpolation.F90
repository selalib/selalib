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

!le centrage est fait dans compute
!compute donne un tableau avec les wj, qui permettent de calculer le polynome
subroutine compute_lagrange_interpolation_1D(xi,lagrange)
type(sll_lagrange_interpolation_1D), pointer :: lagrange
sll_int32 :: i,j,k,bc_type,indice_decalage
sll_real64,dimension(1:lagrange%num_points),intent(in) :: xi
sll_real64,dimension(1:lagrange%d*2-1) :: xiadd
sll_real64,dimension(1:lagrange%d*2) :: xilocal
sll_real64 :: h
sll_real64,dimension(1:2*lagrange%d,1:lagrange%num_points) :: wj

h=xi(2)-xi(1)
if(lagrange%alpha<0)then
indice_decalage=lagrange%alpha/h-1
else
indice_decalage=lagrange%alpha/h
end if
!print*,"h = ",h,"  alpha =",lagrange%alpha,"  indice =",indice_decalage
!print*,""

bc_type=lagrange%bc_type
select case(bc_type)
case (PERIODIC_LAGRANGE)
 print*,"not implemented"
case (HERMITE_LAGRANGE)
print*,"Cas hermite"

!remplissage de xiadd, le tableau avec les noeuds fictifs supplémentaires
do i=1,lagrange%d-1
xiadd(i)=xi(1)
end do 
do i=lagrange%d,2*lagrange%d-1
xiadd(i)=xi(lagrange%num_points)
end do
print*,"xiadd"
print*,xiadd
print*,"fin de xiadd"

!remplissage du tableau wj
 wj=1.0_f64
 do k=1,lagrange%num_points
  do j=1,2*lagrange%d
   !creation du wj
   if(j+indicice_decalage<d)then !on prend des indices à gauche pour faire les wj
    do i=1,lagrange%d*2
    if(j+indice_decalage<1)then !on prend xi(j) à gauche du tableau xi
     !xi(j)=xiadd(debut)
     wj(j,k)=wj(j,k)
    else if(j+indice_decalage>num_points) !on prend xi(j) à droite du tableau xi
     !xi(j)=xiadd(lagrange%d+j+indice_decalage)   simple : xi(j)=xiadd(fin)
     wj(j,k)=wj(j,k)
    else !on prend xi(j) dans le tableau xi
     !xi(j)=xi(j) 
     wj(j,k)=wj(j,k)
    end if
    end do
   elseif(j+indice_decalage+d>num_points)then !on prend des indices à droite pour faire les wj

    if(j+indice_decalage<1)then !on prend xi(j) à gauche du tableau xi

    else if(j+indice_decalage>num_points) !on prend xi(j) à droite du tableau xi

    else !on prend xi(j) dans le tableau xi

    end if

   else !tous les indices sont dans le tableau xi

    wj(j,k)=wj(j,k)*(xi(j+indice_decalage)-xi(j+indice_decalage+1))
    do i=2,lagrange%d
     wj(j,k)=wj(j,k)*(xi(j+indice_decalage)-xi(j+indice_decalage+1-i))
     wj(j,k)=wj(j,k)*(xi(j+indice_decalage)-xi(j+indice_decalage+i)
    end do

   end if
   wj(j,k)=1.0_f64/wj(j,k)
   !fin de création du wj
  end do
 end do
case default
   print *, 'ERROR: compute_lagrange_interpolation_1D(): not recognized boundary condition'
   STOP
end select
print*,wj(:,1)
lagrange%wj=wj


!   do i=1,2*lagrange%d
!    if(i/=j) then
!     wj(j,k)=wj(j,k)*(xi(j)-xi(i))
!    end if
!   end do
!   wj(j,k)=1.0_f64/wj(j,k)


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
 print*,"not implemented"
case (HERMITE_LAGRANGE)
print*,"cas hermite pour interpolation"
 do k=1,lagrange%num_points
  do j=1,lagrange%d+1
   sum1=0.0_f64
   sum2=0.0_f64
   do i=1,lagrange%d
    sum1=sum1+lagrange%fi(i)*lagrange%wj(i,k)/(x(j)-lagrange%xi(i))
    sum2=sum2+lagrange%wj(i,k)/(x(j)-lagrange%xi(i))
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

