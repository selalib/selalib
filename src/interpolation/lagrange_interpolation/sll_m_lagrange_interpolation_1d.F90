module sll_m_lagrange_interpolation_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_m_boundary_condition_descriptors

implicit none

 type :: sll_lagrange_interpolation_1D
   sll_int32                            :: d !half of stencil
   sll_int32                            :: nb_cell 
   sll_int32                            :: bc_type
   sll_int32                            :: index_gap
   sll_real64                           :: alpha
   sll_real64                           :: xmin
   sll_real64                           :: xmax
   sll_real64                           :: deta !< \a deta is the grid spacing
   sll_real64, dimension(:), pointer  :: wj
   sll_real64, dimension(:), pointer  :: wj_scale
   sll_real64, dimension(:), pointer    :: data_out !result=p(x) where p is the polynomial of interpolation
   sll_int32                            :: periodic_last !< \a periodic_last indicates if the input data repeats the first point at the end if we have periodic data. It takes the values 0 (not repeated) or 1 (repeated).  Default : 1.
 end type sll_lagrange_interpolation_1D

! integer, parameter :: SLL_PERIODIC = 0, SLL_HERMITE = 3

interface delete
  module procedure delete_lagrange_interpolation_1D
end interface

contains  !*****************************************************************************


function new_lagrange_interpolation_1D(num_points,xmin,xmax,bc_type,d, periodic_last)
 type(sll_lagrange_interpolation_1D), pointer :: new_lagrange_interpolation_1D
 sll_int32 ::ierr
 sll_int32,intent(in) :: d, num_points,bc_type
 sll_int32, intent(in), optional :: periodic_last !< \a periodic_last indicates if the input data repeats the first point at the end if we have periodic data. It takes the values 0 (not repeated) or 1 (repeated).  Default : 1.
 sll_real64 :: xmin,xmax
 
 SLL_ALLOCATE( new_lagrange_interpolation_1D, ierr )
 
 
 SLL_ALLOCATE(new_lagrange_interpolation_1D%wj(2*d),ierr)
 SLL_ALLOCATE(new_lagrange_interpolation_1D%wj_scale(2*d),ierr)
 SLL_ALLOCATE(new_lagrange_interpolation_1D%data_out(num_points),ierr)
 new_lagrange_interpolation_1D%d=d
 new_lagrange_interpolation_1D%xmin=xmin
 new_lagrange_interpolation_1D%xmax=xmax
 new_lagrange_interpolation_1D%nb_cell=num_points-1
 new_lagrange_interpolation_1D%bc_type=bc_type
 new_lagrange_interpolation_1D%deta = (xmax-xmin)/(new_lagrange_interpolation_1D%nb_cell)
 if (present(periodic_last)) then
    new_lagrange_interpolation_1D%periodic_last = periodic_last
 else
    new_lagrange_interpolation_1D%periodic_last = 1
 end if


end function new_lagrange_interpolation_1D


!> This function computes the weights w_j for the barycentric formula
subroutine compute_lagrange_interpolation_1D(lagrange)
 type(sll_lagrange_interpolation_1D), pointer :: lagrange !< \a lagrange is the lagrange interpolator object
 sll_int32 :: i,j
 sll_int32,dimension(1:4*lagrange%d-2) :: table
 sll_real64,dimension(1:2*lagrange%d) :: wj
 
 
 do i=1,2*lagrange%d-1
 table(i)=2*lagrange%d-1-(i-1)
 table(i+2*lagrange%d-1)=i
 end do
 
 wj=1.0_f64
 do i=1,lagrange%d
  do j=1,2*lagrange%d-1
   wj(i)=wj(i)*table(i+j-1)
  end do
  wj(i)=((-1.0_f64)**(lagrange%d+i))*wj(i)
 end do
 do i=1,lagrange%d
  wj(i+lagrange%d)=-wj(lagrange%d-i+1)
 end do
 wj=1.0_f64/wj 
 
 lagrange%wj=wj

end subroutine compute_lagrange_interpolation_1D

!> This function computes the value at the grid points displacement by -alpha
subroutine interpolate_from_interpolant_array(fi,alpha,lagrange)
type(sll_lagrange_interpolation_1D), pointer :: lagrange !< \a lagrange is the lagrange interpolator object
sll_real64, intent(in) :: alpha !< \a alpha is the (negative) displacement
sll_int32 ::i,j,index_gap
sll_real64 :: sum1,sum2,beta,h
sll_real64,dimension(1:lagrange%nb_cell+1),intent(in) :: fi !< \a fi are the values at the grid points

lagrange%alpha=-alpha

h= lagrange%deta! grid size
! How many cells do we displace? alpha = index_gap*h+beta*h
index_gap=floor(lagrange%alpha/h)
beta =  lagrange%alpha/h-real(index_gap,f64)

! Take care of the case where alpha/h is negative and below machine precision
if (beta == 1.0_f64) then
   beta = 0.0_f64
   index_gap = index_gap+1
end if

! If the displacement is precisely a multiple of h, we need to avoid division by zero
if (beta==0.0_f64) then
   select case(lagrange%bc_type)
   case (SLL_PERIODIC)
      do j=1,lagrange%nb_cell+lagrange%periodic_last
         lagrange%data_out(j) = fi(modulo(index_gap+j-1,lagrange%nb_cell)+1);
      end do
   case(SLL_HERMITE)
      do j=1,lagrange%nb_cell+1
         lagrange%data_out(j) = fi(min(max(1,index_gap+j),lagrange%nb_cell+1));
      end do
   end select

else

sum2=0.0_f64
do j=1,2*lagrange%d
 lagrange%wj_scale(j)=lagrange%wj(j)/(beta+real(lagrange%d-j,f64))
 sum2=sum2+lagrange%wj_scale(j)
end do

select case(lagrange%bc_type)
case (SLL_PERIODIC)
 
 do i=1,lagrange%nb_cell+lagrange%periodic_last
 sum1=0.0_f64
  do j=1,lagrange%d*2
   sum1=sum1+lagrange%wj_scale(j)*fi(modulo(index_gap+(i-1)+(j-1)-(lagrange%d-1),lagrange%nb_cell)+1)
  end do
 lagrange%data_out(i)=sum1/sum2
 end do

case(SLL_HERMITE)

 do i=1,lagrange%nb_cell+1
 sum1=0.0_f64
  do j=1,lagrange%d*2
   if(index_gap+(i-1)+(j-1)-(lagrange%d-1)<0)then
    sum1=sum1+lagrange%wj_scale(j)*fi(1)
   else if(index_gap+(i-1)+(j-1)-(lagrange%d-1) > lagrange%nb_cell)then
    sum1=sum1+lagrange%wj_scale(j)*fi(lagrange%nb_cell+1)
   else
    sum1=sum1+lagrange%wj_scale(j)*fi(index_gap+(i-1)+(j-1)-(lagrange%d-1)+1)
   end if
 end do
 lagrange%data_out(i)=sum1/sum2
 end do

case default
  print *, 'ERROR: compute_lagrange_interpolation_1D(): not recognized boundary condition'
  STOP
end select
end if

end subroutine interpolate_from_interpolant_array 
 
 
subroutine delete_lagrange_interpolation_1D( sll_m_lagrange_interpolation )
 type(sll_lagrange_interpolation_1D), pointer :: sll_m_lagrange_interpolation
 sll_int32                    :: ierr
   SLL_ASSERT( associated(sll_m_lagrange_interpolation) )
   SLL_DEALLOCATE( sll_m_lagrange_interpolation%wj, ierr )
   SLL_DEALLOCATE( sll_m_lagrange_interpolation, ierr )
end subroutine delete_lagrange_interpolation_1D

end module sll_m_lagrange_interpolation_1d

