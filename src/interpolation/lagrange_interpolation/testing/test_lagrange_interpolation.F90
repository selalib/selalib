program test_lagrange_interpolation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_lagrange_interpolation_1d, only: &
    sll_s_compute_lagrange_interpolation_1d, &
    sll_s_interpolate_from_interpolant_array, &
    sll_f_new_lagrange_interpolation_1d, &
    sll_t_lagrange_interpolation_1d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_int32  :: i,d,num_points
sll_real64 :: diff,alpha,xmin,xmax,l
sll_real64,dimension(:),allocatable ::xi,fi,coord
type(sll_t_lagrange_interpolation_1d),pointer ::l_i

d=2
num_points=100
alpha=0.2_f64
allocate(xi(1:num_points))
allocate(fi(1:num_points))
allocate(coord(1:num_points))
!data initialization
xmin=0.0_f64
xmax=num_points-1.0_f64
l=xmax-xmin
do i=1,num_points
 xi(i)=real(i-1,f64)
 fi(i)=f(xi(i),num_points)
 coord(i)=xi(i)+alpha
end do 
diff=0.0_f64

l_i => sll_f_new_lagrange_interpolation_1d(num_points,xmin,xmax,sll_p_periodic,d)
call sll_s_compute_lagrange_interpolation_1d(l_i)
call sll_s_interpolate_from_interpolant_array(fi,-alpha,l_i)
do i=1,num_points
 !print*,"interpolated value = ", l_i%data_out(i), " , Correct value = ",f(coord(i),num_points)
 diff=max(diff,abs(f(coord(i),num_points)-l_i%data_out(i)))
end do


if(diff<0.0001) then
 print *, ""
 print *, "Lagrange interpolation unit test: PASSED"
 print *, "error = ",diff
else
 print *, ""
 print *, "Lagrange interpolation unit test: FAILED"
 print *, "error = ",diff
end if

deallocate(xi)
deallocate(fi)
deallocate(coord)

contains 

function f( x, num_points )
  sll_real64, intent(in) :: x
  sll_int32,  intent(in) :: num_points
  sll_real64 :: f
  f = cos(2*sll_p_pi*x/(num_points-1.0))
end function

end program test_lagrange_interpolation

