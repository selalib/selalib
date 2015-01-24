program lagrange_test
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
 use sll_lagrange_interpolation
 use sll_constants
implicit none
sll_int32  :: i,d,num_points
sll_real64 :: diff,alpha,xmin,xmax,l
sll_real64,dimension(:),allocatable ::xi,fi,coord
type(sll_lagrange_interpolation_1D),pointer ::l_i

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
 xi(i)=i-1
 fi(i)=f(xi(i),num_points)
 coord(i)=xi(i)+alpha
end do 
diff=0.0_f64

l_i => new_lagrange_interpolation_1D(num_points,xmin,xmax,PERIODIC_LAGRANGE,d)
call compute_lagrange_interpolation_1D(alpha,l_i)
call interpolate_array_values(fi,l_i)
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

function f(x,num_points)
sll_int32 :: num_points
sll_real64 :: x,f
f=cos(2*sll_pi*x/(num_points-1.0))
end function

end program

