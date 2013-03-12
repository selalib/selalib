program lagrange_test
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
 use sll_lagrange_interpolation
 use numeric_constants
implicit none
sll_int32  :: i,j,d,num_points
sll_real64 :: res,diff,alpha,xmin,xmax,l
sll_real64,dimension(:),allocatable ::xi,fi,coord
type(sll_lagrange_interpolation_1D),pointer ::l_i

print*,sll_pi

d=6
num_points=256
alpha=0.2
allocate(xi(1:num_points))
allocate(fi(1:num_points))
allocate(coord(1:num_points))
!data initialization
xmin=0.0_f64
xmax=num_points-1.0_f64
l=xmax-xmin
do i=1,num_points
 xi(i)=i-1
 fi(i)=f(xi(i))
 coord(i)=xi(i)+alpha
end do 
diff=0.0_f64

!l_i => new_lagrange_interpolation_1D(num_points,xmin,xmax,HERMITE_LAGRANGE,d)
l_i => new_lagrange_interpolation_1D(num_points,xmin,xmax,PERIODIC_LAGRANGE,d)
call compute_lagrange_interpolation_1D(fi,alpha,l_i)
call interpolate_array_values(l_i)
do i=1,num_points
 !print*,"interpolated value = ", l_i%data_out(i), " , Correct value = ",f(coord(i))
 diff=max(diff,abs(f(coord(i))-l_i%data_out(i)))
end do

if(diff<0.1) then
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

function f(x)
sll_real64 :: x,f
f=sin(x+1.0_f64)+x*x
end function

function f2(x,l)
sll_real64 :: x,f2,l
f2=1.0_8 / (2 + sin(6._f64*sll_pi/l*x))
end function

end program

