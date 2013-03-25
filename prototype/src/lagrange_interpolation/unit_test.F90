program lagrange_test
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
 use sll_lagrange_interpolation
 use numeric_constants
 use sll_lagrange_interpolation2
implicit none
sll_int32  :: i,j,d,num_points
sll_real64 :: res,diff,diff2,alpha,xmin,xmax,l
sll_real64,dimension(:),allocatable ::xi,fi,coord
type(sll_lagrange_interpolation_1D),pointer ::l_i
type(sll_lagrange_interpolation_1D2),pointer ::l_i2

d=2
num_points=5
alpha=-2.2_f64
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
 !fi(i)=f2(xi(i),l)
 coord(i)=xi(i)+alpha
end do 
diff=0.0_f64
diff2=0.0_f64

!l_i => new_lagrange_interpolation_1D(num_points,xmin,xmax,HERMITE_LAGRANGE,d)
l_i => new_lagrange_interpolation_1D(num_points,xmin,xmax,0,d)!PERIODIC_LAGRANGE,d)
l_i2 => new_lagrange_interpolation_1D2(num_points,xmin,xmax,0,d)!PERIODIC_LAGRANGE,d)
call compute_lagrange_interpolation_1D(alpha,l_i)
call compute_lagrange_interpolation_1D2(fi,alpha,l_i2)
call interpolate_array_values(fi,l_i)
call interpolate_array_values2(l_i2)
do i=1,num_points
 print*,"1 interpolated value = ", l_i%data_out(i), " , Correct value = ",f(coord(i),num_points)
 print*,"2 interpolated value = ", l_i2%data_out(i), " , Correct value = ",f(coord(i),num_points)
 diff=max(diff,abs(f(coord(i),num_points)-l_i%data_out(i)))
 diff2=max(diff2,abs(f(coord(i),num_points)-l_i2%data_out(i)))
end do

if(diff<0.1) then
 print *, ""
 print *, "Lagrange interpolation unit test: PASSED"
 print *, "error = ",diff,diff2
else
 print *, ""
 print *, "Lagrange interpolation unit test: FAILED"
 print *, "error = ",diff,diff2
end if

deallocate(xi)
deallocate(fi)
deallocate(coord)

contains 

function f(x,num_points)
sll_int32 :: num_points
sll_real64 :: x,f
!f=x*x+2*x
f=cos(2*sll_pi*x/(num_points-1.0))
end function

function f2(x,l)
sll_real64 :: x,f2,l
f2=1.0_8 / (2 + sin(6._f64*sll_pi/l*x))
end function

end program

