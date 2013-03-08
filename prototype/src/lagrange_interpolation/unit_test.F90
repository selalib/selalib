program lagrange_test
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
 use sll_lagrange_interpolation
implicit none
sll_int32  :: i,j,d,num_points
sll_real64 :: res,diff,alpha
sll_real64,dimension(:),allocatable ::xi,fi,coord
type(sll_lagrange_interpolation_1D),pointer ::l_i

d=4
num_points=50
alpha=0.2
allocate(xi(1:num_points))
allocate(fi(1:num_points))
allocate(coord(1:num_points))
!data initialization
do i=1,num_points
 xi(i)=i-1
 fi(i)=f(xi(i))
 coord(i)=alpha+i-1
end do 
diff=0.0_f64

!test de l'indice en fonction de alpha
!l_i => new_lagrange_interpolation_1D(xi,fi,d,num_points,alpha,HERMITE_LAGRANGE)
l_i => new_lagrange_interpolation_1D(xi,fi,d,num_points,alpha,PERIODIC_LAGRANGE)
call compute_lagrange_interpolation_1D(xi,l_i)
call interpolate_array_values(coord,l_i)
do i=1,num_points
 print*,"interpolated value = ", l_i%data_out(i), " , Correct value = ",f(coord(i))
 diff=max(diff,abs(f(coord(i))-l_i%data_out(i)))
end do

if(diff<1e-10) then
 print *, ' '
 print *, 'Lagrange interpolation unit test: PASSED'
 print *, ' '
else
 print *, ' '
 print *, 'Lagrange interpolation unit test: FAILED'
 print *, ' '
end if

deallocate(xi)
deallocate(fi)
deallocate(coord)

contains 

function f(x)
sll_real64 :: x,f
f=2*x+1.0_f64
end function

end program

