program lagrange_test
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
 use sll_lagrange_interpolation
implicit none
sll_int32  :: i,j,degree
sll_real64 :: res,diff
sll_real64,dimension(:),allocatable ::xi,fi,coord
type(sll_lagrange_interpolation_1D),pointer ::l_i

degree=3
allocate(xi(1:degree+1))
allocate(fi(1:degree+1))
allocate(coord(1:degree+2))
!data initialization
coord(1)=1.3_f64
coord(2)=2.5_f64
coord(3)=3.2_f64
coord(4)=4.5_f64
coord(5)=0.32_f64
do i=1,4
 xi(i)=i
 fi(i)=f(xi(i))
end do 
diff=0.0_f64
l_i => new_lagrange_interpolation_1D(xi,fi,degree,5,real(0,8))
call compute_lagrange_interpolation_1D(xi,l_i)
call interpolate_array_values(coord,l_i)
do i=1,5
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
f=3*x*x+2*x+1.0_f64
end function

end program

