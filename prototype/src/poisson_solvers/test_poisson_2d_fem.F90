program test_poisson_2d_fem
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_constants.h"

use sll_poisson_2d_fem
implicit none

sll_int32  :: i, j
sll_real64 :: dimx, dimy
sll_int32  :: error
sll_real64 :: mode
sll_real64, dimension(:),   pointer :: x
sll_real64, dimension(:),   pointer :: y
sll_real64, dimension(:,:), pointer :: ex
sll_real64, dimension(:,:), pointer :: ey
sll_real64, dimension(:,:), pointer :: rho
sll_real64 :: errmax

call test()

contains

subroutine test()
type( poisson_fem ) :: poisson
sll_real64 :: dx, dy
sll_int32  :: nc_x, nc_y

nc_x = 64
nc_y = 64

SLL_CLEAR_ALLOCATE(ex(1:nc_x+1,1:nc_y+1),error)  
SLL_CLEAR_ALLOCATE(ey(1:nc_x+1,1:nc_y+1),error) 
SLL_CLEAR_ALLOCATE(rho(1:nc_x+1,1:nc_y+1),error)  

SLL_ALLOCATE(x(1:nc_x+1),error)  
SLL_ALLOCATE(y(1:nc_y+1),error) 

dimx = 1.0
dimy = 1.0

dx = dimx / nc_x
dy = dimy / nc_y

!Create an irregular mesh
do i=1,nc_x+1
   x(i) = (i-1)*dx * 0.5 * ((i-1)*dx+1)
enddo
do j=1,nc_y+1
   y(j) = (j-1)*dy * 0.5 * ((j-1)*dy+1)
enddo

mode = 4*sll_pi
do j = 1, nc_y
   do i = 1, nc_x
      rho(i,j) = 2_f64 * mode**2 * sin(mode*x(i))*sin(mode*y(j))
   end do
end do

call initialize(poisson, x, y, nc_x+1, nc_y+1)
call solve(poisson, ex, ey, rho)

errmax = 0.
do j = 2, nc_y
   do i = 2, nc_x
      errmax = max(errmax,abs(rho(i,j)-sin(mode*x(i))*sin(mode*y(j))))
   end do
end do
print*, 'error = ', errmax

do j = 1, nc_y+1
   do i = 1, nc_x+1
      write(11,*) x(i), y(j), rho(i,j), sin(mode*x(i))*sin(mode*y(j))
   end do
   write(11,*)
end do
close(11)


end subroutine test


end program test_poisson_2d_fem
