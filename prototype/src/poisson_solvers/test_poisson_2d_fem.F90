program test_poisson_2d_fem
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_constants.h"

#define nx 64
#define ny 64

use sll_poisson_2d_fem
implicit none

sll_int32  :: i, j
sll_real64 :: dimx, dimy
sll_int32  :: error, mode
sll_real64, dimension(:),   pointer :: x
sll_real64, dimension(:),   pointer :: y
sll_real64, dimension(:,:), pointer :: ex
sll_real64, dimension(:,:), pointer :: ey
sll_real64, dimension(:,:), pointer :: rho

call test()

contains

subroutine test()
type( poisson_fem ) :: poisson
sll_real64 :: dx, dy

ncx = 64
ncy = 64

SLL_CLEAR_ALLOCATE(ex(1:ncx+1),error)  
SLL_CLEAR_ALLOCATE(ey(1:ncy+1),error) 
SLL_CLEAR_ALLOCATE(rho(1:ncx+1),error)  

SLL_ALLOCATE(x(-1:nx+1),error)  
SLL_ALLOCATE(y(-1:ny+1),error) 

dimx = 2 * sll_pi
dimy = 2 * sll_pi

dx = dimx / ncx
dy = dimy / ncy

x(0) = 0.
y(0) = 0.

do i=1,ncx+1
   x(i) = (i*dx) *(i*dx+1)/(1+dimx)
enddo
do j=1,ncy+1
   y(j) = (j*dy) *(j*dy+1)/(1+dimy)
enddo


mode = 2
do j = 1, ncy+1
   do i = 1, ncx+1
      rho(i,j) = 2_f64 * mode**2 * sin(mode*x(i))*sin(mode*y(j))
   end do
end do

call initialize(poisson, x, y, ncx, ncy)
call solve(poisson, ex, ey, rho, ncx, ncy)

do j = 1, ncy+1
   do i = 1, ncx+1
      write(11,*) x(i), y(j), rho(i,j),  sin(mode*x(i))*sin(mode*y(j))
   end do
   write(11,*)
end do
close(11)

end subroutine test


end program test_poisson_2d_fem
