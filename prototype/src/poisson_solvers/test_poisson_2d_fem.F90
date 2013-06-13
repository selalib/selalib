program test_poisson_2d_fem
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_constants.h"

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
sll_real64 :: errmax

call test()

contains

subroutine test()
type( poisson_fem ) :: poisson
sll_real64 :: dx, dy
sll_int32  :: nx, ny

nx = 64
ny = 64

SLL_CLEAR_ALLOCATE(ex(1:nx,1:ny),error)  
SLL_CLEAR_ALLOCATE(ey(1:nx,1:ny),error) 
SLL_CLEAR_ALLOCATE(rho(1:nx,1:ny),error)  

SLL_ALLOCATE(x(-1:nx+1),error)  
SLL_ALLOCATE(y(-1:ny+1),error) 

dimx = 2 * sll_pi
dimy = 2 * sll_pi

dx = dimx / nx
dy = dimy / ny

x(0) = 0.
y(0) = 0.

do i=1,nx
   x(i) = (i*dx) *(i*dx+1)/(1+dimx)
enddo
do j=1,ny
   y(j) = (j*dy) *(j*dy+1)/(1+dimy)
enddo


mode = 2
do j = 1, ny
   do i = 1, nx
      rho(i,j) = 2_f64 * mode**2 * sin(mode*x(i))*sin(mode*y(j))
   end do
end do

call initialize(poisson, x, y, nx, ny)
call solve(poisson, ex, ey, rho, nx, ny)

errmax = 0.
do j = 1, ny
   do i = 1, nx
      write(11,*) x(i), y(j), rho(i,j),  sin(mode*x(i))*sin(mode*y(j))
      errmax = max(errmax,abs(rho(i,j)-sin(mode*x(i))*sin(mode*y(j))))
   end do
   write(11,*)
end do
close(11)
print*, 'error = ', errmax

end subroutine test


end program test_poisson_2d_fem
