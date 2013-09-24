program test_poisson_2d_fem
#include "sll_poisson_solvers_macros.h"
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_constants.h"

use sll_poisson_2d_fem
use sll_poisson_2d_periodic_fem
implicit none

sll_int32  :: i, j
sll_real64 :: dimx, dimy
sll_real64 :: dx, dy
sll_int32  :: nc_x, nc_y
sll_int32  :: error
sll_real64 :: mode
sll_real64, dimension(:),   pointer :: x
sll_real64, dimension(:),   pointer :: y
sll_real64, dimension(:,:), pointer :: ex
sll_real64, dimension(:,:), pointer :: ey
sll_real64, dimension(:,:), pointer :: rho
sll_real64 :: errmax

nc_x = 32
nc_y = 32

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
   x(i) = (i-1)*dx !* 0.5 * ((i-1)*dx+1)
enddo
do j=1,nc_y+1
   y(j) = (j-1)*dy !* 0.5 * ((j-1)*dy+1)
enddo

mode = 4*sll_pi

call test_compact()
call test_periodic()

contains

subroutine test_compact()
type( poisson_fem ) :: poisson

do j = 1, nc_y+1
   do i = 1, nc_x+1
      rho(i,j) = 2_f64 * mode**2 * sin(mode*x(i))*sin(mode*y(j))
   end do
end do

call initialize(poisson, x, y, nc_x+1, nc_y+1)
call solve(poisson, ex, ey, rho)

errmax = 0.
do j = 1, nc_y+1
   do i = 1, nc_x+1
      errmax = errmax +abs(rho(i,j)-sin(mode*x(i))*sin(mode*y(j)))
      write(11,*) x(i), y(j), rho(i,j), sin(mode*x(i))*sin(mode*y(j))
      write(12,*) x(i), y(j), rho(i,j)-sin(mode*x(i))*sin(mode*y(j))
   end do
   write(12,*); write(11,*)
end do
close(11); close(12)
print*, 'compact, error = ', errmax / (nc_x+1) / (nc_y+1)

end subroutine test_compact

subroutine test_periodic()
type( poisson_2d_periodic_fem ) :: poisson

do j = 1, nc_y+1
   do i = 1, nc_x+1
      rho(i,j) = 2_f64 * mode**2 * sin(mode*x(i))*sin(mode*y(j))
   end do
end do

call initialize(poisson, x, y, nc_x+1, nc_y+1)
call solve(poisson, ex, ey, rho)

errmax = 0.
do j = 1, nc_y
   do i = 1, nc_x
      errmax = errmax+abs(rho(i,j)-sin(mode*x(i))*sin(mode*y(j)))
      write(13,*) x(i), y(j), rho(i,j),sin(mode*x(i))*sin(mode*y(j))
      write(14,*) x(i), y(j), rho(i,j)-sin(mode*x(i))*sin(mode*y(j))
   end do
   write(13,*) ; write(14,*) 
end do
close(13); close(14)
print*, 'periodic, error = ', errmax / (nc_x*nc_y)

end subroutine test_periodic

end program test_poisson_2d_fem

