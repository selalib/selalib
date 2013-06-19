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

do j = 1, nc_y
   do i = 1, nc_x
      rho(i,j) = 2_f64 * mode**2 * sin(mode*x(i))*sin(mode*y(j))
   end do
end do

call initialize(poisson, x, y, nc_x+1, nc_y+1)
call solve(poisson, ex, ey, rho)

errmax = 0.
do j = 1, nc_y+1
   do i = 1, nc_x+1
      errmax = max(errmax,abs(rho(i,j)-sin(mode*x(i))*sin(mode*y(j))))
      write(12,*) x(i), y(j), errmax
   end do
   write(12,*) 
end do
close(12)
print*, 'compact, error = ', errmax

do j = 1, nc_y+1
   do i = 1, nc_x+1
      write(11,*) x(i), y(j), rho(i,j), sin(mode*x(i))*sin(mode*y(j))
   end do
   write(11,*)
end do
close(11)


end subroutine test_compact


subroutine test_periodic()
type( poisson_2d_periodic_fem ) :: poisson

do j = 1, nc_y
   do i = 1, nc_x
      rho(i,j) = 2_f64 * mode**2 * cos(mode*x(i))*cos(mode*y(j))
   end do
end do

call initialize(poisson, x, y, nc_x+1, nc_y+1)
call solve(poisson, ex, ey, rho)

errmax = 0.
do j = 2, nc_y
   do i = 2, nc_x
      errmax = max(errmax,abs(rho(i,j)-cos(mode*x(i))*cos(mode*y(j))))
      write(13,*) x(i), y(j), errmax
   end do
   write(13,*) 
end do
close(13)
print*, 'periodic, error = ', errmax

do j = 1, nc_y+1
   do i = 1, nc_x+1
      write(14,*) x(i), y(j), rho(i,j), cos(mode*x(i))*cos(mode*y(j))
   end do
   write(14,*)
end do
close(14)

end subroutine test_periodic


end program test_poisson_2d_fem

subroutine write_mtv_file( x, y )
real(8), dimension(:) :: x
real(8), dimension(:) :: y
integer :: iel, isom, nx, ny
real(8) :: x1, y1

nx = size(x)-1
ny = size(y)-1
open(10, file="mesh.mtv")
write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Elements number ' "
   
do i=1,nx
   do j=1,ny
      write(10,*) x(i  ), y(j  ), 0.
      write(10,*) x(i+1), y(j  ), 0.
      write(10,*) x(i+1), y(j+1), 0.
      write(10,*) x(i  ), y(j+1), 0.
      write(10,*) x(i  ), y(j  ), 0.
      write(10,*)
   end do
end do

!Numeros des elements
do i=2,nx-1
   do j=2,ny-1
      iel = i-1+(j-2)*(nx-1)
      x1 = 0.5*(x(i)+x(i+1))
      y1 = 0.5*(y(j)+y(j+1))
      write(10,"(a)"   ,  advance="no")"@text x1="
      write(10,"(g15.3)", advance="no") x1
      write(10,"(a)"   ,  advance="no")" y1="
      write(10,"(g15.3)", advance="no") y1
      write(10,"(a)"   ,  advance="no")" z1=0. lc=4 ll='"
      write(10,"(i4)"  ,  advance="no") iel
      write(10,"(a)")"'"
   end do
end do

!Numeros des noeud
do i=2,nx
   do j=2,ny
      isom = i-1+(j-2)*(nx-1)
      write(10,"(a)"   ,  advance="no")"@text x1="
      write(10,"(g15.3)", advance="no") x(i)
      write(10,"(a)"   ,  advance="no")" y1="
      write(10,"(g15.3)", advance="no") y(j)
      write(10,"(a)"   ,  advance="no")" z1=0. lc=5 ll='"
      write(10,"(i4)"  ,  advance="no") isom
      write(10,"(a)")"'"
   end do
end do
   
write(10,*)"$END"
   
write(10,*)"$END"
close(10)

end subroutine write_mtv_file

