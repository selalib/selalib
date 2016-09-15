program test_interpolation_m4
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use m_zone, only: readin, dimx, dimy, dx, dy
use m_particules
use sll_m_gnuplot
use sll_m_poisson_2d_base
use sll_m_poisson_2d_periodic

implicit none

type(mesh_fields)                     :: f
type(particle)                        :: p
class(sll_c_poisson_2d_base), pointer :: poisson

sll_int32  :: npm
sll_int32  :: i, j, k, l, error
sll_real64 :: xp, yp
sll_real64 :: xmin, xmax, ymin, ymax
sll_real64, allocatable :: x(:,:), y(:,:)

nx  = 40 
ny  = 40 
allocate(x(0:nx,0:ny))
allocate(y(0:nx,0:ny))

npm = 5

pi = 4.0_f64 * atan(1._f64)

xmin = 0.0_f64; xmax = 1.0_f64
ymin = 0.0_f64; ymax = 1.0_f64

dimx = xmax - xmin
dimy = ymax - ymin

dx = dimx / real(nx,f64)
dy = dimy / real(ny,f64)

do j = 0, ny
  do i = 0, nx
    x(i,j) = real(i,f64)*dx    
    y(i,j) = real(j,f64)*dx    
  end do
end do

SLL_CLEAR_ALLOCATE(f%ex(0:nx,0:ny), error)
SLL_CLEAR_ALLOCATE(f%ey(0:nx,0:ny), error)
SLL_CLEAR_ALLOCATE(f%r0(0:nx,0:ny), error) 

print*, 'initialize particles '

nbpart = npm*nx*ny
SLL_ALLOCATE(p%dpx(nbpart),error)
SLL_ALLOCATE(p%dpy(nbpart),error)
SLL_ALLOCATE(p%idx(nbpart),error)
SLL_ALLOCATE(p%idy(nbpart),error)
SLL_ALLOCATE(p%epx(nbpart),error)
SLL_ALLOCATE(p%epy(nbpart),error)
SLL_ALLOCATE(p%p(nbpart),error)

call random_number(p%dpx)
call random_number(p%dpy)

k = 0
do j = 0, ny-1
do i = 0, nx-1
do l = 1, npm

 k = k+1
 p%idx(k) = i
 p%idy(k) = j
 xp = (real(i,f64)+p%dpx(k))*dx 
 yp = (real(j,f64)+p%dpy(k))*dy 

 p%p(k) = (1._f64+cos(2._f64*pi*xp)*sin(2._f64*pi*yp)) / real(nbpart,f64)

 write(10,*) xp, yp, p%p(k)

end do
end do
end do

print*, 'rho total ', sum(p%p)

print*, 'compute rho'
call calcul_rho( p, f )
print*, 'rho error ', sum(abs(f%r0-cos(2*pi*x)*sin(2*pi*y)))/real(nx*ny,f64)

!gnuplot -p rho.gnu (to plot the initial rho)
call sll_o_gnuplot_2d(xmin, xmax, nx+1, &
                      ymin, ymax, ny+1, &
                      f%r0, 'rho', 1, error)

print*, 'compute rho M4'
call calcul_rho_m4( p, f )
print*, 'rho error ', sum(abs(f%r0-cos(2*pi*x)*sin(2*pi*y)))/real(nx*ny,f64)
!gnuplot -p rho_m4.gnu (to plot the initial rho)
call sll_o_gnuplot_2d(xmin, xmax, nx+1, &
                      ymin, ymax, ny+1, &
                      f%r0, 'rho_m4', 1, error)

print*, 'compute rho M6'
call calcul_rho_m6( p, f )
print*, 'rho error ', sum(abs(f%r0-cos(2*pi*x)*sin(2*pi*y)))/real(nx*ny,f64)
!gnuplot -p rho_m6.gnu (to plot the initial rho)
call sll_o_gnuplot_2d(xmin, xmax, nx+1, &
                      ymin, ymax, ny+1, &
                      f%r0, 'rho_m6', 1, error)

poisson => sll_f_new_poisson_2d_periodic(xmin,xmax,nx,ymin,ymax,ny)
call poisson%compute_e_from_rho( f%ex(0:nx,0:ny), &
         f%ey(0:nx,0:ny), f%r0(0:nx,0:ny))

do j = 0, ny
do i = 0, nx

 xp = xmin + real(i,f64)*dx 
 yp = ymin + real(j,f64)*dy

 write(11,*) sngl(xp), sngl(yp), sngl(f%ex(i,j)), sngl(f%ey(i,j))

end do
write(11,*) 
end do

call interpol_eb_m6( f, p )

do j = 0, ny
do i = 0, nx

 xp = xmin + real(i,f64)*dx
 yp = ymin + real(j,f64)*dy

 write(12,*) sngl(xp), sngl(yp), sngl(f%ex(i,j)), sngl(f%ey(i,j))

end do
write(12,*) 
end do

print*, "x momentum error =", abs(sum(p%epx*p%p)- &
        sum(f%ex(0:nx-1,0:ny-1)*f%r0(0:nx-1,0:ny-1))*dx*dy)
print*, "y momentum error =", abs(sum(p%epy*p%p)- &
        sum(f%ey(0:nx-1,0:ny-1)*f%r0(0:nx-1,0:ny-1))*dx*dy)

end program test_interpolation_m4
