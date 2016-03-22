program test_interpolation_m4
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use zone, only: readin, dimx, dimy, dx, dy
use particules
use sll_m_gnuplot

implicit none

type(tm_mesh_fields) :: f
type(particle)       :: p

sll_int32  :: i, j, k, l, error
sll_real64 :: xp, yp
sll_real64 :: xmin, xmax, ymin, ymax

nx  = 100 
ny  = 100 

pi = 4.0_f64 * atan(1.)

xmin = -5.0_f64; xmax = 5.0_f64
ymin = -5.0_f64; ymax = 5.0_f64

dimx = xmax - xmin
dimy = ymax - ymin

dx = dimx / nx
dy = dimy / ny

SLL_CLEAR_ALLOCATE(f%ex(0:nx,0:ny), error)
SLL_CLEAR_ALLOCATE(f%ey(0:nx,0:ny), error)
SLL_CLEAR_ALLOCATE(f%bz(0:nx,0:ny), error)
SLL_CLEAR_ALLOCATE(f%r0(0:nx,0:ny), error) 

print*, 'initialize particles '

nbpart = 10*(nx-10)*(ny-10)
SLL_ALLOCATE(p%dpx(nbpart),error)
SLL_ALLOCATE(p%dpy(nbpart),error)
SLL_ALLOCATE(p%idx(nbpart),error)
SLL_ALLOCATE(p%idy(nbpart),error)
SLL_ALLOCATE(p%epx(nbpart),error)
SLL_ALLOCATE(p%epy(nbpart),error)
SLL_ALLOCATE(p%bpz(nbpart),error)
SLL_ALLOCATE(p%p(nbpart),error)

call random_number(p%dpx)
call random_number(p%dpy)

k = 0
do j = 5, ny-6
do i = 5, nx-6
do l = 1, 10

 k = k+1
 p%idx(k) = i
 p%idy(k) = j
 xp = xmin + (i+p%dpx(k))*dx
 yp = ymin + (j+p%dpy(k))*dy

 p%p(k) = dimx * dimy * exp(-(xp*xp+yp*yp)) / real(nbpart,f64)

 write(10,*) xp, yp, p%p(k)

end do
end do
end do

print*, sum(p%p)

print*, 'compute rho'
call calcul_rho( p, f )
!gnuplot -p rho.gnu (to plot the initial rho)
call sll_o_gnuplot_2d(xmin, xmax, nx+1, &
                      ymin, ymax, ny+1, &
                      f%r0, 'rho', 1, error)

call calcul_rho_m4( p, f )
!gnuplot -p rho.gnu (to plot the initial rho)
call sll_o_gnuplot_2d(xmin, xmax, nx+1, &
                      ymin, ymax, ny+1, &
                      f%r0, 'rho_m4', 1, error)

do j = 0, ny
do i = 0, nx

 xp = xmin + i*dx 
 yp = ymin + j*dy 

 f%ex(i,j) = exp(-(xp*xp+yp*yp)) 
 f%ey(i,j) = exp(-(xp*xp+yp*yp)) 
 f%bz(i,j) = exp(-(xp*xp+yp*yp)) 

 write(11,*) sngl(xp), sngl(yp), sngl(f%ex(i,j)), sngl(f%ey(i,j)), sngl(f%bz(i,j))
end do
write(11,*) 
end do

call interpol_eb( f, p )

do k = 1, nbpart
   write(12,*) sngl(xmin+(p%idx(k)+p%dpx(k))*dx), &
               sngl(ymin+(p%idy(k)+p%dpy(k))*dy), &
               sngl(p%epx(k)),             &
               sngl(p%epy(k)),             &
               sngl(p%bpz(k))
end do


end program test_interpolation_m4
