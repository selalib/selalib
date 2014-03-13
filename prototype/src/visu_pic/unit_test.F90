program test_visu_pic
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use sll_constants
use sll_visu_pic
use sll_utilities
use biot_savart

implicit none

call plot_test_2d()

call test_animation_2d()

print*,"PASSED"

contains

subroutine plot_test_2d()
sll_real64, allocatable, dimension(:,:) :: density
sll_int32  :: iplot
sll_real64 :: xmin, xmax, vmin, vmax
sll_real64 :: time
sll_int32  :: i
sll_int32  :: error
sll_int32  :: nbpart
sll_int32  :: nx, nv
sll_real64 :: t, angle, r
sll_real64, dimension(:), allocatable :: x
sll_real64, dimension(:), allocatable :: v
sll_real64, dimension(:), allocatable :: w

nbpart = 10000
SLL_ALLOCATE(x(nbpart),error)
SLL_ALLOCATE(v(nbpart),error)
SLL_ALLOCATE(w(nbpart),error)

do i = 1, nbpart 
   t = float(i) / float(nbpart-1)
   angle = t * (sll_pi * 2.) * 50.
   r = t * 2.
   x(i) = r * cos(angle)
   v(i) = r * sin(angle)
end do

w = sqrt(x*x+v*v)

xmin = -4; xmax = 4
vmin = -4; vmax = 4
nx = 64
nv = 64
SLL_ALLOCATE(density(nx,nv), error)

iplot = 1
time = 0._f64
call particles_center_gnuplot( "pic_xv", x, v, xmin, xmax, vmin, vmax, iplot, time )
call distribution_gnuplot( "df_xv", x, v, xmin, xmax, nx, vmin, vmax, nv, iplot, time)  
call particles_center_gnuplot_inline( x, v, xmin, xmax, vmin, vmax, time )
call plot_format_points3d( "pic_xv", x, v, iplot)
call plot_format_xmdv( "pic_xv", x, v, iplot, xmin, xmax, vmin, vmax)
call distribution_xdmf("df_xv", x, v, w, xmin, xmax, nx, vmin, vmax, nv, iplot)  

call distribution_tsc_gnuplot('df_tsc', x, v, w, &
                             xmin, xmax, nx,     &
                             vmin, vmax, nv, iplot)  

call distribution_m4_gnuplot('df_m4', x, v, w, &
                             xmin, xmax, nx,     &
                             vmin, vmax, nv, iplot)  
end subroutine plot_test_2d

subroutine test_animation_2d()
integer :: nbpart, iplot, istep, imov, nstep
sll_real64, dimension(:), pointer :: xp
sll_real64, dimension(:), pointer :: yp
sll_real64, dimension(:), pointer :: op
sll_real64, dimension(:), pointer :: up
sll_real64, dimension(:), pointer :: vp
sll_real64 :: time
sll_real64 :: t0, t1, tcpu
sll_real64 :: xmin, xmax, ymin, ymax
sll_real64 :: dt, delta
sll_int32  :: error

call cpu_time(tcpu)
t0 = getRealTimer()
call initialize( nstep, imov, xp, yp, op, delta, dt, nbpart )
SLL_ALLOCATE(up(nbpart), error)
SLL_ALLOCATE(vp(nbpart), error)
iplot   = 0
time = 0.0

xmin = -3; xmax = 3.
ymin = -2; ymax = 2.

print*, nstep
do istep = 1, nstep       !loop over time
   
   call vitesse(nbpart, xp, yp, op, up, vp, delta, time)

   !call centres(nbpart, xp, yp, op, time)

   call deplace(nbpart, xp, yp, up, vp, dt)

   time = time + dt
   if ( mod(istep, 10) == 0) then
      iplot = iplot + 1
      call particles_center_gnuplot( "pic_xy", xp, yp, &
           xmin, xmax, ymin, ymax, iplot, time )
      call plot_format_points3d( "pic_xy", xp, yp, op, iplot)
      call plot_format_xmdv( "pic_xy", xp, yp, iplot, xmin, xmax, ymin, ymax)
   end if

end do      !next time step

call cpu_time(tcpu)
t1 = getRealTimer()
write(*,"(5x,' CPU time = ', G15.3)") tcpu

end subroutine test_animation_2d


end program test_visu_pic
