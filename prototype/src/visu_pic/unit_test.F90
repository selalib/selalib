program test_visu_pic
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use numeric_constants
use sll_visu_pic
use sll_misc_utils

implicit none

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
   r = t * 10.
   x(i) = r * cos(angle)
   v(i) = r * sin(angle)
end do

w = sqrt(x*x+v*v)

call plot_test_2d()

print*,"PASSED"

contains

subroutine plot_test_2d()
sll_real64, allocatable, dimension(:,:) :: density
sll_int32 :: error
!sll_int32  :: file_id
sll_int32  :: iplot
sll_real64 :: xmin, xmax, vmin, vmax
sll_real64 :: time

xmin = -12; xmax = 12
vmin = -12; vmax = 12
nx = 32
nv = 64
SLL_ALLOCATE(density(nx,nv), error)

iplot = 1
time = 0._f64
call particles_center_gnuplot( "pic_xv", x, v, xmin, xmax, vmin, vmax, iplot, time )
call distribution_gnuplot( "df_xv", x, v, xmin, xmax, nx, vmin, vmax, nv, iplot, time)  
call xv_particles_center_gnuplot_inline( x, v, xmin, xmax, vmin, vmax, time )
call plot_format_points3d( "pic_xv", x, v, iplot)
call plot_format_xmdv( "pic_xv", x, v, iplot, xmin, xmax, vmin, vmax)
call distribution_xdmf("df_xv", x, v, w, xmin, xmax, nx, vmin, vmax, nv, iplot)  

end subroutine plot_test_2d

end program test_visu_pic
