program bsl_1d_quintic_compact
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use sll_utilities, only: int2string
use sll_constants
use sll_quintic_spline_interpolator_1d
#ifndef STDF95
use sll_module_interpolators_1d_base
#endif

implicit none

sll_int32, parameter  :: COMPACT_SPLINE = 0
sll_int32  :: nc_x, nc_v
sll_int32  :: i, j, it, n_steps
sll_real64 :: x_min, x_max, v_min, v_max
sll_real64 :: delta_t, error
sll_int32  :: info
sll_real64 :: delta_x
sll_real64 :: delta_v

sll_real64 :: x, v, xc, vc, time, r
sll_real64, dimension(:,:), allocatable :: df
sll_real64, dimension(:,:), allocatable :: advfield_x
sll_real64, dimension(:,:), allocatable :: advfield_v

#ifdef STDF95
type(quintic_spline_1d_interpolator), pointer  :: interp_x
type(quintic_spline_1d_interpolator), pointer  :: interp_v
#else
class(sll_interpolator_1d_base), pointer     :: interp_x
class(sll_interpolator_1d_base), pointer     :: interp_v
#endif

type(quintic_spline_1d_interpolator), target   :: spline_x
type(quintic_spline_1d_interpolator), target   :: spline_v

print*,'*********************'
print*,' 1D case             '
print*,' 1D in x and 1D in v '
print*,'*********************'

print*, 'set domain size'
x_min =  -5.0_f64; x_max =  5.0_f64
v_min =  -5.0_f64; v_max =  5.0_f64 

r  = 2.0_f64
xc = r
vc = 0.0_f64

nc_x = 100; nc_v = 100

delta_x = (x_max-x_min)/nc_x
delta_v = (v_max-v_min)/nc_v

SLL_ALLOCATE(df(nc_x+1,nc_v+1), info)
SLL_ALLOCATE(advfield_x(nc_x+1,nc_v+1), info)
SLL_ALLOCATE(advfield_v(nc_x+1,nc_v+1), info)

do j = 1, nc_v+1
   do i = 1, nc_x+1
      x = x_min + (i-1)*delta_x 
      v = v_min + (j-1)*delta_v
      df(i,j) = exp(-((x-xc)*(x-xc)+(v-vc)*(v-vc)))
      advfield_x(i,j) = -v
      advfield_v(i,j) =  x
   end do
end do

print*, 'initialize 2d distribution function f(x,v) gaussian'
Print*, 'checking advection of a Gaussian in a uniform field'
#ifdef STDF95
call quintic_spline_initialize(spline_x, nc_x+1, x_min, x_max, COMPACT_SPLINE ) 
call quintic_spline_initialize(spline_v, nc_v+1, v_min, v_max, COMPACT_SPLINE )
#else  
call spline_x%initialize(nc_x+1, x_min, x_max, COMPACT_SPLINE )
call spline_v%initialize(nc_v+1, v_min, v_max, COMPACT_SPLINE )
#endif

interp_x => spline_x
interp_v => spline_v

! run BSL method using 10 time steps
n_steps = 200
delta_t = 0.05
time    = 0.0_f64
call advection_x(0.5*delta_t)
time    = time + 0.5*delta_t
do it = 1, n_steps

   call advection_v(delta_t)
   time = time + 0.5*delta_t
   call plot_df( it )
   call advection_x(delta_t)
   time = time + 0.5*delta_t

   error = 0.0
   xc = r*cos(time)
   vc = r*sin(time)

end do

do j = 1, nc_v+1
   do i = 1, nc_x+1
      x = x_min + (i-1)*delta_x - xc
      v = v_min + (j-1)*delta_v - vc
      error = max(error,abs(df(i,j)-exp(-(x*x+v*v))))
   end do
end do

print*, ' 100x100 nodes, ', it, ' time steps. Error= ', error
print *, 'Successful, exiting program.'
print *, 'PASSED'

contains

   subroutine advection_x(dt)
   sll_real64, intent(in) :: dt
   sll_real64 :: eta

   do j = 1, nc_v
#ifdef STDF95
     call quintic_spline_compute_interpolants(interp_x, df(:,j) )
#else
     call interp_x%compute_interpolants( df(:,j) )
#endif
     do i = 1, nc_x
        eta = x_min + (i-1)*delta_x - dt*advfield_x(i,j)
        eta = max(eta, x_min)
        eta = min(eta, x_max)
#ifdef STDF95
        df(i,j) = interpolate_value_qs1d(interp_x,eta)
#else
        df(i,j) = interp_x%interpolate_value(eta)
#endif
     end do
   end do

   end subroutine advection_x

   subroutine advection_v(dt)
   sll_real64, intent(in) :: dt
   sll_real64 :: eta

   do i = 1, nc_x
#ifdef STDF95
     call quintic_spline_compute_interpolants(interp_v, df(i,:) )
#else
     call interp_v%compute_interpolants( df(i,:) )
#endif
      do j = 1, nc_v
        eta = v_min + (j-1)*delta_v - dt*advfield_v(i,j)
        eta = max(eta, v_min)
        eta = min(eta, v_max)
#ifdef STDF95
        df(i,j) = interpolate_value_qs1d(interp_v,eta)
#else
        df(i,j) = interp_v%interpolate_value(eta)
#endif
     end do
   end do

   end subroutine advection_v

   subroutine plot_df(iplot)

   integer :: iplot, i, j
   character(len=4) :: cplot
 
   call int2string(iplot,cplot)

   open(11, file="df-"//cplot//".dat")
   do j = 1, size(df,2)
      do i = 1, size(df,1)
         x = x_min + (i-1)*(x_max-x_min)/(nc_x)
         v = v_min + (j-1)*(v_max-v_min)/(nc_v)
         write(11,*) sngl(x),sngl(v),sngl(df(i,j))
      end do
      write(11,*)
   end do
   close(11)
   
   open( 90, file = 'df.gnu', position="append" )
   if ( iplot == 1 ) then
      rewind(90)
      !write(90,*)"set cbrange[-1:1]"
      !write(90,*)"set pm3d"
      write(90,*)"set surf"
      write(90,*)"set term x11"
   end if

   write(90,*)"set title 'step = ",iplot,"'"
   write(90,"(a)")"splot 'df-"//cplot//".dat' u 1:2:3 w lines"
   close(90)

   end subroutine plot_df

end program bsl_1d_quintic_compact
