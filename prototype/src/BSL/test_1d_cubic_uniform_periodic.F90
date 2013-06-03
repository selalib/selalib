program bsl_1d_cubic_periodic
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_utilities, only: int2string
  use sll_constants
  use sll_cubic_spline_interpolator_1d
#ifndef STDF95
  use sll_module_interpolators_1d_base
#endif

  implicit none

  sll_int32  :: nc_x, nc_v
  sll_int32  :: i, j, it, n_steps
  sll_real64 :: x_min, x_max, v_min, v_max
  sll_real64 :: delta_t, error
  sll_int32  :: info

  sll_real64 :: x, v
  sll_real64, dimension(:,:) , allocatable :: df

  sll_real64 :: advfield_x, advfield_v

#ifdef STDF95
  type(cubic_spline_1d_interpolator), pointer  :: interp_x
  type(cubic_spline_1d_interpolator), pointer  :: interp_v
#else
  class(sll_interpolator_1d_base), pointer    :: interp_x
  class(sll_interpolator_1d_base), pointer    :: interp_v
#endif

  type(cubic_spline_1d_interpolator), target  :: spline_x
  type(cubic_spline_1d_interpolator), target  :: spline_v

  print*,'*********************'
  print*,' 1D case             '
  print*,' 1D in x and 1D in v '
  print*,'*********************'

  print*, 'set domain size'
  x_min =  -5.0_f64; x_max =  5.0_f64
  v_min =  -5.0_f64; v_max =  5.0_f64 

  nc_x = 100; nc_v = 100
  SLL_ALLOCATE(df(nc_x+1,nc_v+1), info)

  do j = 1, nc_v+1
     do i = 1, nc_x+1
        x = x_min + (i-1)*(x_max-x_min)/nc_x
        v = v_min + (j-1)*(v_max-v_min)/nc_v
        df(i,j) =  exp(-(x*x+v*v))
     end do
  end do

  advfield_x = 1_f64 
  advfield_v = 0.0 

  print*, 'initialize 2d distribution function f(x,v) gaussian'

  Print*, 'checking advection of a Gaussian in a uniform field'
#ifdef STDF95
  call cubic_spline_1d_interpolator_initialize(spline_x, nc_x+1, x_min, x_max, SLL_PERIODIC )
  call cubic_spline_1d_interpolator_initialize(spline_v, nc_v+1, v_min, v_max, SLL_PERIODIC )
#else  
  call spline_x%initialize(nc_x+1, x_min, x_max, SLL_PERIODIC )
  call spline_v%initialize(nc_v+1, v_min, v_max, SLL_PERIODIC )
#endif

  interp_x => spline_x
  interp_v => spline_v

  ! run BSL method using 10 time steps
  n_steps = 100
  delta_t = 10.0_f64/n_steps
  do it = 1, n_steps

     call plot_df( it )

     call advection_x(0.5*delta_t)
     call advection_v(    delta_t)
     call advection_x(0.5*delta_t)

  end do

  ! compute error when Gaussian arrives at center (t=1)
  error = 0.0
  do j = 1, nc_v+1
     do i = 1, nc_x+1
        x = x_min + (i-1)*(x_max-x_min)/nc_x
        v = v_min + (j-1)*(v_max-v_min)/nc_v
        error = max(error,abs(df(i,j)-exp(-(x*x+v*v))))
     end do
  end do

  print*, ' 100 nodes, ', it, ' time steps. Error= ', error

  print *, 'Successful, exiting program.'
  print *, 'PASSED'

contains

   subroutine advection_x(dt)
   sll_real64, intent(in) :: dt

     do j = 1, nc_v
#ifdef STDF95
        df(:,j) = cubic_spline_interpolate_array_at_displacement(interp_x,nc_x+1,df(:,j),dt*advfield_x)
#else
        df(:,j) = interp_x%interpolate_array_disp(nc_x+1,df(:,j),dt*advfield_x)
#endif
     end do

   end subroutine advection_x

   subroutine advection_v(dt)
   sll_real64, intent(in) :: dt

     do i = 1, nc_x
#ifdef STDF95
        df(i,:) = cubic_spline_interpolate_array_at_displacement(interp_v,nc_v+1,df(i,:),dt*advfield_v)
#else
        df(i,:) = interp_v%interpolate_array_disp(nc_v+1,df(i,:),dt*advfield_v)
#endif
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

end program bsl_1d_cubic_periodic
