program unit_test_4d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use numeric_constants
use sll_cubic_spline_interpolator_2d

implicit none

  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  real(8) :: OMP_GET_WTIME

  sll_int32  :: n_x, n_vx, n_y, n_vy
  sll_int32  :: i, j, k, l, it, n_steps
  sll_real64 :: x_min, x_max, y_min, y_max
  sll_real64 :: vx_min, vx_max, vy_min, vy_max
  sll_real64 :: dt, error
  sll_int32  :: info
  sll_real64 :: eta1, eta2
  sll_real64 :: delta_x, delta_y, delta_vx, delta_vy
  sll_real64 :: dx, dy
  sll_real64 :: dvx, dvy
  sll_real64 :: t0,t1,t2

  sll_real64, dimension(:)  , allocatable :: x, y
  sll_real64, dimension(:)  , allocatable :: vx, vy

  sll_real64, dimension(:,:,:,:), allocatable :: df

  class(sll_interpolator_2d_base), pointer    :: interp_xy
  class(sll_interpolator_2d_base), pointer    :: interp_vxvy

  type(cubic_spline_2d_interpolator), target  :: spline_xy
  type(cubic_spline_2d_interpolator), target  :: spline_vxvy

  print*,'*******************************'
  print*,' 4D case                       '
  print*,' 2D in (x,y) and 2D in (vx,vy) '
  print*,'*******************************'

  n_x  = 32; n_y  = 32
  n_vx = 32; n_vy = 32

  SLL_CLEAR_ALLOCATE(x(n_x)         , info)
  SLL_CLEAR_ALLOCATE(y(n_y)         , info)
  SLL_CLEAR_ALLOCATE(vx(n_vx)       , info)
  SLL_CLEAR_ALLOCATE(vy(n_vy)       , info)
  SLL_CLEAR_ALLOCATE(df(n_x,n_y,n_vx,n_vy), info)

  !$OMP PARALLEL NUM_THREADS(1) &
  !$OMP DEFAULT(SHARED)           &
  !$OMP PRIVATE(spline_xy, spline_vxvy, interp_xy, interp_vxvy) 

  print*, 'set domain size'
  x_min  =  -5.0_f64; x_max  =  5.0_f64
  y_min  =  -5.0_f64; y_max  =  5.0_f64
  vx_min =  -5.0_f64; vx_max =  5.0_f64 
  vy_min =  -5.0_f64; vy_max =  5.0_f64 
  
  delta_x  = (x_max-x_min)/(n_x-1)
  delta_y  = (y_max-y_min)/(n_y-1)
  delta_vx = (vx_max-vx_min)/(n_vx-1)
  delta_vy = (vy_max-vy_min)/(n_vy-1)

  call spline_xy%initialize(n_x, n_y, x_min, x_max, y_min, y_max, &
                            PERIODIC_SPLINE, PERIODIC_SPLINE )
  interp_xy   => spline_xy

  call spline_vxvy%initialize(n_vx, n_vy, vx_min, vx_max, vy_min, vy_max, &
                              PERIODIC_SPLINE, PERIODIC_SPLINE )
  interp_vxvy => spline_vxvy

  do i = 1, n_x
     x(i)  = x_min  + (i-1)*delta_x
  end do
  do j = 1, n_y
     y(j)  = y_min  + (j-1)*delta_y
  end do
  do k = 1, n_vx
     vx(k) = vx_min + (k-1)*delta_vx
  end do
  do l = 1, n_vy
     vy(l) = vy_min + (l-1)*delta_vy
  end do

  !$OMP BARRIER
  !$OMP SINGLE
  do l = 1, n_vy
     do k = 1, n_vx
        do j = 1, n_y
           do i = 1, n_x
              df(i,j,k,l) = exp(-((vx(k)-2)**2+vy(l)**2))
           end do
        end do
     end do
  end do
 
  !$OMP END SINGLE

  ! run BSL method using 10 time steps and second order splitting
  n_steps = 1000
  dt = 0.01

  do it = 1, n_steps

     !$OMP SINGLE
     call plot_df(it)
     !$OMP END SINGLE


#ifndef DEBUG
     !$OMP BARRIER
     t0=OMP_GET_WTIME()
#endif

     !$OMP DO 
     do l = 1, n_vy
     do k = 1, n_vx

       call interp_xy%compute_interpolants(df(:,:,k,l))

       do j = 1, n_y
        do i = 1, n_x
           dx = dt*vx(k)
           dy = dt*vy(l)
           eta1 = x_min + (i-1)*delta_x
           eta2 = y_min + (j-1)*delta_y
           eta1 = x_min + modulo(eta1-x_min-dx,x_max-x_min)
           eta2 = y_min + modulo(eta2-y_min-dy,y_max-y_min)
           df(i,j,k,l) = interp_xy%interpolate_value(eta1,eta2)
        end do
        end do
  
     end do
     end do

#ifndef DEBUG
     !$OMP BARRIER
     t1=OMP_GET_WTIME()
#endif

     do j = 1, n_y
     do i = 1, n_x

        call interp_vxvy%compute_interpolants(df(i,j,:,:))

        do l = 1, n_vy
        do k = 1, n_vx
           dvx = -dt*vy(l)
           dvy =  dt*vx(k)
           eta1 = vx_min + (k-1)*delta_vx
           eta2 = vy_min + (l-1)*delta_vy
           eta1 = vx_min + modulo(eta1-vx_min-dvx,vx_max-vx_min)
           eta2 = vy_min + modulo(eta2-vy_min-dvy,vy_max-vy_min)
           df(i,j,k,l) = interp_vxvy%interpolate_value(eta1,eta2)
        end do
        end do

     end do
     end do

#ifndef DEBUG
     !$OMP BARRIER
     t2=OMP_GET_WTIME()

     !$OMP CRITICAL
     print *, it,OMP_GET_NUM_THREADS(),OMP_GET_THREAD_NUM(),t1-t0,t2-t1
     !$OMP END CRITICAL
#endif

  end do

!$OMP END PARALLEL

  do l=1, n_vy
     do k=1, n_vx
        do j=1, n_y
           do i=1, n_y
              error = max(error,abs(df(i,j,k,l)-exp(-(vx(k)**2+vy(l)**2))))
           end do
        end do
     end do
  end do

  print*, ' 100 nodes, ', it, ' time steps. Error= ', error


  print *, 'Successful, exiting program.'
  print *, 'PASSED'


contains

   subroutine plot_df(iplot)

   integer :: iplot
   character(len=4) :: cplot
 
   call int2string(iplot,cplot)

   open(11, file="df-"//cplot//".dat")
   do l = 1, size(df,4)
      do k = 1, size(df,3)
         write(11,*) sngl(vx(k)),sngl(vy(l)),sngl(sum(df(:,:,k,l)))
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

end program unit_test_4d
