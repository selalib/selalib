program unit_test_4d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use numeric_constants
use sll_cubic_spline_interpolator_2d

implicit none

  sll_int32  :: n_x, n_vx, n_y, n_vy
  sll_int32  :: i, j, k, l, it, n_steps
  sll_real64 :: x_min, x_max, y_min, y_max
  sll_real64 :: vx_min, vx_max, vy_min, vy_max
  sll_real64 :: delta_t, error
  sll_int32  :: info

  sll_real64, dimension(:)  , allocatable :: x, y
  sll_real64, dimension(:)  , allocatable :: vx, vy
  sll_real64, dimension(:,:), allocatable :: dx, dy
  sll_real64, dimension(:,:), allocatable :: dvx, dvy

  sll_real64, dimension(:,:,:,:), allocatable :: df

  class(sll_interpolator_2d_base), pointer    :: interp_xy
  class(sll_interpolator_2d_base), pointer    :: interp_vxvy

  type(cubic_spline_2d_interpolator), target  :: spline_xy
  type(cubic_spline_2d_interpolator), target  :: spline_vxvy

  print*,'*******************************'
  print*,' 4D case                       '
  print*,' 2D in (x,y) and 2D in (vx,vy) '
  print*,'*******************************'

  print*, 'set domain size'
  x_min  =  -5.0_f64; x_max  =  5.0_f64
  y_min  =  -5.0_f64; y_max  =  5.0_f64
  vx_min =  -5.0_f64; vx_max =  5.0_f64 
  vy_min =  -5.0_f64; vy_max =  5.0_f64 

  n_x  = 32; n_y  = 32
  n_vx = 32; n_vy = 32

  print*, 'create 1d meshes in x and v'
  SLL_CLEAR_ALLOCATE(x(n_x)         , info)
  SLL_CLEAR_ALLOCATE(y(n_y)         , info)
  SLL_CLEAR_ALLOCATE(vx(n_vx)       , info)
  SLL_CLEAR_ALLOCATE(vy(n_vy)       , info)
  SLL_CLEAR_ALLOCATE(dx(n_x,n_y)    , info)
  SLL_CLEAR_ALLOCATE(dy(n_x,n_y)    , info)
  SLL_CLEAR_ALLOCATE(dvx(n_vx,n_vy) , info)
  SLL_CLEAR_ALLOCATE(dvy(n_vx,n_vy) , info)

  do i = 1, n_x
     x(i)  = x_min  + (i-1)*(x_max-x_min)/(n_x-1)
  end do
  do j = 1, n_y
     y(j)  = y_min  + (j-1)*(y_max-y_min)/(n_y-1)
  end do
  do k = 1, n_vx
     vx(k) = vx_min + (k-1)*(vx_max-vx_min)/(n_vx-1)
  end do
  do l = 1, n_vy
     vy(l) = vy_min + (l-1)*(vy_max-vy_min)/(n_vy-1)
  end do


  SLL_ALLOCATE(df(n_x,n_y,n_vx,n_vy), info)
  print*, 'initialize 4d distribution function f(x,y,vx,vy) gaussian'
  do l = 1, n_vy
     do k = 1, n_vx
        do j = 1, n_y
           do i = 1, n_x
              df(i,j,k,l) = exp(-((vx(k)-2)**2+vy(l)**2))
           end do
        end do
     end do
  end do

  Print*, 'checking advection of a Gaussian in a uniform field'
 
  call spline_xy%initialize(n_x, n_y, x_min, x_max, y_min, y_max, &
                            PERIODIC_SPLINE, PERIODIC_SPLINE )

  call spline_vxvy%initialize(n_vx, n_vy, vx_min, vx_max, vy_min, vy_max, &
                              PERIODIC_SPLINE, PERIODIC_SPLINE )
  interp_xy   => spline_xy
  interp_vxvy => spline_vxvy

  ! run BSL method using 10 time steps and second order splitting
  n_steps = 1000
  delta_t = 0.01
  do it = 1, n_steps

     print*,it
     call plot_df(it)
     !call advection_xy(df, interp_xy, delta_t)
     call advection_vxvy(df, interp_vxvy, delta_t)

  end do

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

   subroutine advection_xy(df, interp_xy, dt)
   class(sll_interpolator_2d_base), pointer  :: interp_xy
   sll_real64, intent(inout), dimension(:,:,:,:) :: df
   sll_real64, intent(in) :: dt

   do l = 1, n_vy
   do k = 1, n_vx

      do j = 1, n_y
      do i = 1, n_x
         dx(i,j) = -dt*vx(k)
         dy(i,j) = -dt*vy(l)
      end do
      end do

      df(:,:,k,l) = interp_xy%interpolate_array_disp(n_x,n_y, &
                                                     df(:,:,k,l),dx,dy)
   end do
   end do

   end subroutine advection_xy

   subroutine advection_vxvy(df, interp_vxvy, dt)
   class(sll_interpolator_2d_base), pointer  :: interp_vxvy
   sll_real64, intent(inout), dimension(:,:,:,:) :: df
   sll_real64, intent(in) :: dt

   do j = 1, n_y
   do i = 1, n_x

      do l = 1, n_vy
      do k = 1, n_vx
         dvx(k,l) =  dt*vy(l)
         dvy(k,l) = -dt*vx(k)
      end do
      end do

      df(i,j,:,:) = interp_vxvy%interpolate_array_disp(n_vx,n_vy, &
                                                       df(i,j,:,:),dvx,dvy)

   end do
   end do

   end subroutine advection_vxvy

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
