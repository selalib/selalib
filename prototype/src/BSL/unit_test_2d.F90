program unit_test_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use numeric_constants
use sll_module_interpolators_2d_base
use sll_cubic_spline_interpolator_2d

implicit none

  sll_int32  :: n_x, n_vx, n_y, n_vy
  sll_int32  :: i, j, k, l, it, n_steps
  sll_real64 :: x_min, x_max, y_min, y_max
  sll_real64 :: vx_min, vx_max, vy_min, vy_max
  sll_real64 :: delta_t, error
  sll_int32  :: info

  sll_real64, dimension(:)      , allocatable :: x, y
  sll_real64, dimension(:)      , allocatable :: vx, vy
  sll_real64, dimension(:,:,:,:), allocatable :: df

  sll_real64, dimension(:), allocatable :: f_x, f_y

  class(sll_interpolator_2d_base), pointer    :: interp_xy
  class(sll_interpolator_2d_base), pointer    :: interp_vxvy

  type(cubic_spline_2d_interpolator), target  :: spline_xy
  type(cubic_spline_2d_interpolator), target  :: spline_vxvy

  print*,'*******************************'
  print*,' 2D case                       '
  print*,' 2D in (x,y) and 2D in (vx,vy) '
  print*,'*******************************'

  print*, 'set domain size'
  x_min =  -5.0_f64; x_max =  5.0_f64
  y_min =  -5.0_f64; y_max =  5.0_f64
  vx_min =  -5.0_f64; vx_max =  5.0_f64 
  vy_min =  -5.0_f64; vy_max =  5.0_f64 

  n_x  = 40; n_y  = 40
  n_vx = 40; n_vy = 40

  print*, 'create 1d meshes in x and v'
  SLL_ALLOCATE(x(n_x),   info)
  SLL_ALLOCATE(y(n_y),   info)
  SLL_ALLOCATE(f_x(n_x), info)
  SLL_ALLOCATE(f_y(n_y), info)
  SLL_ALLOCATE(vx(n_vx), info)
  SLL_ALLOCATE(vy(n_vy), info)

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
              df(i,j,k,l) = exp(-(vx(k)**2+vy(l)**2))
           end do
        end do
     end do
  end do

  f_x = 1.0_f64 
  f_y = 0.0_f64 

  Print*, 'checking advection of a Gaussian in a uniform field'
  
  call spline_xy%initialize(n_x, n_y, x_min, x_max, y_min, y_max, &
                            PERIODIC_SPLINE, PERIODIC_SPLINE )

  call spline_vxvy%initialize(n_vx, n_vy, vx_min, vx_max, vy_min, vy_max, &
                              PERIODIC_SPLINE, PERIODIC_SPLINE )
  interp_xy   => spline_xy
  interp_vxvy => spline_vxvy

  ! run BSL method using 10 time steps and second order splitting
  n_steps = 100
  delta_t = 10.0_f64/n_steps
  do it = 1, n_steps

     call advection_xy(df, interp_xy, delta_t)
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
   !type(cubic_spline_2d_interpolator)  :: interp_xy
   class(sll_interpolator_2d_base), pointer  :: interp_xy
   sll_real64, intent(inout), dimension(:,:,:,:) :: df
   sll_real64, intent(in) :: dt
   sll_real64 :: dx, dy

   do l = 1, n_vy
     do k = 1, n_vx
        call interp_xy%compute_interpolants(df(:,:,k,l))
        do j = 1, n_y
           dy = y_min + modulo(y(j)-y_min-dt*vy(l),y_max-y_min)
           do i = 1, n_x
              dx = x_min + modulo(x(i)-x_min-dt*vx(k),x_max-x_min)
              if( dx < x_min .or. dx > x_max) stop 'erreur x'
              if( dy < y_min .or. dy > y_max) stop 'erreur y'
              df(i,j,k,l) = interp_xy%interpolate_value(dx,dy)
           end do
        end do
     end do
   end do

   end subroutine advection_xy

   subroutine advection_vxvy(df, interp_vxvy, dt)
   !type(cubic_spline_2d_interpolator)  :: interp_vxvy
   class(sll_interpolator_2d_base), pointer  :: interp_vxvy
   sll_real64, intent(inout), dimension(:,:,:,:) :: df
   sll_real64, intent(in) :: dt
   sll_real64 :: dvx, dvy

   do j = 1, n_y
      do i = 1, n_x
         call interp_vxvy%compute_interpolants(df(i,j,:,:))
         do k = 1, n_vx
            dvx = vx_min + modulo(vx(k)-vx_min-dt*f_x(i),vx_max-vx_min)
            do l = 1, n_vy
               dvy = vy_min + modulo(vy(l)-vy_min-dt*f_y(j),vy_max-vy_min)
               df(i,j,k,l) = interp_vxvy%interpolate_value(dvx,dvy)
            end do
         end do
      end do
   end do

   end subroutine advection_vxvy


end program unit_test_2d
