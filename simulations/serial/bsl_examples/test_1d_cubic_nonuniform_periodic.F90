program bsl_1d_cubic_nonuniform_periodic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic

   use sll_m_cubic_spline_interpolator_1d_nonuniform, only: &
      sll_t_cubic_spline_interpolator_1d_nonuniform

   use sll_m_interpolators_1d_base, only: &
      sll_c_interpolator_1d

   use sll_m_utilities, only: &
      sll_s_int2string

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   sll_int32  :: nc_x, nc_v
   sll_int32  :: i, j, it, n_steps
   sll_real64 :: x_min, x_max, v_min, v_max
   sll_real64 :: delta_t, error
   sll_int32  :: info

   sll_real64, dimension(:), allocatable :: x
   sll_real64, dimension(:), allocatable :: v
   sll_real64, dimension(:, :), allocatable :: df

   sll_real64 :: advfield_x, advfield_v

   class(sll_c_interpolator_1d), pointer     :: interp_x
   class(sll_c_interpolator_1d), pointer     :: interp_v

   type(sll_t_cubic_spline_interpolator_1d_nonuniform), target   :: spline_x
   type(sll_t_cubic_spline_interpolator_1d_nonuniform), target   :: spline_v

   print *, '***********************'
   print *, '* 1D case             *'
   print *, '* 1D in x and 1D in v *'
   print *, '***********************'

   print *, 'set domain size'
   x_min = -5.0_f64; x_max = 5.0_f64
   v_min = -5.0_f64; v_max = 5.0_f64

   nc_x = 100; nc_v = 100
   SLL_ALLOCATE(df(nc_x + 1, nc_v + 1), info)

   do j = 1, nc_v + 1
      do i = 1, nc_x + 1
         x(i) = x_min + (i - 1)*(i - 1)*(x_max - x_min)/nc_x/nc_x
         v(j) = v_min + (j - 1)*(j - 1)*(v_max - v_min)/nc_v/nc_v
         df(i, j) = exp(-(x(i)*x(i) + v(j)*v(j)))
      end do
   end do

   advfield_x = 1.0_f64
   advfield_v = 0.0_f64

   print *, 'initialize 2d distribution function f(x,v) sll_m_gaussian'
   print *, 'checking advection of a Gaussian in a uniform field'
   call spline_x%init(nc_x + 1, x_min, x_max, sll_p_periodic)
   call spline_v%init(nc_v + 1, v_min, v_max, sll_p_periodic)

   interp_x => spline_x
   interp_v => spline_v

! run BSL method using 10 time steps
   n_steps = 100
   delta_t = 10.0_f64/n_steps
   do it = 1, n_steps

      call plot_df(it)
      call advection_x(0.5*delta_t)
      call advection_v(delta_t)
      call advection_x(0.5*delta_t)

   end do

! compute error when Gaussian arrives at center (t=1)
   error = 0.0_f64
   do j = 1, nc_v + 1
      do i = 1, nc_x + 1
         error = max(error, abs(df(i, j) - exp(-(x(i)*x(i) + v(j)*v(j)))))
      end do
   end do

   print *, ' 100 nodes, ', it, ' time steps. Error= ', error
   print *, 'Successful, exiting program.'
   print *, 'PASSED'

contains

   subroutine advection_x(dt)
      sll_real64, intent(in) :: dt

      do j = 1, nc_v
         call interp_x%interpolate_array_disp_inplace(nc_x + 1, df(:, j), -dt*advfield_x)
      end do

   end subroutine advection_x

   subroutine advection_v(dt)
      sll_real64, intent(in) :: dt

      do i = 1, nc_x
         call interp_v%interpolate_array_disp_inplace(nc_v + 1, df(i, :), -dt*advfield_v)
      end do

   end subroutine advection_v

   subroutine plot_df(iplot)

      integer :: iplot, i, j
      character(len=4) :: cplot

      call sll_s_int2string(iplot, cplot)

      open (11, file="df-"//cplot//".dat")
      do j = 1, size(df, 2)
         do i = 1, size(df, 1)
            x = x_min + (i - 1)*(x_max - x_min)/(nc_x)
            v = v_min + (j - 1)*(v_max - v_min)/(nc_v)
            write (11, *) sngl(x), sngl(v), sngl(df(i, j))
         end do
         write (11, *)
      end do
      close (11)

      open (90, file='df.gnu', position="append")
      if (iplot == 1) then
         rewind (90)
         !write(90,*)"set cbrange[-1:1]"
         !write(90,*)"set pm3d"
         write (90, *) "set surf"
         write (90, *) "set term x11"
      end if

      write (90, *) "set title 'step = ", iplot, "'"
      write (90, "(a)") "splot 'df-"//cplot//".dat' u 1:2:3 w lines"
      close (90)

   end subroutine plot_df

end program bsl_1d_cubic_nonuniform_periodic
