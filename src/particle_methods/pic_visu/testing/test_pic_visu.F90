program test_pic_visu
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use biot_savart, only: &
      deplace, &
      getrealtimer, &
      initialize, &
      vitesse

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_pic_visu, only: &
      sll_o_distribution_gnuplot, &
      sll_s_distribution_m4_gnuplot, &
      sll_s_distribution_tsc_gnuplot, &
      sll_s_distribution_xdmf, &
      sll_o_particles_center_gnuplot, &
      sll_s_particles_center_gnuplot_inline, &
      sll_o_plot_format_points3d, &
      sll_s_plot_format_xmdv

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   sll_int32  :: nbpart

   nbpart = 100

   call plot_test_2d()

   call test_animation_2d()

   print *, "PASSED"

contains

   subroutine plot_test_2d()
      sll_real64, allocatable, dimension(:, :) :: density
      sll_int32  :: iplot
      sll_real64 :: xmin, xmax, vmin, vmax
      sll_real64 :: time
      sll_int32  :: i
      sll_int32  :: error
      sll_int32  :: nx, nv
      sll_real64 :: t, angle, r
      sll_real64, dimension(:), allocatable :: x
      sll_real64, dimension(:), allocatable :: v
      sll_real64, dimension(:), allocatable :: w

      SLL_ALLOCATE(x(nbpart), error)
      SLL_ALLOCATE(v(nbpart), error)
      SLL_ALLOCATE(w(nbpart), error)

      do i = 1, nbpart
         t = real(i, f64)/real(nbpart - 1, f64)
         angle = t*(sll_p_pi*2.)*50._f64
         r = t*2.
         x(i) = r*cos(angle)
         v(i) = r*sin(angle)
      end do

      w = sqrt(x*x + v*v)

      xmin = -4.0_f64; xmax = 4.0_f64
      vmin = -4.0_f64; vmax = 4.0_f64
      nx = 64
      nv = 64
      SLL_ALLOCATE(density(nx, nv), error)

      iplot = 1
      time = 0._f64
      call sll_o_particles_center_gnuplot("pic_xv", x, v, xmin, xmax, vmin, vmax, iplot, time)
      call sll_o_distribution_gnuplot("df_xv", x, v, xmin, xmax, nx, vmin, vmax, nv, iplot, time)
      call sll_s_particles_center_gnuplot_inline(x, v, xmin, xmax, vmin, vmax, time)
      call sll_o_plot_format_points3d("pic_xv", x, v, iplot)
      call sll_s_plot_format_xmdv("pic_xv", x, v, iplot, xmin, xmax, vmin, vmax)
      call sll_s_distribution_xdmf("df_xv", x, v, w, xmin, xmax, nx, vmin, vmax, nv, iplot)

      call sll_s_distribution_tsc_gnuplot('df_tsc', x, v, w, &
                                          xmin, xmax, nx, &
                                          vmin, vmax, nv, iplot)

      call sll_s_distribution_m4_gnuplot('df_m4', x, v, w, &
                                         xmin, xmax, nx, &
                                         vmin, vmax, nv, iplot)
   end subroutine plot_test_2d

   subroutine test_animation_2d()
      integer :: iplot, istep, imov, nstep = 10
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
      call initialize(imov, xp, yp, op, delta, nbpart)
      SLL_ALLOCATE(up(nbpart), error)
      SLL_ALLOCATE(vp(nbpart), error)

      iplot = 0
      time = 0.0_f64

      xmin = -3.0_f64; xmax = 3.0_f64
      ymin = -2.0_f64; ymax = 2.0_f64

      dt = 0.02_f64

      do istep = 1, nstep       !loop over time

         call vitesse(nbpart, xp, yp, op, up, vp, delta, time)

         !call centres(nbpart, xp, yp, op, time)

         call deplace(nbpart, xp, yp, up, vp, dt)

         time = time + dt
         if (mod(istep, 10) == 0) then
            iplot = iplot + 1
            call sll_o_particles_center_gnuplot("pic_xy", xp, yp, &
                                                xmin, xmax, ymin, ymax, iplot, time)
            call sll_o_plot_format_points3d("pic_xy", xp, yp, op, iplot)
            call sll_s_plot_format_xmdv("pic_xy", xp, yp, iplot, xmin, xmax, ymin, ymax)
         end if

      end do      !next time step

      call cpu_time(tcpu)
      t1 = getRealTimer()
      write (*, "(5x,' CPU time = ', G15.3)") tcpu

   end subroutine test_animation_2d

end program test_pic_visu
