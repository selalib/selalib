!**************************************************************
!  Copyright INRIA
!  Authors :
!     CALVI project team
!
!  This code SeLaLib (for Semi-Lagrangian-Library)
!  is a parallel library for simulating the plasma turbulence
!  in a tokamak.
!
!  This software is governed by the CeCILL-B license
!  under French law and abiding by the rules of distribution
!  of free software.  You can  use, modify and redistribute
!  the software under the terms of the CeCILL-B license as
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info".
!**************************************************************

!> @ingroup particle_methods
!> \brief
!> This module provides some routines for plotting during PIC simulations.

module sll_m_pic_visu

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

   use sll_m_ascii_io, only: &
      sll_s_ascii_file_create

   use sll_m_gnuplot, only: &
      sll_o_gnuplot_2d

   use sll_m_utilities, only: &
      sll_s_int2string, &
      sll_s_new_file_id

   use sll_m_xdmf, only: &
      sll_s_xdmf_corect2d_nodes

   implicit none

   public :: &
      sll_s_compute_df_cic, &
      sll_o_distribution_gnuplot, &
      sll_s_distribution_m4_gnuplot, &
      sll_s_distribution_tsc_gnuplot, &
      sll_s_distribution_xdmf, &
      sll_s_electricpotential_gnuplot_inline, &
      sll_s_energies_electrostatic_gnuplot_inline, &
      sll_o_particles_center_gnuplot, &
      sll_s_particles_center_gnuplot_inline, &
      sll_o_plot_format_points3d, &
      sll_s_plot_format_xmdv

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> plot particles centers with gnuplot
   interface sll_o_particles_center_gnuplot
      module procedure xv_particles_center_gnuplot
   end interface sll_o_particles_center_gnuplot

!> plot particles distribution with gnuplot
   interface sll_o_distribution_gnuplot
      module procedure distribution_xv_gnuplot
   end interface sll_o_distribution_gnuplot

!> write point3d file to plot particles characteristics
   interface sll_o_plot_format_points3d
      module procedure pq_plot_format_points3d
      module procedure pqr_plot_format_points3d
      module procedure pqrs_plot_format_points3d
   end interface sll_o_plot_format_points3d

contains

!> To plot particles run : gnuplot -persitant plot_name.gnu
   subroutine xv_particles_center_gnuplot(plot_name, &
                                          x, v, xmin, xmax, vmin, vmax, iplot, time)

      character(len=*), intent(in) :: plot_name
      sll_real64, dimension(:), intent(in) :: x
      sll_real64, dimension(:), intent(in) :: v
      sll_real64, intent(in) :: xmin
      sll_real64, intent(in) :: xmax
      sll_real64, intent(in) :: vmin
      sll_real64, intent(in) :: vmax
      sll_int32, intent(in) :: iplot
      sll_real64, optional   :: time

      character(len=4) :: fin
      sll_int32        :: file_id
      sll_int32        :: nbpart
      sll_int32        :: k

      call sll_s_int2string(iplot, fin)

      open (newunit=file_id, file=plot_name//'.gnu', position="append")
      if (iplot <= 1) then
         rewind (file_id)
         write (file_id, "('set xr[',g15.3,':',g15.3,']')") xmin, xmax
         write (file_id, "('set yr[',g15.3,':',g15.3,']')") vmin, vmax
      end if
      if (present(time)) then
         write (file_id, "(A18,G10.3,A1)") "set title 'Time = ", time, "'"
      end if
      write (file_id, "(a)") "plot '"//plot_name//"_"//fin//".dat' w d "
      close (file_id)

      nbpart = size(x)
      SLL_ASSERT(nbpart == size(v))
      open (newunit=file_id, file=plot_name//"_"//fin//'.dat')
      do k = 1, nbpart
         write (file_id, *) sngl(x(k)), sngl(v(k))
      end do
      close (file_id)

   end subroutine xv_particles_center_gnuplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine compute_df(df, x, v, xmin, xmax, nx, vmin, vmax, nv)
      sll_real64, dimension(:), intent(in) :: x, v
      sll_real64 :: delta_x, delta_v, xmin, xmax, vmin, vmax
      sll_int32, intent(in) :: nx
      sll_int32, intent(in) :: nv
      sll_real64, dimension(nx, nv) :: df
      sll_real64 :: weight
      sll_int32  :: i, j, k

      delta_x = (xmax - xmin)/(nx - 1)
      delta_v = (vmax - vmin)/(nv - 1)

      weight = 1._f64!/(delta_x*delta_v)   ! needs improvement

      df = 0.0_f64
      do k = 1, size(x)
         i = floor((x(k) - xmin)/delta_x) + 1
         j = floor((v(k) - vmin)/delta_v) + 1
         df(i, j) = df(i, j) + weight
      end do

   end subroutine compute_df

   subroutine distribution_xv_gnuplot(plot_name, x, v, xmin, xmax, nx, &
                                      vmin, vmax, nv, iplot, time)

      character(len=*), intent(in) :: plot_name
      sll_int32, intent(in) :: nx
      sll_int32, intent(in) :: nv
      sll_real64, dimension(:), intent(in) :: x, v
      sll_int32 :: iplot
      sll_real64, dimension(nx, nv) :: df
      sll_real64 :: time
      sll_real64 :: delta_x, delta_v, xmin, xmax, vmin, vmax
      character(len=4) :: fin
!sll_int32 :: file_id
      sll_int32 :: error

      delta_x = (xmax - xmin)/nx
      delta_v = (vmax - vmin)/nv
      call sll_s_int2string(iplot, fin)

      SLL_ASSERT(size(x) == size(v))

      write (*, "(A7,G10.3)") "Time = ", time
      call compute_df(df, x, v, xmin, xmax, nx, vmin, vmax, nv)

!call sll_s_new_file_id(file_id, error)
!open( file_id, file = plot_name//'.gnu', position="append" )
!if ( iplot == 1 ) rewind(file_id)
!write(file_id,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
!write(file_id,*)"splot  '"//plot_name//"_"//fin//".dat' w l"
!close(file_id)

      call sll_o_gnuplot_2d(xmin, xmax, nx, vmin, vmax, nv, df, plot_name, iplot, error)

   end subroutine distribution_xv_gnuplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Call this function inline with | gnuplot
   subroutine sll_s_particles_center_gnuplot_inline(x, v, xmin, xmax, ymin, ymax, time)
      sll_real64, dimension(:), intent(in) :: x, v
      sll_real64, intent(in) :: time
      sll_real64 :: xmin, xmax, ymin, ymax
      sll_int32  :: nbpart
      sll_int32  :: k

      print *, 'set term x11 2'
      print *, "set title 'time=", time, "'"
      print *, "unset logscale"
      nbpart = size(x)
      SLL_ASSERT(nbpart == size(v))

100   format('plot [', f7.3, ':', f7.3, '][', f7.3, ':', f7.3, '] ''-'' w d')

      write (*, 100) xmin, xmax, ymin, ymax
      do k = 1, nbpart
         print *, x(k), v(k)
      end do
      print *, 'e'

   end subroutine sll_s_particles_center_gnuplot_inline

!> point3D format http://www.visitusers.org/index.php?title=Reading_point_data
!> This format is designed to plot x,y,z particles with one weight (only for
!> characteristics)
!> This format is readable by VisIt
   subroutine pq_plot_format_points3d(plot_name, x, v, iplot)
      character(len=*), intent(in) :: plot_name
      sll_real64, dimension(:), intent(in) :: x, v
      sll_int32, intent(in) :: iplot
      character(len=4) :: fin
      sll_int32 :: file_id, error
      sll_int32 :: nbpart
      sll_int32  :: k

      call sll_s_int2string(iplot, fin)

      call sll_s_ascii_file_create(plot_name//"_"//fin//".3D", file_id, error)
      nbpart = size(x)
      SLL_ASSERT(size(v) == nbpart)
      write (file_id, "(a)") 'x v zero zero'
      do k = 1, nbpart
         write (file_id, "(4e15.3)") x(k), v(k), 0., 0.
      end do
      close (file_id)

   end subroutine pq_plot_format_points3d

!> point3D format http://www.visitusers.org/index.php?title=Reading_point_data
!> This format is designed to plot x,y,z particles
!> This format is readable by VisIt
   subroutine pqr_plot_format_points3d(plot_name, x, y, z, iplot)
      character(len=*), intent(in) :: plot_name
      sll_real64, dimension(:), intent(in) :: x, y, z
      sll_int32, intent(in) :: iplot
      character(len=4) :: fin
      sll_int32 :: file_id, error
      sll_int32 :: nbpart
      sll_int32  :: k

      call sll_s_int2string(iplot, fin)

      call sll_s_ascii_file_create(plot_name//"_"//fin//".3D", file_id, error)
      nbpart = size(x)
      SLL_ASSERT(size(y) == nbpart)
      SLL_ASSERT(size(z) == nbpart)
      write (file_id, "(a)") 'x y z zero'
      do k = 1, nbpart
         write (file_id, "(4e15.3)") x(k), y(k), z(k), 0.
      end do
      close (file_id)

   end subroutine pqr_plot_format_points3d

!> point3D format http://www.visitusers.org/index.php?title=Reading_point_data
!> This format is designed to plot x,y,z particles with one weight (only four
!> characteristics)
!> This format is readable by VisIt
   subroutine pqrs_plot_format_points3d(plot_name, p, q, r, s, iplot)
      character(len=*), intent(in) :: plot_name
      sll_real64, dimension(:), intent(in) :: p, q, r, s
      sll_int32, intent(in) :: iplot
      character(len=4) :: fin
      sll_int32 :: file_id, error
      sll_int32 :: nbpart
      sll_int32  :: k

      call sll_s_int2string(iplot, fin)

      call sll_s_ascii_file_create(plot_name//"_"//fin//".3D", file_id, error)
      nbpart = size(p)
      SLL_ASSERT(size(q) == nbpart)
      SLL_ASSERT(size(r) == nbpart)
      SLL_ASSERT(size(s) == nbpart)
      write (file_id, "(a)") 'x y z val'
      do k = 1, nbpart
         write (file_id, "(4e15.3)") p(k), q(k), r(k), s(k)
      end do
      close (file_id)

   end subroutine pqrs_plot_format_points3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> xmdv format http://davis.wpi.edu/xmdv/fileformats.html
!> This format is designed for particles with a lot of
!> characteristics. We have to give max and min values.
!> This format is readable by VisIt
   subroutine sll_s_plot_format_xmdv(plot_name, x, v, iplot, xmin, xmax, vmin, vmax)
      character(len=*), intent(in) :: plot_name
      sll_real64, dimension(:), intent(in) :: x
      sll_real64, dimension(:), intent(in) :: v
      sll_int32, intent(in) :: iplot
      character(len=4) :: fin
      sll_int32  :: file_id, error
      sll_real64 :: xmin, xmax, vmin, vmax
      sll_int32  :: nbpart
      sll_int32  :: k

      nbpart = size(x)

      SLL_ASSERT(size(v) == nbpart)
      call sll_s_int2string(iplot, fin)

      call sll_s_ascii_file_create(plot_name//"_"//fin//".okc", file_id, error)

      write (file_id, "(2i7)") 2, nbpart, 1   ! two characteristics, the number 1 is unused
      write (file_id, "(a)") 'x'
      write (file_id, "(a)") 'v'
      write (file_id, "(2f8.3,1x,i7)") xmin, xmax, 1
      write (file_id, "(2f8.3,1x,i7)") vmin, vmax, 1
      do k = 1, nbpart
         write (file_id, "(2e15.3)") x(k), v(k)
      end do
      close (file_id)

   end subroutine sll_s_plot_format_xmdv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>VisIt readable output for particles density
!>Data file format could be XML, HDF5 or Binary (not fully implemented yet)
   subroutine sll_s_distribution_xdmf(plot_name, x, v, w, &
                                      xmin, xmax, nx, &
                                      vmin, vmax, nv, iplot)
      character(len=*), intent(in) :: plot_name
      sll_real64, dimension(:), intent(in) :: x
      sll_real64, dimension(:), intent(in) :: v
      sll_real64, dimension(:), intent(in) :: w
      sll_int32, intent(in) :: nx
      sll_int32, intent(in) :: nv
      sll_int32 :: iplot
      sll_real64, dimension(nx, nv) :: df
      sll_real64 :: delta_x, delta_v, xmin, xmax, vmin, vmax
      character(len=4) :: fin

      call sll_s_int2string(iplot, fin)

      call sll_s_compute_df_cic(x, v, w, xmin, xmax, nx, vmin, vmax, nv, df)
      delta_x = (xmax - xmin)/(nx - 1)
      delta_v = (vmax - vmin)/(nv - 1)
      call sll_s_xdmf_corect2d_nodes(plot_name//'_'//fin, df, "density", xmin, delta_x, vmin, delta_v)

   end subroutine sll_s_distribution_xdmf

!>Compute grid field from particles distribution with the NGP scheme
!> (Nearest Grid Point)
   subroutine compute_df_ngp(xp, yp, wp, xmin, xmax, nx, ymin, ymax, ny, df)
      sll_real64, dimension(:), intent(in) :: xp, yp, wp
      sll_real64, intent(in) :: xmin, xmax, ymin, ymax
      sll_int32 :: ip, jp, kp
      sll_real64 :: xt, yt
      sll_int32 :: nx, ny
      sll_real64, dimension(nx, ny), intent(out) :: df
      sll_int32 :: nbpart

      nbpart = size(xp)
      SLL_ASSERT(nbpart == size(yp))
      SLL_ASSERT(nbpart == size(wp))
      df = 0.0_f64

      do kp = 1, nbpart

         xt = (xp(kp) - xmin)/(xmax - xmin)*nx
         yt = (yp(kp) - ymin)/(ymax - ymin)*ny

         ip = floor(xt)
         jp = floor(yt)

         SLL_ASSERT(ip >= 1 .and. ip <= nx .and. jp >= 1 .and. jp <= ny)

         df(ip, jp) = df(ip, jp) + wp(kp)

      end do

   end subroutine compute_df_ngp

!>Compute grid field from particles distribution with the CIC scheme (Cloud In
!Cell)
   subroutine sll_s_compute_df_cic(xp, yp, wp, xmin, xmax, nx, ymin, ymax, ny, df)
      sll_real64, dimension(:), intent(in) :: xp, yp, wp
      sll_real64, intent(in) :: xmin, xmax, ymin, ymax
      sll_int32 :: ip, jp, kp
      sll_real64 :: a1, a2, a3, a4, xt, yt
      sll_int32 :: nx, ny
      sll_real64, dimension(nx, ny), intent(out) :: df
      sll_int32 :: nbpart

      nbpart = size(xp)
      SLL_ASSERT(nbpart == size(yp))
      SLL_ASSERT(nbpart == size(wp))
      df = 0.0_f64

      do kp = 1, nbpart

         xt = (xp(kp) - xmin)/(xmax - xmin)*nx
         yt = (yp(kp) - ymin)/(ymax - ymin)*ny

         ip = floor(xt)
         jp = floor(yt)

         a1 = (ip + 1 - xt)*(jp + 1 - yt)
         a2 = (xt - ip)*(jp + 1 - yt)
         a3 = (ip + 1 - xt)*(yt - jp)
         a4 = (xt - ip)*(yt - jp)

         ip = ip + 1
         jp = jp + 1
         if (ip == nx + 1) ip = nx
         if (jp == ny + 1) jp = ny
!print *,ip, xt, jp, floor(yt), nx, ny
         SLL_ASSERT(ip > 0 .and. ip <= nx .and. jp > 0 .and. jp <= ny)

!   df(ip  ,jp  )=df(ip  ,jp  )+a1*wp(kp)
!   df(ip+1,jp  )=df(ip+1,jp  )+a2*wp(kp)
!   df(ip  ,jp+1)=df(ip  ,jp+1)+a3*wp(kp)
!   df(ip+1,jp+1)=df(ip+1,jp+1)+a4*wp(kp)

         if (jp > 0 .and. jp <= ny) then
            if (ip > 0 .and. ip <= nx) df(ip, jp) = df(ip, jp) + a1*wp(kp)
            if (ip + 1 > 0 .and. ip + 1 <= nx) df(ip + 1, jp) = df(ip + 1, jp) + a2*wp(kp)
         end if
         if (jp + 1 > 0 .and. jp + 1 <= ny) then
            if (ip > 0 .and. ip <= nx) df(ip, jp + 1) = df(ip, jp + 1) + a3*wp(kp)
            if (ip + 1 > 0 .and. ip + 1 <= nx) df(ip + 1, jp + 1) = df(ip + 1, jp + 1) + a4*wp(kp)
         end if
      end do

   end subroutine sll_s_compute_df_cic

   subroutine compute_df_m4(xp, yp, wp, xmin, xmax, nx, ymin, ymax, ny, df)
      sll_real64, intent(in) :: xp(:)
      sll_real64, intent(in) :: yp(:)
      sll_real64, intent(in) :: wp(:)
      sll_real64, intent(in) :: xmin
      sll_real64, intent(in) :: xmax
      sll_int32, intent(in) :: nx
      sll_real64, intent(in) :: ymin
      sll_real64, intent(in) :: ymax
      sll_int32, intent(in) :: ny
      sll_real64, intent(out) :: df(nx, ny)

      sll_int32  :: ip, jp, kp
      sll_real64 :: xt, x, cx, cm1x, cm2x, cp1x, cp2x
      sll_real64 :: yt, y, cy, cm1y, cm2y, cp1y, cp2y

! Initialize output array to zero
      df = 0.0_f64

      do kp = 1, size(xp)

         xt = (xp(kp) - xmin)/(xmax - xmin)*nx
         yt = (yp(kp) - ymin)/(ymax - ymin)*ny

         ip = floor(xt)
         jp = floor(yt)

         x = xt - ip
         y = yt - jp

         cm2x = f_m4(2.+x)
         cp2x = f_m4(2.-x)
         cm1x = f_m4(1.+x)
         cp1x = f_m4(1.-x)
         cx = f_m4(x)
         cy = f_m4(y)
         cm2y = f_m4(2.+y)
         cp2y = f_m4(2.-y)
         cm1y = f_m4(1.+y)
         cp1y = f_m4(1.-y)

         if (ip > 2 .and. jp > 2 .and. ip < nx - 1 .and. jp < ny - 1) then

            df(ip - 2, jp - 2) = df(ip - 2, jp - 2) + cm2x*cm2y*wp(kp)
            df(ip - 2, jp - 1) = df(ip - 2, jp - 1) + cm2x*cm1y*wp(kp)
            df(ip - 2, jp) = df(ip - 2, jp) + cm2x*cy*wp(kp)
            df(ip - 2, jp + 1) = df(ip - 2, jp + 1) + cm2x*cp1y*wp(kp)
            df(ip - 2, jp + 2) = df(ip - 2, jp + 2) + cm2x*cp2y*wp(kp)

            df(ip - 1, jp - 2) = df(ip - 1, jp - 2) + cm1x*cm2y*wp(kp)
            df(ip - 1, jp - 1) = df(ip - 1, jp - 1) + cm1x*cm1y*wp(kp)
            df(ip - 1, jp) = df(ip - 1, jp) + cm1x*cy*wp(kp)
            df(ip - 1, jp + 1) = df(ip - 1, jp + 1) + cm1x*cp1y*wp(kp)
            df(ip - 1, jp + 2) = df(ip - 1, jp + 2) + cm1x*cp2y*wp(kp)

            df(ip, jp - 2) = df(ip, jp - 2) + cx*cm2y*wp(kp)
            df(ip, jp - 1) = df(ip, jp - 1) + cx*cm1y*wp(kp)
            df(ip, jp) = df(ip, jp) + cx*cy*wp(kp)
            df(ip, jp + 1) = df(ip, jp + 1) + cx*cp1y*wp(kp)
            df(ip, jp + 2) = df(ip, jp + 2) + cx*cp2y*wp(kp)

            df(ip + 1, jp - 2) = df(ip + 1, jp - 2) + cp1x*cm2y*wp(kp)
            df(ip + 1, jp - 1) = df(ip + 1, jp - 1) + cp1x*cm1y*wp(kp)
            df(ip + 1, jp) = df(ip + 1, jp) + cp1x*cy*wp(kp)
            df(ip + 1, jp + 1) = df(ip + 1, jp + 1) + cp1x*cp1y*wp(kp)
            df(ip + 1, jp + 2) = df(ip + 1, jp + 2) + cp1x*cp2y*wp(kp)

            df(ip + 2, jp - 2) = df(ip + 2, jp - 2) + cp2x*cm2y*wp(kp)
            df(ip + 2, jp - 1) = df(ip + 2, jp - 1) + cp2x*cm1y*wp(kp)
            df(ip + 2, jp) = df(ip + 2, jp) + cp2x*cy*wp(kp)
            df(ip + 2, jp + 1) = df(ip + 2, jp + 1) + cp2x*cp1y*wp(kp)
            df(ip + 2, jp + 2) = df(ip + 2, jp + 2) + cp2x*cp2y*wp(kp)

         end if

      end do

   end subroutine compute_df_m4

!> M4 function from Monhagan (SPH method)
   sll_real64 function f_m4(x)
      sll_real64, intent(in) :: x
      if (x .gt. 2.) then
         f_m4 = 0._f64
      else if (x .ge. 1. .and. x .le. 2.) then
         f_m4 = 0.5*(2.-x)**2*(1.-x)
      else if (x .le. 1.) then
         f_m4 = 1.-2.5*x**2 + 1.5*(abs(x))**3
      end if

      return
   end function f_m4

!>GNUplot readable output for particles density
   subroutine sll_s_distribution_m4_gnuplot(plot_name, x, v, w, &
                                            xmin, xmax, nx, &
                                            vmin, vmax, nv, iplot)
      character(len=*), intent(in) :: plot_name
      sll_real64, intent(in) :: x(:)
      sll_real64, intent(in) :: v(:)
      sll_real64, intent(in) :: w(:)
      sll_real64, intent(in) :: xmin
      sll_real64, intent(in) :: xmax
      sll_int32, intent(in) :: nx
      sll_real64, intent(in) :: vmin
      sll_real64, intent(in) :: vmax
      sll_int32, intent(in) :: nv
      sll_int32, intent(in) :: iplot

      sll_int32        :: error
      sll_real64       :: df(nx, nv)
      character(len=4) :: fin

      call sll_s_int2string(iplot, fin)

      call compute_df_m4(x, v, w, xmin, xmax, nx, vmin, vmax, nv, df)

      call sll_o_gnuplot_2d(xmin, xmax, nx, vmin, vmax, nv, df, plot_name, iplot, error)

   end subroutine sll_s_distribution_m4_gnuplot

#define FONCTION1( X ) (0.75_f64-(X)*(X))

#define FONCTION2( X ) (0.5_f64 * ( 1.5_f64 - (X) )**2)

#define BSPLINES(X,Y) \
   c_1x = FONCTION2(1 + X); \
   c1x = FONCTION1(X); \
   c2x = FONCTION2(1 - X); \
   c_1y = FONCTION2(1 + Y); \
   c1y = FONCTION1(Y); \
   c2y = FONCTION2(1 - Y)

   subroutine compute_df_tsc(xp, yp, wp, xmin, xmax, nx, ymin, ymax, ny, df)
      sll_real64, dimension(:), intent(in) :: xp, yp, wp
      sll_real64, intent(in) :: xmin, xmax, ymin, ymax
      sll_int32 :: ip, jp, kp
      sll_real64 :: xt, yt
      sll_int32 :: nx, ny
      sll_real64, dimension(nx, ny), intent(out) :: df
      sll_real64 :: x, c1x, c_1x, c2x
      sll_real64 :: y, c1y, c_1y, c2y

      df = 0.0_f64

      do kp = 1, size(xp)

         xt = (xp(kp) - xmin)/(xmax - xmin)*nx
         yt = (yp(kp) - ymin)/(ymax - ymin)*ny

         ip = floor(xt)
         jp = floor(yt)

         x = xt - ip
         y = yt - jp

         BSPLINES(x, y)

         if (ip > 1 .and. ip < nx .and. jp > 1 .and. jp < ny) then

            df(ip, jp) = df(ip, jp) + c1x*c1y*wp(kp)
            df(ip, jp - 1) = df(ip, jp - 1) + c1x*c_1y*wp(kp)
            df(ip, jp + 1) = df(ip, jp + 1) + c1x*c2y*wp(kp)
            df(ip - 1, jp) = df(ip - 1, jp) + c_1x*c1y*wp(kp)
            df(ip - 1, jp - 1) = df(ip - 1, jp - 1) + c_1x*c_1y*wp(kp)
            df(ip - 1, jp + 1) = df(ip - 1, jp + 1) + c_1x*c2y*wp(kp)
            df(ip + 1, jp) = df(ip + 1, jp) + c2x*c1y*wp(kp)
            df(ip + 1, jp - 1) = df(ip + 1, jp - 1) + c2x*c_1y*wp(kp)
            df(ip + 1, jp + 1) = df(ip + 1, jp + 1) + c2x*c2y*wp(kp)

         end if

      end do

   end subroutine compute_df_tsc

!>GNUplot readable output for particles density
   subroutine sll_s_distribution_tsc_gnuplot(plot_name, x, v, w, &
                                             xmin, xmax, nx, &
                                             vmin, vmax, nv, iplot)
      character(len=*), intent(in) :: plot_name
      sll_real64, dimension(:), intent(in) :: x
      sll_real64, dimension(:), intent(in) :: v
      sll_real64, dimension(:), intent(in) :: w
      sll_int32, intent(in) :: nx
      sll_int32, intent(in) :: nv
      sll_int32 :: iplot, error
      sll_real64, dimension(nx, nv) :: df
      sll_real64 :: xmin, xmax, vmin, vmax
      character(len=4) :: fin

      call sll_s_int2string(iplot, fin)

      call compute_df_tsc(x, v, w, xmin, xmax, nx, vmin, vmax, nv, df)

      call sll_o_gnuplot_2d(xmin, xmax, nx, vmin, vmax, nv, df, plot_name, iplot, error)

   end subroutine sll_s_distribution_tsc_gnuplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Call this function inline with | gnuplot
!> Plots different energies for Electrostatic Particle-in-Cell
!> Pseudo/Virtual file plot in gnuplot using plot "-"
   subroutine sll_s_energies_electrostatic_gnuplot_inline(kinetic_energy, &
                                                          electrostatic_energy, impulse, timestepwidth)
      sll_real64, dimension(:), intent(in) :: kinetic_energy, electrostatic_energy, impulse
      sll_int32  :: timesteps, k
      sll_real64, dimension(size(kinetic_energy)) :: total_energy
      sll_real64, intent(in) :: timestepwidth

      timesteps = size(kinetic_energy)
      SLL_ASSERT(timesteps == size(electrostatic_energy))
      SLL_ASSERT(timesteps == size(total_energy))

      total_energy = electrostatic_energy + kinetic_energy

      if (timesteps > 2) then
         print *, "   "
!!Plot Energy errors
         print *, "set term x11 3"
         print *, "set logscale y"
         print *, "set title 'Relative Energies'"
         print *, "plot '-'  title 'kinetic' with lines,\ "
         if (maxval(electrostatic_energy) /= 0.0_f64) print *, "'-'  title 'electrostatic' with lines,\ "
         print *, "'-'   title 'total' with lines;"
         do k = 2, timesteps
            print *, (k - 1)*timestepwidth, kinetic_energy(k)/maxval(kinetic_energy)
         end do
         print *, "e"
         if (maxval(electrostatic_energy) /= 0.0_f64) then
         do k = 2, timesteps
            print *, (k - 1)*timestepwidth, electrostatic_energy(k)/maxval(electrostatic_energy)
         end do
         end if
         print *, "e"
         do k = 2, timesteps
            print *, (k - 1)*timestepwidth, total_energy(k)/maxval(total_energy)
         end do
         print *, "e"
         print *, "   "

!Plot Energy errors
         print *, "set term x11 4"
         print *, "set autoscale y"
         print *, "set title 'Total Energy Error'"
         print *, "set logscale y"

         print *, "plot  '-' using 1:2  title 'Relative Error' with lines"
         do k = 2, timesteps
            print *, (k - 1)*timestepwidth, abs(((total_energy(k) - total_energy(1))/total_energy(1)))
         end do
         print *, "e"
         print *, "   "

!!plot Impulse Error
         print *, "set term x11 6"
         print *, "set autoscale y"
         print *, "set title 'Relative Impulse Error'"
         print *, "set logscale y"
         print *, "plot '-' using 1:2 title 'Relative Impulse Error' with lines"
         do k = 2, timesteps
            if (impulse(1) /= 0) then
               print *, (k - 1)*timestepwidth, abs(impulse(k) - impulse(1))/abs(impulse(1))
            else if (maxval(abs(impulse)) /= 0.0_f64) then
               print *, (k - 1)*timestepwidth, abs(impulse(k) - impulse(1))/maxval(abs(impulse))
            else
               print *, (k - 1)*timestepwidth, abs(impulse(k) - impulse(1))
            end if
         end do
         print *, "e"

!!!plot Electrostatic energy to see landau damping
!print*, "set term x11 5"
!print*, "set autoscale y"
!print*, "set title 'Relative Electrostatic Energy'"
!print*, "set logscale y"
!print*, "plot  '-' using 1:2  title 'Electrostatic Energy' with lines"
!do k = 2, timesteps
!   print*, k,  electrostatic_energy(k)
!end do
!print*, "e"

      end if
   end subroutine sll_s_energies_electrostatic_gnuplot_inline

!> Call this function inline with | gnuplot
!> Plots different energies for Electrostatic Particle-in-Cell
!> Pseudo/Virtual file plot in gnuplot using plot "-"
   subroutine sll_s_electricpotential_gnuplot_inline(values, eval_knots)
      sll_real64, dimension(:), intent(in) :: values, eval_knots
      sll_int32  :: nknots, k

      nknots = size(eval_knots)
      SLL_ASSERT(nknots == size(eval_knots))

!Plot Energy errors
      print *, "set term x11 7"
      print *, "set autoscale y"
      print *, "set title 'Electric Potential'"
      print *, "unset logscale y"
      print *, "plot  '-' using 1:2  title 'Potential' with lines"
      do k = 1, nknots
         print *, eval_knots(k), values(k)
      end do
      print *, "e"

   end subroutine sll_s_electricpotential_gnuplot_inline

end module sll_m_pic_visu
