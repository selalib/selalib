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

!> \author Pierre Navaro
!> \brief  
!> This module provides some routines for plotting during PIC simulations.
module sll_visu_pic
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_file_io.h"

implicit none

!> plot particles centers with gnuplot
interface particles_center_gnuplot
   module procedure xv_particles_center_gnuplot
end interface particles_center_gnuplot

!> plot particles distribution with gnuplot
interface distribution_gnuplot
   module procedure distribution_xv_gnuplot
end interface distribution_gnuplot

!> write point3d file to plot particles characteristics
interface plot_format_points3d
   module procedure pq_plot_format_points3d
   module procedure pqr_plot_format_points3d
   module procedure pqrs_plot_format_points3d
end interface plot_format_points3d

sll_int32, private :: i, j, k

contains

!> To plot particles run : gnuplot -persitant plot_name.gnu
subroutine xv_particles_center_gnuplot( plot_name, x, v, xmin, xmax, vmin, vmax, iplot, time )
character(len=*), intent(in) :: plot_name
sll_real64, dimension(:), intent(in) :: x, v
sll_int32 :: iplot, error
sll_real64 :: time, xmin, xmax, vmin, vmax
character(len=4) :: fin
sll_int32 :: file_id
sll_int32 :: nbpart

call int2string(iplot, fin)

call sll_new_file_id(file_id, error)
open( file_id, file = plot_name//'.gnu', position="append" )
if ( iplot <= 1 ) then
   rewind(file_id)
   write(file_id,"('set xr[',g15.3,':',g15.3,']')") xmin, xmax
   write(file_id,"('set yr[',g15.3,':',g15.3,']')") vmin, vmax
end if
write(file_id,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
write(file_id,"(a)")"plot '"//plot_name//"_"//fin//".dat' w d "
close(file_id)

nbpart = size(x)
SLL_ASSERT(nbpart == size(v))
open(file_id, file = plot_name//"_"//fin//'.dat' )
do i=1,nbpart
   write(file_id,*) sngl(x(i)),sngl(v(i)) 
end do
close(file_id)

end subroutine xv_particles_center_gnuplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_df( df, x, v, xmin, xmax, nx, vmin, vmax, nv)
sll_real64, dimension(:), intent(in) :: x, v
sll_real64 :: delta_x, delta_v, xmin, xmax, vmin, vmax
sll_int32, intent(in) :: nx
sll_int32, intent(in) :: nv
sll_real64, dimension(nx,nv) :: df
sll_real64 :: weight

delta_x = (xmax-xmin)/(nx-1)
delta_v = (vmax-vmin)/(nv-1)

weight = 1._f64!/(delta_x*delta_v)   ! needs improvement

df = 0.d0
do k=1,size(x)
   do i=1,nx
      if (xmin+(i-1)*delta_x <= x(k) .and. x(k) < xmin+i*delta_x) then
         do j=1,nv
            if (vmin+(j-1)*delta_v <= v(k) .and. v(k) < vmin+j*delta_v) then
               df(i,j) = df(i,j) + weight
            endif
         enddo
      end if
   enddo
enddo

end subroutine compute_df

subroutine distribution_xv_gnuplot( plot_name, x, v, xmin, xmax, nx, &
                                    vmin, vmax, nv, iplot, time)  

character(len=*), intent(in) :: plot_name
sll_int32, intent(in) :: nx
sll_int32, intent(in) :: nv
sll_real64, dimension(:), intent(in) :: x, v
sll_int32 :: iplot
sll_real64, dimension(nx,nv) :: df
sll_real64 :: time
sll_real64 :: delta_x, delta_v, xmin, xmax, vmin, vmax
character(len=4) :: fin
sll_int32 :: file_id, error

delta_x = (xmax-xmin)/nx
delta_v = (vmax-vmin)/nv
call int2string(iplot,fin)

SLL_ASSERT(size(x)==size(v))

call compute_df( df, x, v, xmin, xmax, nx, vmin, vmax, nv)

call sll_new_file_id(file_id, error)
open( file_id, file = plot_name//'.gnu', position="append" )
if ( iplot == 1 ) rewind(file_id)
write(file_id,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
write(file_id,*)"splot  '"//plot_name//"_"//fin//".dat' w l"
close(file_id)

call sll_gnuplot_2d(xmin, xmax, nx, vmin, vmax, nv, df, plot_name, iplot, error)  

end subroutine distribution_xv_gnuplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Call this function inline with | gnuplot
subroutine particles_center_gnuplot_inline( x, v, xmin, xmax, ymin, ymax, time )
sll_real64, dimension(:), intent(in) :: x, v
sll_real64, intent(in) :: time
sll_real64 :: xmin, xmax, ymin, ymax
sll_int32  :: nbpart

print*, "set title'time=",time,"'"
print*, 'set term x11'
nbpart = size(x)
SLL_ASSERT( nbpart == size(v))
write(*,100) xmin,xmax,ymin,ymax
do k = 1, nbpart
   print*, x(k), v(k)
end do
print*, 'e'

100 format('p [',f7.3,':',f7.3,'][',f7.3,':',f7.3,'] ''-'' w d')

end subroutine particles_center_gnuplot_inline

!> point3D format http://www.visitusers.org/index.php?title=Reading_point_data 
!> This format is designed to plot x,y,z particles with one weight (only four
!> characteristics)
!> This format is readable by VisIt
subroutine pq_plot_format_points3d( plot_name, x, v, iplot)
character(len=*), intent(in) :: plot_name
sll_real64, dimension(:), intent(in) :: x, v
sll_int32, intent(in) :: iplot
character(len=4) :: fin
sll_int32 :: file_id, error
sll_int32 :: nbpart

call int2string(iplot, fin)

call sll_ascii_file_create(plot_name//"_"//fin//".3D", file_id, error)
nbpart = size(x)
SLL_ASSERT( size(v) == nbpart)
write(file_id,"(a)") 'x v zero zero'
do k = 1, nbpart
  write(file_id,"(4e15.3)")x(k),v(k), 0., 0.
end do
close(file_id)

end subroutine pq_plot_format_points3d

!> point3D format http://www.visitusers.org/index.php?title=Reading_point_data 
!> This format is designed to plot x,y,z particles with one weight (only four
!> characteristics)
!> This format is readable by VisIt
subroutine pqr_plot_format_points3d( plot_name, x, y, z, iplot)
character(len=*), intent(in) :: plot_name
sll_real64, dimension(:), intent(in) :: x, y, z
sll_int32, intent(in) :: iplot
character(len=4) :: fin
sll_int32 :: file_id, error
sll_int32 :: nbpart

call int2string(iplot, fin)

call sll_ascii_file_create(plot_name//"_"//fin//".3D", file_id, error)
nbpart = size(x)
SLL_ASSERT( size(y) == nbpart)
SLL_ASSERT( size(z) == nbpart)
write(file_id,"(a)") 'x y z zero'
do k = 1, nbpart
  write(file_id,"(4e15.3)")x(k),y(k),z(k), 0.
end do
close(file_id)

end subroutine pqr_plot_format_points3d

!> point3D format http://www.visitusers.org/index.php?title=Reading_point_data 
!> This format is designed to plot x,y,z particles with one weight (only four
!> characteristics)
!> This format is readable by VisIt
subroutine pqrs_plot_format_points3d( plot_name, p, q, r, s, iplot)
character(len=*), intent(in) :: plot_name
sll_real64, dimension(:), intent(in) :: p, q, r, s
sll_int32, intent(in) :: iplot
character(len=4) :: fin
sll_int32 :: file_id, error
sll_int32 :: nbpart

call int2string(iplot, fin)

call sll_ascii_file_create(plot_name//"_"//fin//".3D", file_id, error)
nbpart = size(p)
SLL_ASSERT( size(q) == nbpart)
SLL_ASSERT( size(r) == nbpart)
SLL_ASSERT( size(s) == nbpart)
write(file_id,"(a)") 'x y z val'
do k = 1, nbpart
  write(file_id,"(4e15.3)")p(k),q(k),r(k),s(k)
end do
close(file_id)

end subroutine pqrs_plot_format_points3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> xmdv format http://davis.wpi.edu/xmdv/fileformats.html
!> This format is designed for particles with a lot of
!> characteristics. We have to give max and min values.
!> This format is readable by VisIt
subroutine plot_format_xmdv( plot_name, x, v, iplot, xmin, xmax, vmin, vmax)
character(len=*), intent(in) :: plot_name
sll_real64, dimension(:), intent(in) :: x
sll_real64, dimension(:), intent(in) :: v
sll_int32, intent(in) :: iplot
character(len=4) :: fin
sll_int32 :: file_id, error
sll_real64 :: xmin, xmax, vmin, vmax
sll_int32 :: nbpart

nbpart = size(x)

SLL_ASSERT( size(v) == nbpart)
call int2string(iplot, fin)

call sll_ascii_file_create(plot_name//"_"//fin//".okc", file_id, error)

write(file_id,"(2i7)") 2, nbpart, 1   ! two characteristics, the number 1 is unused
write(file_id,"(a)") 'x'
write(file_id,"(a)") 'v'
write(file_id,"(2f8.3,1x,i7)") xmin, xmax, 1
write(file_id,"(2f8.3,1x,i7)") vmin, vmax, 1
do k = 1, nbpart
  write(file_id,"(2e15.3)")x(k),v(k)
end do
close(file_id)

end subroutine plot_format_xmdv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>VisIt readable output for particles density
!>Data file format could be XML, HDF5 or Binary (not fully implemented yet)
subroutine distribution_xdmf(plot_name, x, v, w, xmin, xmax, nx, vmin, vmax, nv, iplot)  
character(len=*), intent(in) :: plot_name
sll_real64, dimension(:), intent(in) :: x
sll_real64, dimension(:), intent(in) :: v
sll_real64, dimension(:), intent(in) :: w
sll_int32, intent(in) :: nx
sll_int32, intent(in) :: nv
sll_int32 :: iplot
sll_real64, dimension(nx,nv) :: df
sll_real64 :: delta_x, delta_v, xmin, xmax, vmin, vmax
character(len=4) :: fin

call int2string(iplot, fin)

call compute_df_cic(x, v, w, xmin, xmax, nx, vmin, vmax, nv, df)  
delta_x = (xmax-xmin)/(nx-1)
delta_v = (vmax-vmin)/(nv-1)
call sll_xdmf_corect2d_nodes( plot_name//'_'//fin, df, "density", xmin, delta_x, vmin, delta_v) 

end subroutine distribution_xdmf

!>Compute grid field from particles distribution with the CIC scheme (Cloud In
!Cell)
subroutine compute_df_cic(xp, yp, wp, xmin, xmax, nx, ymin, ymax, ny, df)  
sll_real64, dimension(:), intent(in) :: xp, yp, wp
sll_real64, intent(in) :: xmin, xmax, ymin, ymax
sll_int32 :: ip, jp, kp
sll_real64 :: a1, a2, a3, a4, xt, yt
sll_int32 :: nx, ny
sll_real64, dimension(nx,ny), intent(out) :: df
sll_int32 :: nbpart

nbpart = size(xp)
SLL_ASSERT(nbpart == size(yp))
SLL_ASSERT(nbpart == size(wp))
df = 0.0

do kp = 1,nbpart

   xt = (xp(kp)-xmin)/(xmax-xmin)*nx
   yt = (yp(kp)-ymin)/(ymax-ymin)*ny

   ip = floor(xt)
   jp = floor(yt)

   SLL_ASSERT(ip > 0 .and. ip < nx .and. jp > 0 .and. jp < ny)

   a1 = (ip+1 - xt) * (jp+1 - yt)
   a2 = (xt   - ip) * (jp+1 - yt)
   a3 = (ip+1 - xt) * (yt   - jp)
   a4 = (xt   - ip) * (yt   - jp)

   df(ip  ,jp  )=df(ip  ,jp  )+a1*wp(kp)
   df(ip+1,jp  )=df(ip+1,jp  )+a2*wp(kp)
   df(ip  ,jp+1)=df(ip  ,jp+1)+a3*wp(kp)
   df(ip+1,jp+1)=df(ip+1,jp+1)+a4*wp(kp)

end do

end subroutine compute_df_cic

end module sll_visu_pic
