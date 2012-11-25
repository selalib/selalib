!> \file sll_visu_pic.F90
!> \namespace sll_visu_pic
!> \authors                    
!> Pierre Navaro (navaro@math.unistra.fr) 
!> \brief  
!> This module provides some routines for plotting during PIC simulations.
!>
module sll_visu_pic
#include "sll_working_precision.h"
#include "sll_assert.h"
use sll_io
use sll_ascii_io
use sll_gnuplot
use sll_xdmf

implicit none

interface particles_center_gnuplot
   module procedure xv_particles_center_gnuplot
end interface particles_center_gnuplot

interface distribution_gnuplot
   module procedure distribution_xv_gnuplot
end interface distribution_gnuplot

sll_int32, private :: i, j, k

contains

subroutine xv_particles_center_gnuplot( x, v, xmin, xmax, vmin, vmax, iplot, time )
sll_real64, dimension(:), intent(in) :: x, v
sll_int32 :: iplot, error
sll_real64 :: time, xmin, xmax, vmin, vmax
character(len=4) :: fin
sll_int32 :: file_id
sll_int32 :: nbpart

call int2string(iplot, fin)

call sll_new_file_id(file_id, error)
open( file_id, file = 'pxv.gnu', position="append" )
if ( iplot <= 1 ) then
   rewind(file_id)
   write(file_id,"('set xr[',g15.3,':',g15.3,']')") xmin, xmax
   write(file_id,"('set yr[',g15.3,':',g15.3,']')") vmin, vmax
end if
write(file_id,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
write(file_id,"(a)")"plot 'pxv_"//fin//"' w d "
close(file_id)

nbpart = size(x)
SLL_ASSERT(nbpart == size(v))
open(file_id, file = 'pxv_'//fin )
do k=1,nbpart
   write(file_id,*) sngl(x),sngl(v) 
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

delta_x = (xmax-xmin)/nx
delta_v = (vmax-vmin)/nv

weight = 1._f64/(delta_x*delta_v)   ! needs improvement

df = 0.d0
do k=1,size(x)
   do i=1,nx
      if (xmin+(i-1)*delta_x <= x(k) .and. x(k) < xmin+i*delta_x) then
         do j=1,nv
            if (vmin+(i-1)*delta_v <= v(k) .and. v(k) < vmin+i*delta_v) then
               df(i,j) = df(i,j) + weight
            endif
         enddo
      end if
   enddo
enddo

end subroutine compute_df

subroutine distribution_xv_gnuplot( x, v, xmin, xmax, nx, &
                                    vmin, vmax, nv, iplot, time)  

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
open( file_id, file = 'df_xv.gnu', position="append" )
if ( iplot == 1 ) rewind(file_id)
write(file_id,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
write(file_id,*)"splot  'df_xv_"//fin//".dat' w l"
close(file_id)

call sll_gnuplot_rect_2d(xmin, xmax, nx, vmin, vmax, nv, df, 'df_xv', iplot, error)  

end subroutine distribution_xv_gnuplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Call this function inline with | gnuplot
subroutine xv_particles_center_gnuplot_inline( x, v, xmin, xmax, ymin, ymax, time )
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

end subroutine xv_particles_center_gnuplot_inline

!> point3D format http://www.visitusers.org/index.php?title=Reading_point_data 
!> This format is designed to plot x,y,z particles with one weight (only four
!> characteristics)
subroutine plot_particles_points3d( x, v, iplot)
sll_real64, dimension(:), intent(in) :: x, v
sll_int32, intent(in) :: iplot
character(len=4) :: fin
sll_int32 :: file_id, error
sll_int32 :: nbpart

call int2string(iplot, fin)

call sll_ascii_file_create("particles_"//fin//".3D", file_id, error)
nbpart = size(x)
SLL_ASSERT( size(v) == nbpart)
write(file_id,"(a)") 'x v val1 val2'
do k = 1, nbpart
  write(file_id,"(4e15.3)")x(k),v(k), 1., 1.
end do
close(file_id)

end subroutine plot_particles_points3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> xmdv format http://davis.wpi.edu/xmdv/fileformats.html
!> This format is designed for particles with a lot of
!> characteristics. We have to give max and min values.
subroutine plot_particles_xmdv( x, v, iplot, xmin, xmax, vmin, vmax)
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

call sll_ascii_file_create("particles_"//fin//".okc", file_id, error)

write(file_id,"(2i7)") 2, nbpart, 1   ! two characteristics, the number 1 is unused
write(file_id,"(a)") 'x'
write(file_id,"(a)") 'v'
write(file_id,"(2f8.3,1x,i7)") xmin, xmax, 1
write(file_id,"(2f8.3,1x,i7)") vmin, vmax, 1
do k = 1, nbpart
  write(file_id,"(2e15.3)")x(k),v(k)
end do
close(file_id)

end subroutine plot_particles_xmdv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_particle_density( x, v, xmin, xmax, nx, vmin, vmax, nv, iplot)  
sll_real64, dimension(:), intent(in) :: x
sll_real64, dimension(:), intent(in) :: v
sll_int32, intent(in) :: nx
sll_int32, intent(in) :: nv
sll_int32 :: iplot
sll_real64, dimension(nx,nv) :: df
sll_real64 :: delta_x, delta_v, xmin, xmax, vmin, vmax
character(len=4) :: fin
sll_int32 :: error

call int2string(iplot, fin)

call compute_df( df, x, v, xmin, xmax, nx, vmin, vmax, nv)
call sll_gnuplot_rect_2d(xmin, xmax, nx, vmin, vmax, nv, df, 'density', iplot, error)  
call sll_xdmf_corect2d_nodes( 'df_'//fin, df, "density", xmin, delta_x, vmin, delta_v) 

end subroutine plot_particle_density

end module sll_visu_pic
