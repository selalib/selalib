module diagno
#include "sll_working_precision.h"
use sll_m_pic_visu
use sll_m_utilities
use sll_m_ascii_io
use sll_m_gnuplot
use sll_m_xdmf
use particules

implicit none

sll_int32, private :: i, j, k

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_phases( ele, iplot, time )

type(particle) :: ele
sll_int32 :: iplot, k
sll_real64 :: time
!sll_real64 :: gama, aux
sll_real64 :: speed
character(len=4) :: fin

call sll_s_int2string(iplot, fin)

open( 11, file = 'part.gnu', position="append" )
open( 12, file = 'xvx.gnu' , position="append" )
open( 13, file = 'yvy.gnu' , position="append" )

if ( iplot .eq. 1 ) then
   rewind(11)
   write(11,*)"set xr[0:",sngl(dimx),"]"
   write(11,*)"set yr[0:",sngl(dimy),"]"
   rewind(12)
   write(12,*)"set xr[0:",sngl(dimx),"]"
   write(12,*)"set yr[-1:1]"
   rewind(13)
   write(13,*)"set xr[0:",sngl(dimx),"]"
   write(13,*)"set yr[-1:1]"
end if

write(11,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
write(11,*)"plot 'part_"//fin//"' w d "
close(11)
write(12,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
write(12,*)"plot 'part_"//fin//"' u 1:3 w d "
close(12)
write(13,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
write(13,*)"plot 'part_"//fin//"' u 2:4 w d "
close(13)

open( 14, file = 'part_'//fin )
do k=1,nbpart
   speed = sqrt( ele%vpx(k)*ele%vpx(k) + &
        &        ele%vpy(k)*ele%vpy(k) )
   write(14,*)  sngl(ele%idx(k)*dx+ele%dpx(k))  &
              , sngl(ele%idy(k)*dy+ele%dpy(k))  &
              , sngl(ele%vpx(k))                &
              , sngl(ele%vpy(k))                &
              , sngl(speed)
end do
close(14)

end subroutine plot_phases

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Calcul de la fct de distribution en fonction de (vx,vy)

subroutine distribution_v(ele, iplot, time)  

sll_int32 :: i, k, iplot
sll_real64, dimension(:,:), allocatable :: df
type(particle) :: ele
sll_real64 :: time, vx, vy, aux, vth=1.0_f64
sll_real64 :: delta_v, vmin, vmax
character(len=4) :: fin
sll_int32, parameter :: nv = 64

call sll_s_int2string(iplot,fin)

allocate(df(nv,nv))

vmin = -6.d0
vmax = +6.d0
delta_v = (vmax-vmin)/nv

df = 0.d0
do k=1,nbpart
  i = floor(ele%vpx(k)/(vmax-vmin)*nv)
  j = floor(ele%vpy(k)/(vmax-vmin)*nv)
  df(i,j) = df(i,j) + ele%p(k)
enddo

open( 27, file = 'df_v.gnu', position="append" )
if ( iplot .eq. 1 ) rewind(27)
write(27,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"

open( 28, file = 'df_v_'//fin )

write(27,*)"splot  'df_v_"//fin//"' w l, 'df_theo' w l"
write(27,*)"pause 1"

do i=1,nv
  do j=1,nv
    vx = vmin+(i-0.5)*delta_v
    vy = vmin+(j-0.5)*delta_v
    write(28,*) vx,vy,df(i,j)/(delta_v*delta_v)
  end do
  write(28,*)
enddo

close(27)
close(28)

if ( iplot == 1 ) then

  open( 37, file = 'df_theo.gnu' )
  rewind(37)
  write(37,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
  write(37,"(A)")"splot  'df_theo' w l"
  close(37)

  open( 38, file = 'df_theo')
  do i=1,nv
    do j=1,nv
      vx = vmin+(i-0.5)*delta_v
      vy = vmin+(j-0.5)*delta_v
      aux = exp(-(vx*vx+vy*vy)/(2*vth*vth))/(2*pi*vth*vth)
      write(38,*) vx, vy, aux*dimx*dimy 
    end do
    write(38,*)
  enddo
  close(38)

endif

end subroutine distribution_v

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Calcul de la fct de distribution en fonction de (x,y)

subroutine distribution_x(ele, iplot, time)  

sll_int32 :: i, k, iplot
sll_real64, dimension(100,100) :: df
type(particle) :: ele
sll_real64 :: time, x, y, delta_x, delta_y
character(len=4) :: fin
sll_real64 :: ppx
sll_real64 :: ppy

call sll_s_int2string(iplot, fin)

df = 0.d0
delta_x = dimx/100
delta_y = dimy/100
do k=1,nbpart
  ppx =(ele%idx(k)+ele%dpx(k))*dx!  ele%idx(k)*dx+ele%dpx(k)
  ppy =(ele%idy(k)+ele%dpy(k))*dy!  ele%idy(k)*dy+ele%dpy(k)
  i   = floor(ppx/dimx*100)
  j   = floor(ppy/dimx*100)
  df(i,j) = df(i,j) + ele%p(k)
enddo


open( 27, file = 'df_x.gnu', position="append" )
if ( iplot .eq. 1 ) rewind(27)
write(27,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"

open( 28, file = 'df_x_'//fin )

write(27,*)"splot  'df_x_"//fin//"' w l"
write(27,*)"pause 1"

do i=1,100  
   do j=1,100 
      x = (i-0.5)*delta_x
      y = (j-0.5)*delta_y
      write(28,*) x,y,df(i,j)  
   end do
   write(28,*)
enddo

close(27)
close(28)

!call sll_o_gnuplot_2d(0._f64, dimx, 40, 0._f64, dimy, 40, df, 'density', iplot, error)  
!call sll_s_xdmf_corect2d_nodes('df_'//fin,df,"density",x_min,delta_x,y_min,delta_y) 

end subroutine distribution_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine modeE( tm, iplot, time )

type(tm_mesh_fields) :: tm
sll_real64 :: time, aux
sll_int32 :: iplot

aux =0.d0
do i=0,nx-1
  do j=0,ny-1
    aux = aux + tm%ex(i,j)*tm%ex(i,j)*dx*dy
  end do
end do
aux = 0.5*log(aux)

open(34,file='modeE.dat',position="append")
if (iplot==1) rewind(34)
write(34,*) time, aux
close(34)

end subroutine modeE

!> point3D http://www.visitusers.org/index.php?title=Reading_point_data 
subroutine plot_particles_points3d( p, iplot)

type(particle), intent(in) :: p
sll_int32, intent(in) :: iplot
character(len=4) :: fin
sll_int32 :: file_id, error

call sll_s_int2string(iplot, fin)

call sll_s_ascii_file_create("particles_"//fin//".3D", file_id, error)

write(file_id,"(a)") 'x y vx vy'
do k = 1, nbpart
  write(file_id,"(4e15.3)")p%idx(k)*dx+p%dpx(k), &
                           p%idy(k)*dy+p%dpy(k), &
                           p%vpx(k), p%vpy(k)
end do
close(file_id)

end subroutine plot_particles_points3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> xmdv format http://davis.wpi.edu/xmdv/fileformats.html
subroutine plot_particles_xmdv( p, iplot, xmin, xmax, ymin, ymax)

type(particle), intent(in) :: p
sll_int32, intent(in) :: iplot
character(len=4) :: fin
sll_int32 :: file_id, error
sll_real64 :: xmin, xmax, ymin, ymax
sll_real64 :: vxmin, vxmax, vymin, vymax
sll_real64 :: pmin, pmax

call sll_s_int2string(iplot, fin)

call sll_s_ascii_file_create("particles_"//fin//".okc", file_id, error)

vxmin = minval(p%vpx)
vxmax = maxval(p%vpx)
vymin = minval(p%vpy)
vymax = maxval(p%vpy)
pmin  = minval(p%p)
pmax  = maxval(p%p)

write(file_id,"(2i7)") 5, nbpart, 1
write(file_id,"(a)") 'x'
write(file_id,"(a)") 'y'
write(file_id,"(a)") 'vx'
write(file_id,"(a)") 'vy'
write(file_id,"(a)") 'weight'
write(file_id,"(2f8.3,1x,i7)") xmin, xmax, 4
write(file_id,"(2f8.3,1x,i7)") ymin, ymax, 4
write(file_id,"(2f8.3,1x,i7)") vxmin, vxmax, 4
write(file_id,"(2f8.3,1x,i7)") vymin, vymax, 4
write(file_id,"(2f8.3,1x,i7)") pmin, pmax, 4
do k = 1, nbpart
  write(file_id,"(6e15.3)")p%dpx(k),p%dpy(:),0.,p%vpx(k),p%vpy(k),p%p(k)
end do
close(file_id)

end subroutine plot_particles_xmdv


end module diagno
