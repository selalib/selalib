module diagno
#include "sll_working_precision.h"
use sll_m_pic_visu
use sll_m_utilities
use sll_m_ascii_io
use sll_m_gnuplot
use sll_m_xdmf
use zone

implicit none

sll_int32, private :: i, j, k

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_part( ele, time, iplot )

type (particle) :: ele
sll_int32 :: iplot
sll_real64 :: time, umod

open(40, file="viry.dat", position="append")
if( iplot == 1 ) rewind(40)
   
umod = sqrt(ele%vit(1,1)**2+ele%vit(1,2)**2)
   
write(40,"(G15.3,1X,5F12.7)") sngl(time)         &
      ,sngl(ele%pos(1,1)), sngl(ele%pos(1,2))    &
      ,sngl(ele%vit(1,1)), sngl(ele%vit(1,2))    &
      ,umod
close(40)

end subroutine plot_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diag_champ_part( ele, time, iplot )

type(particle) :: ele
sll_real64 :: time
sll_int32 :: iplot

if ( nbpart > 0 ) then

   open(17, file="chpart", position="append"  )
   if( iplot == 1 ) rewind(17)
   
   write(17,*) sngl(time), sngl(ele%epx(1)), &
        sngl(ele%epy(1)), sngl(ele%bpz(1))

   close(17)

endif

end subroutine diag_champ_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_phases( ele, iplot, time )

type(particle) :: ele
sll_int32 :: iplot, ipart
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
do ipart=1,nbpart
   speed = sqrt( ele%vit(ipart,1)*ele%vit(ipart,1) + &
        &        ele%vit(ipart,2)*ele%vit(ipart,2) )
   write(14,*) sngl(ele%pos(ipart,1)),sngl(ele%pos(ipart,2))    &
              ,sngl(ele%vit(ipart,1)),sngl(ele%vit(ipart,2)), sngl(speed)
end do
close(14)


end subroutine plot_phases

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Calcul de la fct de distribution en fonction de (vx,vy)

subroutine distribution_v(ele, iplot, time)  

sll_int32 :: i, ipart, iplot
sll_real64, dimension(:,:), allocatable :: df
type(particle) :: ele
sll_real64 :: time, vx, vy, aux, vth=1.0_f64
sll_real64 :: delta_v, vmin, vmax
character(len=4) :: fin
sll_int32, parameter :: nv = 64

call sll_s_int2string(iplot,fin)

allocate(df(nv,nv))

vmin = -6.d0
vmax = 6.d0
delta_v = (vmax-vmin)/nv

df = 0.d0
do ipart=1,nbpart
   do i=1,nv
      do j=1,nv
         if (vmin+(i-1)*delta_v <= ele%vit(ipart,1) .and. &
              & ele%vit(ipart,1) < vmin+i*delta_v  .and. & 
              & vmin+(j-1)*delta_v <= ele%vit(ipart,2) .and. &
              & ele%vit(ipart,2) < vmin+j*delta_v) then
            df(i,j) = df(i,j) + ele%p(ipart)
         endif
      enddo
   enddo
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

if ( iplot .eq. 1 ) then
   open( 37, file = 'df_theo.gnu' )
   rewind(37)
   write(37,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"

   open( 38, file = 'df_theo')

   write(37,*)"splot  'df_theo' w l"
   write(37,*)"pause 1"

   do i=1,nv
      do j=1,nv
         vx = vmin+(i-0.5)*delta_v
         vy = vmin+(j-0.5)*delta_v
         aux = exp(-(vx*vx+vy*vy)/(2*vth*vth))/(2*pi*vth*vth)
         write(38,*) vx, vy, aux*dimx*dimy 
      end do
      write(38,*)
   enddo
   
   close(37)
   close(38)

endif

end subroutine distribution_v

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Calcul de la fct de distribution en fonction de (x,y)

subroutine distribution_x(ele, iplot, time)  

sll_int32 :: i, ipart, iplot
sll_real64, dimension(100,100) :: df
type(particle) :: ele
sll_real64 :: time, x, y, delta_x, delta_y
character(len=4) :: fin

call sll_s_int2string(iplot, fin)

df = 0.d0
delta_x = dimx/100
delta_y = dimy/100
do ipart=1,nbpart
  do i=1,100
    do j=1,100
      if ((i-1)*delta_x <= ele%pos(ipart,1) .and. &
        & ele%pos(ipart,1) < i*delta_x  .and. & 
        & (j-1)*delta_y <= ele%pos(ipart,2) .and. &
        & ele%pos(ipart,2) < j*delta_y) then
        df(i,j) = df(i,j) + ele%p(ipart)
      endif
    enddo
  enddo
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
write(*,*) time, aux

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
  write(file_id,"(4e15.3)")p%pos(k,:),p%vit(k,:)
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

vxmin = minval(p%vit(:,1))
vxmax = maxval(p%vit(:,1))
vymin = minval(p%vit(:,2))
vymax = maxval(p%vit(:,2))
pmin  = minval(p%p(:))
pmax  = maxval(p%p(:))

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
  write(file_id,"(6e15.3)")p%pos(k,:),0.,p%vit(k,:),p%p(k)
end do
close(file_id)

end subroutine plot_particles_xmdv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_particle_density( ele, iplot)  
sll_int32 :: i, ipart, iplot
sll_int32, parameter :: nx = 50, ny = 50
sll_real64, dimension(nx,ny) :: df
type(particle) :: ele
sll_real64 :: delta_x, delta_y, x_min, x_max, y_min, y_max
character(len=4) :: fin
sll_int32 :: error

call sll_s_int2string(iplot, fin)

x_min = 0.0_f64; x_max = dimx
y_min = 0.0_f64; y_max = dimy
df = 0.d0
delta_x = dimx/(nx-1)
delta_y = dimy/(ny-1)
do ipart=1,nbpart
 do i=1,ny
  if ((i-1)*delta_x <= ele%pos(ipart,1) .and. ele%pos(ipart,1) < i*delta_x) then
   do j=1,ny
    if((j-1)*delta_y <= ele%pos(ipart,2) .and. ele%pos(ipart,2) < j*delta_y) then
     df(i,j) = df(i,j) + ele%p(ipart)
    end if
   enddo
  endif
 enddo
enddo

call sll_o_gnuplot_2d(0._f64, dimx, 40, 0._f64, dimy, 40, df, 'density', iplot, error)  
call sll_s_xdmf_corect2d_nodes('df_'//fin,df,"density",x_min,delta_x,y_min,delta_y) 

end subroutine plot_particle_density

end module diagno
