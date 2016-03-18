! rk4 solver
program test_pic2d
#include "sll_working_precision.h"
#include "sll_memory.h"

use zone
use particules
use diagno

use sll_m_poisson_2d_base
use sll_m_poisson_2d_periodic

implicit none

type(tm_mesh_fields) :: f
type(particle)       :: p

sll_real64 :: time,xt(2,120000),vt(2,120000),ep,xtemp(2,120000)
sll_real64 :: xmin,Et1(4,120000),Et2(4,120000),Et3(4,120000),Et4(4,120000)
sll_real64 :: xmax
sll_real64 :: ymin
sll_real64 :: ymax

sll_int32  :: istep
sll_int32  :: iplot
sll_int32  :: iargc
sll_int32  :: n,m
sll_int32  :: i
sll_int32  :: j
sll_int32  :: error

sll_real64 :: aux1, aux2
sll_real64 :: t0, t1

character(len=272) :: argv
class(sll_c_poisson_2d_base), pointer :: poisson

n = iargc()
if (n == 0) stop 'Usage: ./bin/test_pic2d fichier-de-donnees.nml'
do i = 1, n
  call getarg( i, argv)
  write(*,'(i2, 1x, a)') i, argv
end do

call readin( trim(argv) )

SLL_CLEAR_ALLOCATE(f%ex(0:nx,0:ny), error)
SLL_CLEAR_ALLOCATE(f%ey(0:nx,0:ny), error)
SLL_CLEAR_ALLOCATE(f%bz(0:nx,0:ny), error)
SLL_CLEAR_ALLOCATE(f%r0(0:nx,0:ny), error) 

time  = 0.d0
iplot = 0
ep=0.5d0
!if( nstep > nstepmax ) nstep = nstepmax

!********************************************************************
print"('dt = ', g15.3)", dt
print"('eps= ', g15.3)", ep
istep = 1

do i=0,nx
  aux1 = alpha/kx * sin(kx*i*dx)
  aux2 = alpha * cos(kx*i*dx)
  do j=0,ny
    f%ex(i,j) = aux1
    f%r0(i,j) = aux2
  enddo
enddo
      
xmin = 0.0_f64; xmax = dimx
ymin = 0.0_f64; ymax = dimy

call plasma( p ) 

call calcul_rho( p, f )

!gnuplot -p rho.gnu (to plot the initial rho)
!call sll_o_gnuplot_2d(xmin, xmax, nx+1, &
!                      ymin, ymax, ny+1, &
!                      f%r0, 'rho', 1, error)

poisson => sll_f_new_poisson_2d_periodic( xmin, xmax, nx, &
                       ymin, ymax, ny) 

call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)

!gnuplot -p ex.gnu (to plot the initial ex field)
!call sll_o_gnuplot_2d(xmin, xmax, nx+1, &
!                      ymin, ymax, ny+1, &
!                      f%ex, 'ex', 1, error)

!call cpu_time(t0)
do istep = 1, nstep
  call interpol_eb( f, p )
    xtemp(1,:)=p%dpx+xmin+p%idx*dx
    xtemp(2,:)=p%dpy+ymin+p%idy*dy
Et1(1,:)=p%vpx/ep
Et1(2,:)=p%vpy/ep
Et1(3,:)=p%epx/ep+p%vpy/ep**2
Et1(4,:)=p%epy/ep-p%vpx/ep**2
    do m=1,nbpart
        vt(1,m)=p%vpx(m)+dt/2.0d0*Et1(3,m)
        vt(2,m)=p%vpy(m)+dt/2.0d0*Et1(4,m)
        xt(1,m)=xtemp(1,m)+dt/2.0d0*Et1(1,m)
        xt(2,m)=xtemp(2,m)+dt/2.0d0*Et1(2,m)
        call apply_bc()
        p%idx(m) = floor((xt(1,m)-xmin)/dimx*nx)
        p%dpx(m) = real(xt(1,m) -xmin- p%idx(m)*dx, f32)
        p%idy(m) = floor((xt(2,m)-ymin)/dimy*ny)
        p%dpy(m) = real(xt(2,m) -ymin- p%idy(m)*dy, f32)
    enddo
    call interpol_eb( f, p )
Et2(1,:)=vt(1,:)/ep
Et2(2,:)=vt(2,:)/ep
Et2(3,:)=p%epx/ep+vt(2,:)/ep**2
Et2(4,:)=p%epy/ep-vt(1,:)/ep**2
    do m=1,nbpart
        vt(1,m)=p%vpx(m)+dt/2.0d0*Et2(3,m)
        vt(2,m)=p%vpy(m)+dt/2.0d0*Et2(4,m)
        xt(1,m)=xtemp(1,m)+dt/2.0d0*Et2(1,m)
        xt(2,m)=xtemp(2,m)+dt/2.0d0*Et2(2,m)
        call apply_bc()
        p%idx(m) = floor((xt(1,m)-xmin)/dimx*nx)
        p%dpx(m) = real(xt(1,m) -xmin- p%idx(m)*dx, f32)
        p%idy(m) = floor((xt(2,m)-ymin)/dimy*ny)
        p%dpy(m) = real(xt(2,m) -ymin- p%idy(m)*dy, f32)
    enddo
    call interpol_eb( f, p )
Et3(1,:)=vt(1,:)/ep
Et3(2,:)=vt(2,:)/ep
Et3(3,:)=p%epx/ep+vt(2,:)/ep**2
Et3(4,:)=p%epy/ep-vt(1,:)/ep**2
    do m=1,nbpart
        vt(1,m)=p%vpx(m)+dt*Et3(3,m)
        vt(2,m)=p%vpy(m)+dt*Et3(4,m)
        xt(1,m)=xtemp(1,m)+dt*Et3(1,m)
        xt(2,m)=xtemp(2,m)+dt*Et3(2,m)
        call apply_bc()
        p%idx(m) = floor((xt(1,m)-xmin)/dimx*nx)
        p%dpx(m) = real(xt(1,m) -xmin- p%idx(m)*dx, f32)
        p%idy(m) = floor((xt(2,m)-ymin)/dimy*ny)
        p%dpy(m) = real(xt(2,m) -ymin- p%idy(m)*dy, f32)
    enddo
    call interpol_eb( f, p )
Et4(1,:)=vt(1,:)/ep
Et4(2,:)=vt(2,:)/ep
Et4(3,:)=p%epx/ep+vt(2,:)/ep**2
Et4(4,:)=p%epy/ep-vt(1,:)/ep**2
    do m=1,nbpart
        vt(1,m)=p%vpx(m)+dt/6.0d0*(Et1(3,m)+2.0d0*Et2(3,m)+2.0d0*Et3(3,m)+Et4(3,m))
        vt(2,m)=p%vpy(m)+dt/6.0d0*(Et1(4,m)+2.0d0*Et2(4,m)+2.0d0*Et3(4,m)+Et4(4,m))
        xt(1,m)=xtemp(1,m)+dt/6.0d0*(Et1(1,m)+2.0d0*Et2(1,m)+2.0d0*Et3(1,m)+Et4(1,m))
        xt(2,m)=xtemp(2,m)+dt/6.0d0*(Et1(2,m)+2.0d0*Et2(2,m)+2.0d0*Et3(2,m)+Et4(2,m))
        call apply_bc()
        p%idx(m) = floor((xt(1,m)-xmin)/dimx*nx)
        p%dpx(m) = real(xt(1,m) -xmin- p%idx(m)*dx, f32)
        p%idy(m) = floor((xt(2,m)-ymin)/dimy*ny)
        p%dpy(m) = real(xt(2,m) -ymin- p%idy(m)*dy, f32)
    enddo
    p%vpx=vt(1,:)
    p%vpy=vt(2,:)
  call calcul_rho( p, f )
  call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
  time = time + dt
  !call modeE( f, istep, time )
  !write(*,"('istep = ', i6, ' time = ',g15.3)", advance='no') istep, time
end do

!call cpu_time(t1)
!print"('CPU time = ', g15.3)", t1-t0
print*,'PASSED'

open(unit=851,file='fh1.dat')
do i=1,nx
do j=1,ny
write(851,*)f%r0(i,j)!fh_fsl(i,j)!  , !eta2feet(i,j),eta1feet(i,j),
enddo
enddo
close(851)

contains

subroutine apply_bc()
do while ( xt(1,m) > xmax )
xt(1,m) = xt(1,m) - dimx
enddo
do while ( xt(1,m) < xmin )
xt(1,m)= xt(1,m) + dimx
enddo
do while ( xt(2,m) > ymax )
xt(2,m)  = xt(2,m)  - dimy
enddo
do while ( xt(2,m)  < ymin )
xt(2,m) = xt(2,m)  + dimy
enddo
end subroutine apply_bc

end program test_pic2d
