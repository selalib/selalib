!Limit model (1.8)
program test_pic2d
#include "sll_working_precision.h"
#include "sll_memory.h"

use zone
use particules
use diagno
use sll_m_fft
use sll_m_poisson_2d_base
use sll_m_poisson_2d_periodic
use sll_m_constants

implicit none

type(tm_mesh_fields) :: f
type(particle)       :: p

!-----added----
type(sll_t_fft) :: PlnF, PlnB,PlnFy, PlnBy
integer(4)  ,parameter :: npp=200000
real(8)  :: xt(2),vt(2,npp),xu(2,npp),vu(2,npp),lx(128),ly(16)
real(8)  :: temp(npp)
complex(8) ::fn(128,16),ftilde(128),gn(128,16),gtilde(16),tte(128,16)
!----end added--

sll_real64 :: time
sll_real64 :: xmin
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

!---added---
!ep=0.5d0

call sll_s_fft_init_c2c_1d(PlnF,128,fn(:,1),ftilde,sll_p_fft_FORWARD,optimization=sll_p_FFT_MEASURE)
call sll_s_fft_init_c2c_1d(PlnB,128,fn(:,1),ftilde,sll_p_FFT_BACKWARD,optimization=sll_p_FFT_MEASURE)
call sll_s_fft_init_c2c_1d(PlnFy,16,gn(1,:),gtilde,sll_p_fft_FORWARD,optimization=sll_p_FFT_MEASURE)
call sll_s_fft_init_c2c_1d(PlnBy,16,gn(1,:),gtilde,sll_p_FFT_BACKWARD,optimization=sll_p_FFT_MEASURE)

!--end added--

!********************************************************************

!istep = 1

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

m=nx/2
lx=(/ (n, n=0,m-1), (n, n=-m,-1 )/)*2.0d0*sll_p_pi/dimx
m=ny/2
ly=(/ (n, n=0,m-1), (n, n=-m,-1 )/)*2.0d0*sll_p_pi/dimy

call plasma( p ) 
print"('nbpart = ', g15.3)", nbpart
print"('dt = ', g15.3)", dt
call calcul_rho( p, f )

poisson => sll_f_new_poisson_2d_periodic(xmin,xmax,nx,ymin,ymax,ny)
call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)


!----------prepared initial data-----
xu(1,:)=p%dpx+xmin+p%idx*dx
xu(2,:)=p%dpy+ymin+p%idy*dy
vu(1,:)=p%vpx
vu(2,:)=p%vpy
!-------end of preparation-----

!allocate (elec(nstep))
do istep = 1, nstep
    call interpol_eb( f, p )
    !---ode solver 1st order scheme
    do m=1,nbpart
        xt(1)=xu(1,m)+dt*p%epy(m)
        xt(2)=xu(2,m)-dt*p%epx(m)
        call apply_bc()
        p%idx(m) = floor((xt(1)-xmin)/dimx*nx)
        p%dpx(m) = real(xt(1)-xmin- p%idx(m)*dx, f32)
        p%idy(m) = floor((xt(2)-ymin)/dimy*ny)
        p%dpy(m) = real(xt(2)-ymin- p%idy(m)*dy, f32)
        xu(:,m)=xt
    enddo
    tte=f%ex(0:127,0:15)
    do n=1,ny
        call sll_s_fft_exec_c2c_1d(PlnF, tte(:,n), ftilde)
        do m=1,nx
            ftilde(m)=ftilde(m)*sll_p_i1*lx(m)/nx
        enddo
        call sll_s_fft_exec_c2c_1d(PlnB, ftilde,fn(:,n))
    enddo
    tte=f%ey(0:127,0:15)
    do n=1,nx
        call sll_s_fft_exec_c2c_1d(PlnFy, tte(n,:), gtilde)
        do m=1,ny
            gtilde(m)=gtilde(m)*sll_p_i1*ly(m)/ny
        enddo
        call sll_s_fft_exec_c2c_1d(PlnBy, gtilde,   gn(n,:))
    enddo
    f%ex(0:127,0:15)=dreal(fn)
    f%ex(128,:)=f%ex(0,:)
    f%ex(:,10)=f%ex(:,0)
    f%ex(128,16)=f%ex(0,0)
    f%ey(0:127,0:15)=dreal(gn)
    f%ey(128,:)=f%ey(0,:)
    f%ey(:,16)=f%ey(:,0)
    f%ey(128,16)=f%ey(0,0)
    call interpol_eb( f, p )
    temp=p%epx+p%epy
    do n=1,nbpart
        vt(1,n)=(vu(1,n)-dt*temp(n)/2.0d0*vu(2,n))/(1.0d0+(dt*temp(n)/2.0d0)**2)!vu(1,:)-dt*temp/2.0d0*vu(2,:)!
        vt(2,n)=(vu(2,n)+dt*temp(n)/2.0d0*vu(1,n))/(1.0d0+(dt*temp(n)/2.0d0)**2)!vu(2,:)+dt*temp/2.0d0*vu(1,:)!
    enddo
    vu=vt
    call calcul_rho( p, f )
    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
enddo
call sll_s_fft_free(PlnF)
call sll_s_fft_free(PlnB)
call sll_s_fft_free(PlnFy)
call sll_s_fft_free(PlnBy)
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

!deallocate (elec)

contains


subroutine apply_bc()
do while ( xt(1) > xmax )
xt(1) = xt(1) - dimx
enddo
do while ( xt(1) < xmin )
xt(1)= xt(1) + dimx
enddo
do while ( xt(2) > ymax )
xt(2)  = xt(2)  - dimy
enddo
do while ( xt(2)  < ymin )
xt(2) = xt(2)  + dimy
enddo
end subroutine apply_bc

end program test_pic2d
