!NKS-two-scale solver for 4d VP
!2nd order EWI scheme
!two-scale E
program test_pic2d
#include "sll_working_precision.h"
#include "sll_memory.h"

use zone
use particules
use sll_m_fft
use sll_m_poisson_2d_base
use sll_m_poisson_2d_periodic
use sll_m_constants

implicit none

type(tm_mesh_fields) :: f
type(particle)       :: p

type(sll_t_fft) :: PlnF, PlnB
integer(4)  ,parameter :: Ntau=32,npp=204800
real(8)  :: xxt(2)
real(8)  :: ep,dtau,tau(0:Ntau-1),ltau(0:Ntau-1),auxpx(2,npp)
complex(8) :: pl(0:ntau-1),ql(0:ntau-1),gp(2,0:Ntau-1,npp),gm(2,0:Ntau-1,npp),wp(2,npp),wm(2,npp)
complex(8) ::  xtemp1(2,0:ntau-1,npp),xtemp2(2,0:ntau-1,npp),Et(2,0:ntau-1,npp)
complex(8) :: temp(2,0:Ntau-1),temptilde(2,0:Ntau-1),up(2,0:ntau-1,npp),um(2,0:ntau-1,npp)
complex(8) :: vp(2,npp),vm(2,npp),z(2),fex(0:64,0:32,0:ntau-1),fey(0:64,0:32,0:ntau-1)
complex(8) :: up0(2,0:ntau-1,npp),um0(2,0:ntau-1,npp)
!real(8),  allocatable :: energy(:)

sll_real64 :: time
sll_real64 :: xmin
sll_real64 :: xmax
sll_real64 :: ymin
sll_real64 :: ymax
sll_int32  :: istep = 1
sll_int32  :: iplot
sll_int32  :: iargc
sll_int32  :: n,m
sll_int32  :: i
sll_int32  :: j
sll_int32  :: error

sll_real64 :: aux1, aux2
sll_real64 :: t1

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

ep=0.1d0/1.0d0
dtau=2.0d0*sll_p_pi/ntau
call sll_s_fft_init_c2c_1d(PlnF,Ntau,temp(1,:),temptilde(1,:),sll_p_fft_FORWARD)
call sll_s_fft_init_c2c_1d(PlnB,Ntau,temp(1,:),temptilde(1,:),sll_p_FFT_BACKWARD)
m=ntau/2
ltau=(/ (n, n=0,m-1), (n, n=-m,-1 )/)
pl(0)=dt
ql(0)=dt**2/2.0d0
do i=1,Ntau-1
    pl(i) =2.0d0*ep**2*sll_p_i1*(exp(-sll_p_i1*ltau(i)*dt/2.0d0/ep**2)-1.0d0)/ltau(i)
    ql(i)=2.0d0*ep**2*(2.0d0*ep**2*(1.0d0-exp(-sll_p_i1*ltau(i)*dt/2.0d0/ep**2))-sll_p_i1*ltau(i)*dt)/ltau(i)**2
enddo
do i=0,Ntau-1
  tau(i) =i*dtau
enddo

do i=0,nx
  aux1 = alpha/kx * sin(kx*i*dx)
  aux2 = alpha * cos(kx*i*dx)
  do j=0,ny
    f%ex(i,j) = aux1
    f%r0(i,j) = aux2+dsin(ky*j*dy)
    f%ey(i,j) = -dcos(ky*j*dy)/ky!0.0d0!
  enddo
enddo
      
xmin = 0.0_f64; xmax = dimx
ymin = 0.0_f64; ymax = dimy

call plasma( p )
print"('ep     = ', g15.3)", ep
print"('nbpart = ', g15.3)", nbpart
print"('dt     = ', g15.3)", dt
print"('ntau   = ', g15.3)", ntau
print"('nx     = ', g15.3)", nx
print"('ny     = ', g15.3)", ny
print"('kx     = ', g15.3)", kx
print"('ky     = ', g15.3)", ky
print"('xmin, xmax  = ', 2g15.3)", xmin, xmax
print"('ymin, ymax  = ', 2g15.3)", ymin, ymax
print"('alpha  = ', g15.3)", alpha

poisson => sll_f_new_poisson_2d_periodic(xmin,xmax,nx,ymin,ymax,ny)
!call calcul_rho_m6( p, f )
!call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
!call calcul_energy(p,f,energy(0))

auxpx(1,:)=(p%dpx+p%idx)*dx
auxpx(2,:)=(p%dpy+p%idy)*dy
wp(1,:)=2.0d0*(auxpx(1,:)+ep*p%vpy)
wp(2,:)=2.0d0*(auxpx(2,:)-ep*p%vpx)
wm(1,:)=-2.0d0*ep*p%vpy
wm(2,:)=2.0d0*ep*p%vpx
do n=0,ntau-1
    xtemp1(1,n,:)=dcos(tau(n))*wp(1,:)-dsin(tau(n))*wp(2,:)
    xtemp1(2,n,:)=dsin(tau(n))*wp(1,:)+dcos(tau(n))*wp(2,:)
    xtemp2(1,n,:)=dcos(tau(n))*wm(1,:)+dsin(tau(n))*wm(2,:)
    xtemp2(2,n,:)=-dsin(tau(n))*wm(1,:)+dcos(tau(n))*wm(2,:)
enddo
xtemp1=(xtemp1+xtemp2)/2.0d0
do n=0,ntau-1
    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)
    xtemp2(2,n,:)=-dsin(tau(n))*xtemp1(1,n,:)+dcos(tau(n))*xtemp1(2,n,:)
    do m=1,nbpart
        xxt=dreal(xtemp2(:,n,m))
        call apply_bc()
        p%idx(m) = floor(xxt(1)/dimx*nx)
        p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
        p%idy(m) = floor(xxt(2)/dimy*ny)
        p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    enddo
    call interpol_eb_m6( f, p )
    Et(1,n,:)=p%epx !g(0,tau,w(0))
    Et(2,n,:)=p%epy
enddo
do m=1,npp
    xtemp1(1,:,m)=2.0d0*Et(2,:,m)
    xtemp1(2,:,m)=-2.0d0*Et(1,:,m)!g_+
    call sll_s_fft_exec_c2c_1d(PlnF, xtemp1(1,:,m), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, xtemp1(2,:,m), temptilde(2,:))
    do n=1,Ntau-1
        temptilde(:,n)=-sll_p_i1*temptilde(:,n)/ltau(n)/Ntau
    enddo
    temptilde(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(1,:), temp(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(2,:), temp(2,:))!AF+
    up(1,:,m)=wp(1,m)+2.0d0*ep**2*(temp(1,:)-temp(1,0))
    up(2,:,m)=wp(2,m)+2.0d0*ep**2*(temp(2,:)-temp(2,0))!1st ini data of U_+
    !---
    do n=0,ntau-1
        xtemp2(1,n,m)=-2.0d0*(dsin(2.0d0*tau(n))*Et(1,n,m)+dcos(2.0d0*tau(n))*Et(2,n,m))
        xtemp2(2,n,m)=2.0d0*(-dsin(2.0d0*tau(n))*Et(2,n,m)+dcos(2.0d0*tau(n))*Et(1,n,m))!g_-
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, xtemp2(1,:,m), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, xtemp2(2,:,m), temptilde(2,:))
    do n=1,Ntau-1
        temptilde(:,n)=-sll_p_i1*temptilde(:,n)/ltau(n)/Ntau
    enddo
    temptilde(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(1,:), temp(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(2,:), temp(2,:))!AF-
    um(1,:,m)=wm(1,m)+2.0d0*ep**2*(temp(1,:)-temp(1,0))
    um(2,:,m)=wm(2,m)+2.0d0*ep**2*(temp(2,:)-temp(2,0))!1st ini data of U_-
enddo
!--corrected more initial data
do n=0,ntau-1
    xtemp1(1,n,:)=dcos(tau(n))*up(1,n,:)-dsin(tau(n))*up(2,n,:)!v_+ 2scaled
    xtemp1(2,n,:)=dcos(tau(n))*up(2,n,:)+dsin(tau(n))*up(1,n,:)
    xtemp2(1,n,:)=dcos(tau(n))*um(1,n,:)+dsin(tau(n))*um(2,n,:)!v_- 2scaled
    xtemp2(2,n,:)=dcos(tau(n))*um(2,n,:)-dsin(tau(n))*um(1,n,:)
    xtemp1(:,n,:)=(xtemp1(:,n,:)+xtemp2(:,n,:))/2.0d0!z 2scaled
    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)!x 2scaled
    xtemp2(2,n,:)=dcos(tau(n))*xtemp1(2,n,:)-dsin(tau(n))*xtemp1(1,n,:)
    do m=1,nbpart
        xxt=dreal(xtemp2(:,n,m))
        call apply_bc()
        p%idx(m) = floor(xxt(1)/dimx*nx)
        p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
        p%idy(m) = floor(xxt(2)/dimy*ny)
        p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    enddo
    call calcul_rho_m6( p, f )
    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
    fex(:,:,n)=f%ex
    fey(:,:,n)=f%ey!E_1st(0,x)
enddo
do n=0,ntau-1
    xtemp1(1,n,:)=dcos(tau(n))*up(1,n,:)-dsin(tau(n))*up(2,n,:)
    xtemp1(2,n,:)=dsin(tau(n))*up(1,n,:)+dcos(tau(n))*up(2,n,:)
    xtemp2(1,n,:)=dcos(tau(n))*um(1,n,:)+dsin(tau(n))*um(2,n,:)
    xtemp2(2,n,:)=-dsin(tau(n))*um(1,n,:)+dcos(tau(n))*um(2,n,:)
enddo
xtemp1=(xtemp1+xtemp2)/2.0d0
do n=0,ntau-1
    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)
    xtemp2(2,n,:)=-dsin(tau(n))*xtemp1(1,n,:)+dcos(tau(n))*xtemp1(2,n,:)
    do m=1,nbpart
        xxt=dreal(xtemp2(:,n,m))
        call apply_bc()
        p%idx(m) = floor(xxt(1)/dimx*nx)
        p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
        p%idy(m) = floor(xxt(2)/dimy*ny)
        p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    enddo
    f%ex=fex(:,:,n)
    f%ey=fey(:,:,n)
    call interpol_eb_m6( f, p )
    Et(1,n,:)=p%epx !g_1st(0,tau,U_1st(0))
    Et(2,n,:)=p%epy
enddo
do m=1,npp
    xtemp1(1,:,m)=2.0d0*Et(2,:,m)
    xtemp1(2,:,m)=-2.0d0*Et(1,:,m)!g_+
    call sll_s_fft_exec_c2c_1d(PlnF, xtemp1(1,:,m), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, xtemp1(2,:,m), temptilde(2,:))
    up0(:,0,m)=temptilde(:,0)/ntau!Pi g_+
    do n=1,Ntau-1
        temptilde(:,n)=-sll_p_i1*temptilde(:,n)/ltau(n)/Ntau
    enddo
    temptilde(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(1,:), temp(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(2,:), temp(2,:))!AF+
    up(1,:,m)=wp(1,m)+2.0d0*ep**2*(temp(1,:)-temp(1,0))
    up(2,:,m)=wp(2,m)+2.0d0*ep**2*(temp(2,:)-temp(2,0))!3rd ini data of U_+
    !---
    do n=0,ntau-1
        xtemp2(1,n,m)=-2.0d0*(dsin(2.0d0*tau(n))*Et(1,n,m)+dcos(2.0d0*tau(n))*Et(2,n,m))
        xtemp2(2,n,m)=2.0d0*(-dsin(2.0d0*tau(n))*Et(2,n,m)+dcos(2.0d0*tau(n))*Et(1,n,m))!g_-
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, xtemp2(1,:,m), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, xtemp2(2,:,m), temptilde(2,:))
    um0(:,0,m)=temptilde(:,0)/ntau!Pi g_-
    do n=1,Ntau-1
        temptilde(:,n)=-sll_p_i1*temptilde(:,n)/ltau(n)/Ntau
    enddo
    temptilde(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(1,:), temp(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(2,:), temp(2,:))!AF-
    um(1,:,m)=wm(1,m)+2.0d0*ep**2*(temp(1,:)-temp(1,0))
    um(2,:,m)=wm(2,m)+2.0d0*ep**2*(temp(2,:)-temp(2,0))!3rd ini data of U_-
enddo

do n=0,ntau-1
    xtemp1(1,n,:)=dcos(tau(n))*up(1,n,:)-dsin(tau(n))*up(2,n,:)!v_+ 2scaled
    xtemp1(2,n,:)=dcos(tau(n))*up(2,n,:)+dsin(tau(n))*up(1,n,:)
    xtemp2(1,n,:)=dcos(tau(n))*um(1,n,:)+dsin(tau(n))*um(2,n,:)!v_- 2scaled
    xtemp2(2,n,:)=dcos(tau(n))*um(2,n,:)-dsin(tau(n))*um(1,n,:)
    xtemp1(:,n,:)=(xtemp1(:,n,:)+xtemp2(:,n,:))/2.0d0!z 2scaled
    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)!x 2scaled
    xtemp2(2,n,:)=dcos(tau(n))*xtemp1(2,n,:)-dsin(tau(n))*xtemp1(1,n,:)
    do m=1,nbpart
        xxt=dreal(xtemp2(:,n,m))
        call apply_bc()
        p%idx(m) = floor(xxt(1)/dimx*nx)
        p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
        p%idy(m) = floor(xxt(2)/dimy*ny)
        p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    enddo
    call calcul_rho_m6( p, f )
    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
    fex(:,:,n)=f%ex
    fey(:,:,n)=f%ey!E_4(0,x)
enddo

!--time iteration---
time=dt
do n=0,ntau-1
    xtemp1(1,n,:)=dcos(tau(n))*up(1,n,:)-dsin(tau(n))*up(2,n,:)
    xtemp1(2,n,:)=dsin(tau(n))*up(1,n,:)+dcos(tau(n))*up(2,n,:)
    xtemp2(1,n,:)=dcos(tau(n))*um(1,n,:)+dsin(tau(n))*um(2,n,:)
    xtemp2(2,n,:)=-dsin(tau(n))*um(1,n,:)+dcos(tau(n))*um(2,n,:)
enddo
xtemp1=(xtemp1+xtemp2)/2.0d0
do n=0,ntau-1
    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)
    xtemp2(2,n,:)=-dsin(tau(n))*xtemp1(1,n,:)+dcos(tau(n))*xtemp1(2,n,:)
    do m=1,nbpart
        xxt=dreal(xtemp2(:,n,m))
        call apply_bc()
        p%idx(m) = floor(xxt(1)/dimx*nx)
        p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
        p%idy(m) = floor(xxt(2)/dimy*ny)
        p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    enddo
    f%ex=fex(:,:,n)
    f%ey=fey(:,:,n)
    call interpol_eb_m6( f, p )
    Et(1,n,:)=p%epx !g_3rd(0,tau,U_3rd(0))
    Et(2,n,:)=p%epy
enddo
do m=1,npp
    temp(1,:)=2.0d0*Et(2,:,m)
    temp(2,:)=-2.0d0*Et(1,:,m)
    call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), temptilde(2,:))
    xtemp1(:,:,m)=temptilde/ntau!g_+tilde(t=0)
    !---
    do n=0,ntau-1
        temp(1,n)=-2.0d0*(dsin(2.0d0*tau(n))*Et(1,n,m)+dcos(2.0d0*tau(n))*Et(2,n,m))
        temp(2,n)=2.0d0*(-dsin(2.0d0*tau(n))*Et(2,n,m)+dcos(2.0d0*tau(n))*Et(1,n,m))
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), temptilde(2,:))
    xtemp2(:,:,m)=temptilde/ntau!g_-tilde(t=0)
enddo
do m=1,npp
    call sll_s_fft_exec_c2c_1d(PlnF, up(1,:,m), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, up(2,:,m), temptilde(2,:))
    do n=0,ntau-1
        temp(:,n)=exp(-sll_p_i1*ltau(n)*dt/2.0d0/ep**2)*temptilde(:,n)/ntau+pl(n)*xtemp1(:,n,m)!utilde_+^1,predict
    enddo
    call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), up0(1,:,m))!u_+(t1),predict
    call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), up0(2,:,m))
    call sll_s_fft_exec_c2c_1d(PlnF, um(1,:,m), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, um(2,:,m), temptilde(2,:))
    do n=0,ntau-1
        temp(:,n)=exp(-sll_p_i1*ltau(n)*dt/2.0d0/ep**2)*temptilde(:,n)/ntau+pl(n)*xtemp2(:,n,m)!utilde_-^1,predict
    enddo
    call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), um0(1,:,m))!u_-(t1),predict
    call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), um0(2,:,m))
enddo
gp=xtemp1
gm=xtemp2
do n=0,ntau-1
    xtemp1(1,n,:)=dcos(tau(n))*up0(1,n,:)-dsin(tau(n))*up0(2,n,:)!v_+ 2scaled
    xtemp1(2,n,:)=dcos(tau(n))*up0(2,n,:)+dsin(tau(n))*up0(1,n,:)
    xtemp2(1,n,:)=dcos(tau(n))*um0(1,n,:)+dsin(tau(n))*um0(2,n,:)!v_- 2scaled
    xtemp2(2,n,:)=dcos(tau(n))*um0(2,n,:)-dsin(tau(n))*um0(1,n,:)
    xtemp1(:,n,:)=(xtemp1(:,n,:)+xtemp2(:,n,:))/2.0d0!z 2scaled
    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)!x 2scaled
    xtemp2(2,n,:)=dcos(tau(n))*xtemp1(2,n,:)-dsin(tau(n))*xtemp1(1,n,:)
    do m=1,nbpart
        xxt=dreal(xtemp2(:,n,m))
        call apply_bc()
        p%idx(m) = floor(xxt(1)/dimx*nx)
        p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
        p%idy(m) = floor(xxt(2)/dimy*ny)
        p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    enddo
    call calcul_rho_m6( p, f )
    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
    fex(:,:,n)=f%ex
    fey(:,:,n)=f%ey!prediction
enddo
!--correction--
do n=0,ntau-1
    xtemp1(1,n,:)=dcos(tau(n))*up0(1,n,:)-dsin(tau(n))*up0(2,n,:)
    xtemp1(2,n,:)=dsin(tau(n))*up0(1,n,:)+dcos(tau(n))*up0(2,n,:)
    xtemp2(1,n,:)=dcos(tau(n))*um0(1,n,:)+dsin(tau(n))*um0(2,n,:)
    xtemp2(2,n,:)=-dsin(tau(n))*um0(1,n,:)+dcos(tau(n))*um0(2,n,:)
enddo
xtemp1=(xtemp1+xtemp2)/2.0d0
do n=0,ntau-1
    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)
    xtemp2(2,n,:)=-dsin(tau(n))*xtemp1(1,n,:)+dcos(tau(n))*xtemp1(2,n,:)
    do m=1,nbpart
        xxt=dreal(xtemp2(:,n,m))
        call apply_bc()
        p%idx(m) = floor(xxt(1)/dimx*nx)
        p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
        p%idy(m) = floor(xxt(2)/dimy*ny)
        p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    enddo
    f%ex=fex(:,:,n)
    f%ey=fey(:,:,n)
    call interpol_eb_m6( f, p )
    Et(1,n,:)=p%epx !g(t1,tau,U(t1))
    Et(2,n,:)=p%epy
enddo
do m=1,npp
    temp(1,:)=2.0d0*Et(2,:,m)
    temp(2,:)=-2.0d0*Et(1,:,m)
    call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), temptilde(2,:))
    xtemp1(:,:,m)=temptilde/ntau!g_+tilde(t1) predict
    !---
    do n=0,ntau-1
        temp(1,n)=-2.0d0*(dsin(2.0d0*tau(n))*Et(1,n,m)+dcos(2.0d0*tau(n))*Et(2,n,m))
        temp(2,n)=2.0d0*(-dsin(2.0d0*tau(n))*Et(2,n,m)+dcos(2.0d0*tau(n))*Et(1,n,m))
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), temptilde(2,:))
    xtemp2(:,:,m)=temptilde/ntau!g_-tilde(t1) predict
enddo
do m=1,npp
    call sll_s_fft_exec_c2c_1d(PlnF, up(1,:,m), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, up(2,:,m), temptilde(2,:))
    do n=0,ntau-1
        temp(:,n)=exp(-sll_p_i1*ltau(n)*dt/2.0d0/ep**2)*temptilde(:,n)/ntau+pl(n)*xtemp1(:,n,m) &
                    +ql(n)*(xtemp1(:,n,m)-gp(:,n,m))/dt
    enddo
    call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), up(1,:,m))!u_+(t1)
    call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), up(2,:,m))
    call sll_s_fft_exec_c2c_1d(PlnF, um(1,:,m), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, um(2,:,m), temptilde(2,:))
    do n=0,ntau-1
        temp(:,n)=exp(-sll_p_i1*ltau(n)*dt/2.0d0/ep**2)*temptilde(:,n)/ntau+pl(n)*xtemp2(:,n,m) &
                    +ql(n)*(xtemp2(:,n,m)-gm(:,n,m))/dt
    enddo
    call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), um(1,:,m))!u_-(t1)
    call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), um(2,:,m))
enddo
do n=0,ntau-1
    xtemp1(1,n,:)=dcos(tau(n))*up(1,n,:)-dsin(tau(n))*up(2,n,:)!v_+ 2scaled
    xtemp1(2,n,:)=dcos(tau(n))*up(2,n,:)+dsin(tau(n))*up(1,n,:)
    xtemp2(1,n,:)=dcos(tau(n))*um(1,n,:)+dsin(tau(n))*um(2,n,:)!v_- 2scaled
    xtemp2(2,n,:)=dcos(tau(n))*um(2,n,:)-dsin(tau(n))*um(1,n,:)
    xtemp1(:,n,:)=(xtemp1(:,n,:)+xtemp2(:,n,:))/2.0d0!z 2scaled
    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)!x 2scaled
    xtemp2(2,n,:)=dcos(tau(n))*xtemp1(2,n,:)-dsin(tau(n))*xtemp1(1,n,:)
    do m=1,nbpart
        xxt=dreal(xtemp2(:,n,m))
        call apply_bc()
        p%idx(m) = floor(xxt(1)/dimx*nx)
        p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
        p%idy(m) = floor(xxt(2)/dimy*ny)
        p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    enddo
    call calcul_rho_m6( p, f )
    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
    fex(:,:,n)=f%ex
    fey(:,:,n)=f%ey
enddo
!--for energy
call energyuse()
!call calcul_energy(p,f,energy(1))
!--end for energy

do istep = 2, nstep
    do n=0,ntau-1
        xtemp1(1,n,:)=dcos(tau(n))*up(1,n,:)-dsin(tau(n))*up(2,n,:)
        xtemp1(2,n,:)=dsin(tau(n))*up(1,n,:)+dcos(tau(n))*up(2,n,:)
        xtemp2(1,n,:)=dcos(tau(n))*um(1,n,:)+dsin(tau(n))*um(2,n,:)
        xtemp2(2,n,:)=-dsin(tau(n))*um(1,n,:)+dcos(tau(n))*um(2,n,:)
    enddo
    xtemp1=(xtemp1+xtemp2)/2.0d0
    do n=0,ntau-1
        xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)
        xtemp2(2,n,:)=-dsin(tau(n))*xtemp1(1,n,:)+dcos(tau(n))*xtemp1(2,n,:)
        do m=1,nbpart
            xxt=dreal(xtemp2(:,n,m))
            call apply_bc()
            p%idx(m) = floor(xxt(1)/dimx*nx)
            p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
            p%idy(m) = floor(xxt(2)/dimy*ny)
            p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
        enddo
        f%ex=fex(:,:,n)
        f%ey=fey(:,:,n)
        call interpol_eb_m6( f, p )
        Et(1,n,:)=p%epx !g(tn)
        Et(2,n,:)=p%epy
    enddo
    do m=1,npp
        temp(1,:)=2.0d0*Et(2,:,m)
        temp(2,:)=-2.0d0*Et(1,:,m)
        call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), temptilde(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), temptilde(2,:))
        xtemp1(:,:,m)=temptilde/ntau!g_+tilde(tn)
        !---
        do n=0,ntau-1
            temp(1,n)=-2.0d0*(dsin(2.0d0*tau(n))*Et(1,n,m)+dcos(2.0d0*tau(n))*Et(2,n,m))
            temp(2,n)=2.0d0*(-dsin(2.0d0*tau(n))*Et(2,n,m)+dcos(2.0d0*tau(n))*Et(1,n,m))
        enddo
        call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), temptilde(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), temptilde(2,:))
        xtemp2(:,:,m)=temptilde/ntau!g_-tilde(tn)
    enddo
    do m=1,npp
        call sll_s_fft_exec_c2c_1d(PlnF, up(1,:,m), temptilde(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, up(2,:,m), temptilde(2,:))
        do n=0,ntau-1
            temp(:,n)=exp(-sll_p_i1*ltau(n)*dt/2.0d0/ep**2)*temptilde(:,n)/ntau+pl(n)*xtemp1(:,n,m) &
                        +ql(n)*(xtemp1(:,n,m)-gp(:,n,m))/dt
        enddo
        call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), up(1,:,m))!u_+(tn)
        call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), up(2,:,m))
        call sll_s_fft_exec_c2c_1d(PlnF, um(1,:,m), temptilde(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, um(2,:,m), temptilde(2,:))
        do n=0,ntau-1
            temp(:,n)=exp(-sll_p_i1*ltau(n)*dt/2.0d0/ep**2)*temptilde(:,n)/ntau+pl(n)*xtemp2(:,n,m) &
                        +ql(n)*(xtemp2(:,n,m)-gm(:,n,m))/dt
        enddo
        call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), um(1,:,m))!u_-(tn)
        call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), um(2,:,m))
    enddo
    gp=xtemp1
    gm=xtemp2
    !--updata E--
    time=dt*istep
    do n=0,ntau-1
        xtemp1(1,n,:)=dcos(tau(n))*up(1,n,:)-dsin(tau(n))*up(2,n,:)!v_+ 2scaled
        xtemp1(2,n,:)=dcos(tau(n))*up(2,n,:)+dsin(tau(n))*up(1,n,:)
        xtemp2(1,n,:)=dcos(tau(n))*um(1,n,:)+dsin(tau(n))*um(2,n,:)!v_- 2scaled
        xtemp2(2,n,:)=dcos(tau(n))*um(2,n,:)-dsin(tau(n))*um(1,n,:)
        xtemp1(:,n,:)=(xtemp1(:,n,:)+xtemp2(:,n,:))/2.0d0!z 2scaled
        xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)!x 2scaled
        xtemp2(2,n,:)=dcos(tau(n))*xtemp1(2,n,:)-dsin(tau(n))*xtemp1(1,n,:)
        do m=1,nbpart
            xxt=dreal(xtemp2(:,n,m))
            call apply_bc()
            p%idx(m) = floor(xxt(1)/dimx*nx)
            p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
            p%idy(m) = floor(xxt(2)/dimy*ny)
            p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
        enddo
        call calcul_rho_m6( p, f )
        call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
        fex(:,:,n)=f%ex
        fey(:,:,n)=f%ey
    enddo
    !--for energy
    call energyuse()
!    call calcul_energy(p,f,energy(istep))
    !--end for energy
enddo

do m=1,npp
    call sll_s_fft_exec_c2c_1d(PlnF, up(1,:,m),temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, up(2,:,m),temptilde(2,:))
    wp(:,m)=0.0d0
    do n=0,ntau-1
        wp(:,m)=wp(:,m)+temptilde(:,n)/Ntau*cdexp(sll_p_i1*ltau(n)*time/2.0d0/ep**2)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, um(1,:,m),temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, um(2,:,m),temptilde(2,:))
    wm(:,m)=0.0d0
    do n=0,ntau-1
        wm(:,m)=wm(:,m)+temptilde(:,n)/Ntau*cdexp(sll_p_i1*ltau(n)*time/2.0d0/ep**2)
    enddo
    vp(1,m)=dcos(time/2.0d0/ep**2)*wp(1,m)-dsin(time/2.0d0/ep**2)*wp(2,m)
    vp(2,m)=dcos(time/2.0d0/ep**2)*wp(2,m)+dsin(time/2.0d0/ep**2)*wp(1,m)
    vm(1,m)=dcos(time/2.0d0/ep**2)*wm(1,m)+dsin(time/2.0d0/ep**2)*wm(2,m)
    vm(2,m)=dcos(time/2.0d0/ep**2)*wm(2,m)-dsin(time/2.0d0/ep**2)*wm(1,m)
    z=(vp(:,m)+vm(:,m))/2.0d0
    xxt(1)=dreal(dcos(time/2.0d0/ep**2)*z(1)+dsin(time/2.0d0/ep**2)*z(2))
    xxt(2)=dreal(dcos(time/2.0d0/ep**2)*z(2)-dsin(time/2.0d0/ep**2)*z(1))
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
enddo
call calcul_rho_m6( p, f )
call sll_s_fft_free(PlnF)
call sll_s_fft_free(PlnB)


open(unit=851,file='fh64.dat')
!do n=0,nstep
!write(851,*)energy(n)
!enddo
!deallocate (energy)
do i=1,nx
do j=1,ny
write(851,*)f%r0(i,j)
enddo
enddo
close(851)


contains

subroutine apply_bc()
do while ( xxt(1) > xmax )
xxt(1) = xxt(1) - dimx
enddo
do while ( xxt(1) < xmin )
xxt(1)= xxt(1) + dimx
enddo
do while ( xxt(2) > ymax )
xxt(2)  = xxt(2)  - dimy
enddo
do while ( xxt(2)  < ymin )
xxt(2) = xxt(2)  + dimy
enddo
end subroutine apply_bc

subroutine energyuse()
do m=1,npp
    call sll_s_fft_exec_c2c_1d(PlnF, up(1,:,m),temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, up(2,:,m),temptilde(2,:))
    wp(:,m)=0.0d0
    do n=0,ntau-1
        wp(:,m)=wp(:,m)+temptilde(:,n)/Ntau*cdexp(sll_p_i1*ltau(n)*time/2.0d0/ep**2)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, um(1,:,m),temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, um(2,:,m),temptilde(2,:))
    wm(:,m)=0.0d0
    do n=0,ntau-1
        wm(:,m)=wm(:,m)+temptilde(:,n)/Ntau*cdexp(sll_p_i1*ltau(n)*time/2.0d0/ep**2)
    enddo
    vp(1,m)=dcos(time/2.0d0/ep**2)*wp(1,m)-dsin(time/2.0d0/ep**2)*wp(2,m)
    vp(2,m)=dcos(time/2.0d0/ep**2)*wp(2,m)+dsin(time/2.0d0/ep**2)*wp(1,m)
    vm(1,m)=dcos(time/2.0d0/ep**2)*wm(1,m)+dsin(time/2.0d0/ep**2)*wm(2,m)
    vm(2,m)=dcos(time/2.0d0/ep**2)*wm(2,m)-dsin(time/2.0d0/ep**2)*wm(1,m)
    z=(vp(:,m)+vm(:,m))/2.0d0
    xxt(1)=dreal(dcos(time/2.0d0/ep**2)*z(1)+dsin(time/2.0d0/ep**2)*z(2))
    xxt(2)=dreal(dcos(time/2.0d0/ep**2)*z(2)-dsin(time/2.0d0/ep**2)*z(1))
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    p%vpx(m)=(dsin(time/2.0d0/ep**2)*vm(1,m)+dcos(time/2.0d0/ep**2)*vm(2,m))/2.0d0/ep
    p%vpy(m)=-(-dsin(time/2.0d0/ep**2)*vm(2,m)+dcos(time/2.0d0/ep**2)*vm(1,m))/2.0d0/ep
enddo
call calcul_rho_m6( p, f )
call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
open(10, file="energy.dat", position="append")
if (istep ==1) rewind(10)
write(10,*) time, 0.5*log(sum(f%ex**2)), 0.5*log(sum(f%ey**2))
close(10)
end subroutine energyuse


end program test_pic2d

!subroutine morecorrection4th(ep,ntau,tau,npp,up,um,xtemp1)
!use sll_m_poisson_2d_base
!use sll_m_poisson_2d_periodic
!use particules
!use sll_m_fft
!use sll_m_poisson_2d_base
!use sll_m_poisson_2d_periodic
!use sll_m_constants
!implicit none
!
!type(tm_mesh_fields) :: f
!type(particle)       :: p
!
!type(sll_t_fft) :: PlFx, PlBx,PlFy, PlBy,PlnF,PlnB
!integer(4),intent(in)    :: Npp,Ntau
!complex(8),intent(in) :: up(2,0:ntau-1,npp),um(2,0:ntau-1,npp)
!real(8),intent(in) :: ep,tau(0:ntau-1)
!complex(8),intent(out) :: xtemp1(2,0:ntau-1,npp)
!integer(4)  :: n,m,l,kl
!real(8) :: ltau(0:Ntau-1),lx(0:119),ly(0:9),xxt(2)
!complex(8) :: tempx(0:119),tempy(0:9),temptau(0:ntau-1),tautilde(0:ntau-1),rhot(0:120,0:10,0:ntau-1)
!complex(8) :: temp(0:120,0:10,0:ntau-1)
!complex(8) :: xtemp(2,0:ntau-1,npp),xtemp2(2,0:ntau-1,npp),dxE(0:119,0:9),dyE(0:119,0:9)
!class(sll_c_poisson_2d_base), pointer :: poisson
!sll_int32  :: error
!ny=10
!nx=120
!m=ntau/2
!ltau=(/ (n, n=0,m-1), (n, n=-m,-1 )/)
!m=nx/2
!lx=(/ (n, n=0,m-1), (n, n=-m,-1 )/)
!m=ny/2
!ly=(/ (n, n=0,m-1), (n, n=-m,-1 )/)
!SLL_CLEAR_ALLOCATE(f%ex(0:nx,0:ny), error)
!SLL_CLEAR_ALLOCATE(f%ey(0:nx,0:ny), error)
!SLL_CLEAR_ALLOCATE(f%r0(0:nx,0:ny), error)
!call sll_s_fft_init_c2c_1d(PlnF,Ntau,temptau,tautilde,sll_p_fft_FORWARD,optimization=sll_p_FFT_MEASURE)
!call sll_s_fft_init_c2c_1d(PlnB,Ntau,tautilde,temptau,sll_p_fft_BACKWARD,optimization=sll_p_FFT_MEASURE)
!call sll_s_fft_init_c2c_1d(PlFx,Nx,dxE(:,1),tempx,sll_p_fft_FORWARD,optimization=sll_p_FFT_MEASURE)
!call sll_s_fft_init_c2c_1d(PlBx,Nx,dxE(:,1),tempx,sll_p_FFT_BACKWARD,optimization=sll_p_FFT_MEASURE)
!call sll_s_fft_init_c2c_1d(PlFy,Ny,dxE(1,:),tempy,sll_p_fft_FORWARD,optimization=sll_p_FFT_MEASURE)
!call sll_s_fft_init_c2c_1d(PlBy,Ny,dxE(1,:),tempy,sll_p_FFT_BACKWARD,optimization=sll_p_FFT_MEASURE)
!poisson => sll_f_new_poisson_2d_periodic(0.0d0,dimx,nx,0.0d0,dimy,ny)
!call plasma( p )
!
!do n=0,ntau-1
!    xtemp1(1,n,:)=dcos(tau(n))*up(1,n,:)-dsin(tau(n))*up(2,n,:)!v_+ 2scaled
!    xtemp1(2,n,:)=dcos(tau(n))*up(2,n,:)+dsin(tau(n))*up(1,n,:)
!    xtemp2(1,n,:)=dcos(tau(n))*um(1,n,:)+dsin(tau(n))*um(2,n,:)!v_- 2scaled
!    xtemp2(2,n,:)=dcos(tau(n))*um(2,n,:)-dsin(tau(n))*um(1,n,:)
!
!!    xtemp(1,n,:)=(dcos(tau(n))*xtemp2(2,n,:)+dsin(tau(n))*xtemp2(1,n,:))/2.0d0/ep!y 2scaled
!!    xtemp(2,n,:)=-(dcos(tau(n))*xtemp2(1,n,:)-dsin(tau(n))*xtemp2(2,n,:))/2.0d0/ep
!
!    xtemp1(:,n,:)=(xtemp1(:,n,:)+xtemp2(:,n,:))/2.0d0!z 2scaled
!    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)!x 2scaled
!    xtemp2(2,n,:)=dcos(tau(n))*xtemp1(2,n,:)-dsin(tau(n))*xtemp1(1,n,:)
!    do m=1,nbpart
!!        p%vpx(m)=xtemp(1,n,m)
!!        p%vpy(m)=xtemp(2,n,m)
!        xxt=dreal(xtemp2(:,n,m))
!        call apply_bc()
!        p%idx(m) = floor(xxt(1)/dimx*nx)
!        p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
!        p%idy(m) = floor(xxt(2)/dimy*ny)
!        p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
!    enddo
!    call calcul_energy( p, f, 0.0d0*tau(n)) !changes made to f%rho,ex,ey
!    dxE=f%ex(0:nx-1,0:ny-1)
!    dyE=f%ey(0:nx-1,0:ny-1)
!    call calcul_rho_m6( p, f )
!    rhot(:,:,n)=f%r0
!    do l=0,ny-1
!        call sll_s_fft_exec_c2c_1d(PlFx, dxE(:,l), tempx)
!        do m=0,nx-1
!            tempx(m)=tempx(m)*sll_p_i1*lx(m)/Nx
!        enddo
!        call sll_s_fft_exec_c2c_1d(PlBx, tempx,dxE(:,l))
!    enddo
!    do l=0,nx-1
!        call sll_s_fft_exec_c2c_1d(PlFy, dyE(l,:), tempy)
!        do m=0,ny-1
!            tempy(m)=tempy(m)*sll_p_i1*ly(m)/Ny
!        enddo
!        call sll_s_fft_exec_c2c_1d(PlBy, tempy,dyE(l,:))
!    enddo
!    temp(0:nx-1,0:ny-1,n)=-(dxE+dyE)/ep
!enddo
!do l=0,nx-1
!    do m=0,ny-1
!        temptau=rhot(l,m,:)
!        call sll_s_fft_exec_c2c_1d(PlnF, temptau, tautilde)
!        do kl=0,Ntau-1
!            tautilde(kl)=tautilde(kl)*sll_p_i1*ltau(kl)/Ntau
!        enddo
!        call sll_s_fft_exec_c2c_1d(PlnB, tautilde, temptau)
!        rhot(l,m,:)=-temptau/2.0d0/ep**2+temp(l,m,:)!final dt\rho
!    enddo
!enddo
!rhot(0:nx-1,ny,:)  = rhot(0:nx-1,0,:)
!rhot(nx,0:ny-1,:)  = rhot(0,0:ny-1,:)
!rhot(nx,ny,:)= rhot(0,0,:)
!
!do n=0,ntau-1
!    xtemp1(1,n,:)=dcos(tau(n))*up(1,n,:)-dsin(tau(n))*up(2,n,:)
!    xtemp1(2,n,:)=dsin(tau(n))*up(1,n,:)+dcos(tau(n))*up(2,n,:)
!    xtemp2(1,n,:)=dcos(tau(n))*um(1,n,:)+dsin(tau(n))*um(2,n,:)
!    xtemp2(2,n,:)=-dsin(tau(n))*um(1,n,:)+dcos(tau(n))*um(2,n,:)
!enddo
!xtemp1=(xtemp1+xtemp2)/2.0d0
!do n=0,ntau-1
!    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)
!    xtemp2(2,n,:)=-dsin(tau(n))*xtemp1(1,n,:)+dcos(tau(n))*xtemp1(2,n,:)
!    do m=1,nbpart
!        xxt=dreal(xtemp2(:,n,m))
!        call apply_bc()
!        p%idx(m) = floor(xxt(1)/dimx*nx)
!        p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
!        p%idy(m) = floor(xxt(2)/dimy*ny)
!        p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
!    enddo
!    f%r0=rhot(:,:,n)
!    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
!    call interpol_eb_m6( f, p )
!    xtemp1(1,n,:)=p%epx !g_t(0,tau,U_3(0))
!    xtemp1(2,n,:)=p%epy
!enddo
!
!call sll_s_fft_free(PlnF)
!call sll_s_fft_free(PlnB)
!call sll_s_fft_free(PlFx)
!call sll_s_fft_free(PlBx)
!call sll_s_fft_free(PlFy)
!call sll_s_fft_free(PlBy)
!
!contains
!
!subroutine apply_bc()
!do while ( xxt(1) > dimx)
!xxt(1) = xxt(1) - dimx
!enddo
!do while ( xxt(1) < 0.0d0 )
!xxt(1)= xxt(1) + dimx
!enddo
!do while ( xxt(2) > dimy )
!xxt(2)  = xxt(2)  - dimy
!enddo
!do while ( xxt(2)  < 0.0d0 )
!xxt(2) = xxt(2)  + dimy
!enddo
!end subroutine apply_bc
!end subroutine morecorrection4th


!!--more correction--
!do n=0,ntau-1
!    xtemp1(1,n,:)=dcos(tau(n))*up(1,n,:)-dsin(tau(n))*up(2,n,:)
!    xtemp1(2,n,:)=dsin(tau(n))*up(1,n,:)+dcos(tau(n))*up(2,n,:)
!    xtemp2(1,n,:)=dcos(tau(n))*um(1,n,:)+dsin(tau(n))*um(2,n,:)
!    xtemp2(2,n,:)=-dsin(tau(n))*um(1,n,:)+dcos(tau(n))*um(2,n,:)
!enddo
!xtemp1=(xtemp1+xtemp2)/2.0d0
!do n=0,ntau-1
!    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)
!    xtemp2(2,n,:)=-dsin(tau(n))*xtemp1(1,n,:)+dcos(tau(n))*xtemp1(2,n,:)
!    do m=1,nbpart
!    xxt=dreal(xtemp2(:,n,m))
!    call apply_bc()
!    p%idx(m) = floor(xxt(1)/dimx*nx)
!    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
!    p%idy(m) = floor(xxt(2)/dimy*ny)
!    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
!    enddo
!    f%ex=fex(:,:,n)
!    f%ey=fey(:,:,n)
!    call interpol_eb_m6( f, p )
!    Et(1,n,:)=p%epx !g_3(0,tau,U_3(0))
!    Et(2,n,:)=p%epy
!enddo
!do m=1,npp
!    xtemp1(1,:,m)=2.0d0*Et(2,:,m)
!    xtemp1(2,:,m)=-2.0d0*Et(1,:,m)!g_+
!    call sll_s_fft_exec_c2c_1d(PlnF, xtemp1(1,:,m), temptilde(1,:))
!    call sll_s_fft_exec_c2c_1d(PlnF, xtemp1(2,:,m), temptilde(2,:))
!    do n=1,Ntau-1
!    temptilde(:,n)=-sll_p_i1*temptilde(:,n)/ltau(n)/Ntau
!    enddo
!    temptilde(:,0)=0.0d0
!    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(1,:), temp(1,:))
!    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(2,:), temp(2,:))!AF+
!    up(1,:,m)=wp(1,m)+2.0d0*ep**2*(temp(1,:)-temp(1,0))
!    up(2,:,m)=wp(2,m)+2.0d0*ep**2*(temp(2,:)-temp(2,0))!4 ini data of U_+
!    !---
!    do n=0,ntau-1
!        xtemp2(1,n,m)=-2.0d0*(dsin(2.0d0*tau(n))*Et(1,n,m)+dcos(2.0d0*tau(n))*Et(2,n,m))
!        xtemp2(2,n,m)=2.0d0*(-dsin(2.0d0*tau(n))*Et(2,n,m)+dcos(2.0d0*tau(n))*Et(1,n,m))!g_-
!    enddo
!    call sll_s_fft_exec_c2c_1d(PlnF, xtemp2(1,:,m), temptilde(1,:))
!    call sll_s_fft_exec_c2c_1d(PlnF, xtemp2(2,:,m), temptilde(2,:))
!    do n=1,Ntau-1
!    temptilde(:,n)=-sll_p_i1*temptilde(:,n)/ltau(n)/Ntau
!    enddo
!    temptilde(:,0)=0.0d0
!    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(1,:), temp(1,:))
!    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(2,:), temp(2,:))!AF-
!    um(1,:,m)=wm(1,m)+2.0d0*ep**2*(temp(1,:)-temp(1,0))
!    um(2,:,m)=wm(2,m)+2.0d0*ep**2*(temp(2,:)-temp(2,0))!4 ini data of U_-
!enddo
!!------
!do n=0,ntau-1
!    xtemp1(1,n,:)=dcos(tau(n))*up(1,n,:)-dsin(tau(n))*up(2,n,:)!v_+ 2scaled
!    xtemp1(2,n,:)=dcos(tau(n))*up(2,n,:)+dsin(tau(n))*up(1,n,:)
!    xtemp2(1,n,:)=dcos(tau(n))*um(1,n,:)+dsin(tau(n))*um(2,n,:)!v_- 2scaled
!    xtemp2(2,n,:)=dcos(tau(n))*um(2,n,:)-dsin(tau(n))*um(1,n,:)
!    xtemp1(:,n,:)=(xtemp1(:,n,:)+xtemp2(:,n,:))/2.0d0!z 2scaled
!    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)!x 2scaled
!    xtemp2(2,n,:)=dcos(tau(n))*xtemp1(2,n,:)-dsin(tau(n))*xtemp1(1,n,:)
!    do m=1,nbpart
!    xxt=dreal(xtemp2(:,n,m))
!    call apply_bc()
!    p%idx(m) = floor(xxt(1)/dimx*nx)
!    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
!    p%idy(m) = floor(xxt(2)/dimy*ny)
!    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
!    enddo
!    call calcul_rho_m6( p, f )
!    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
!    fex(:,:,n)=f%ex
!    fey(:,:,n)=f%ey!E_4(0,x)
!enddo

!----not working--4th order correction, adding derivative----
!call morecorrection4th(ep,ntau,tau,npp,up,um,gp)!xtemp1 store \partial_t g(0,tau,U_3)
!!gp=0.0d0
!do n=0,ntau-1
!    xtemp1(1,n,:)=dcos(tau(n))*up0(1,0,:)-dsin(tau(n))*up0(2,0,:)
!    xtemp1(2,n,:)=dsin(tau(n))*up0(1,0,:)+dcos(tau(n))*up0(2,0,:)
!    xtemp2(1,n,:)=dcos(tau(n))*um0(1,0,:)+dsin(tau(n))*um0(2,0,:)
!    xtemp2(2,n,:)=-dsin(tau(n))*um0(1,0,:)+dcos(tau(n))*um0(2,0,:)
!enddo
!xtemp1=(xtemp1+xtemp2)/2.0d0
!do n=0,ntau-1
!    xtemp2(1,n,:)=dcos(tau(n))*xtemp1(1,n,:)+dsin(tau(n))*xtemp1(2,n,:)
!    xtemp2(2,n,:)=-dsin(tau(n))*xtemp1(1,n,:)+dcos(tau(n))*xtemp1(2,n,:)
!    do  m=1,npp
!        Et(1,n,m)=alpha*dcos(kx*dreal(xtemp2(1,n,m)))
!        Et(2,n,m)=0.0d0
!    enddo
!enddo
!do m=1,npp
!    xtemp1(1,:,m)=2.0d0*(Et(2,:,m)+gp(2,:,m))
!    xtemp1(2,:,m)=-2.0d0*(Et(1,:,m)+gp(1,:,m))!(dtF+dUFPiF) g_+
!    call sll_s_fft_exec_c2c_1d(PlnF, xtemp1(1,:,m), temptilde(1,:))
!    call sll_s_fft_exec_c2c_1d(PlnF, xtemp1(2,:,m), temptilde(2,:))
!    do n=1,Ntau-1
!        temptilde(:,n)=-temptilde(:,n)/ltau(n)**2/Ntau
!    enddo
!    temptilde(:,0)=0.0d0
!    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(1,:), up0(1,:,m))
!    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(2,:), up0(2,:,m))!A^2F+
!    do n=0,ntau-1
!        xtemp2(1,n,m)=-2.0d0*(dsin(2.0d0*tau(n))*(Et(1,n,m)+gp(1,n,m))+dcos(2.0d0*tau(n))*(Et(2,n,m)+gp(2,n,m)))
!        xtemp2(2,n,m)=2.0d0*(-dsin(2.0d0*tau(n))*(Et(2,n,m)+gp(2,n,m))+dcos(2.0d0*tau(n))*(Et(1,n,m)+gp(1,n,m)))!g_-
!    enddo
!    call sll_s_fft_exec_c2c_1d(PlnF, xtemp2(1,:,m), temptilde(1,:))
!    call sll_s_fft_exec_c2c_1d(PlnF, xtemp2(2,:,m), temptilde(2,:))
!    do n=1,Ntau-1
!        temptilde(:,n)=-temptilde(:,n)/ltau(n)**2/Ntau
!    enddo
!    temptilde(:,0)=0.0d0
!    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(1,:), um0(1,:,m))
!    call sll_s_fft_exec_c2c_1d(PlnB, temptilde(2,:), um0(2,:,m))!AF-
!enddo
!do n=0,ntau-1
!    up(:,n,:)=up(:,n,:)-4.0d0*ep**4*(up0(:,n,:)-up0(:,0,:))
!    um(:,n,:)=um(:,n,:)-4.0d0*ep**4*(um0(:,n,:)-um0(:,0,:))
!enddo
!--------

