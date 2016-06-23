!NKS-two-scale solver for 4d VP
!2nd order EWI scheme
!two-scale E
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
use particules_m6, only: calcul_rho_m6, interpol_eb_m6

implicit none

type(tm_mesh_fields)  :: f
type(particle)        :: p

type(sll_t_fft)       :: PlnF
type(sll_t_fft)       :: PlnB
integer(4), parameter :: Ntau=32
integer(4), parameter :: npp=204800
real(8)               :: xxt(2)
real(8)               :: ep
real(8)               :: dtau
real(8)               :: tau(0:Ntau-1)
real(8)               :: ltau(0:Ntau-1)
real(8)               :: auxpx(2,npp)
complex(8)            :: pl(0:ntau-1)
complex(8)            :: ql(0:ntau-1)
complex(8)            :: gp(2,0:Ntau-1,npp)
complex(8)            :: gm(2,0:Ntau-1,npp)
complex(8)            :: wp(2,npp)
complex(8)            :: wm(2,npp)
complex(8)            :: xtemp1(2,0:ntau-1,npp)
complex(8)            :: xtemp2(2,0:ntau-1,npp)
complex(8)            :: Et(2,0:ntau-1,npp)
complex(8)            :: temp(2,0:Ntau-1)
complex(8)            :: temptilde(2,0:Ntau-1)
complex(8)            :: up(2,0:ntau-1,npp)
complex(8)            :: um(2,0:ntau-1,npp)
complex(8)            :: vp(2,npp)
complex(8)            :: vm(2,npp)
complex(8)            :: z(2)
complex(8)            :: fex(0:64,0:32,0:ntau-1)
complex(8)            :: fey(0:64,0:32,0:ntau-1)
complex(8)            :: up0(2,0:ntau-1,npp)
complex(8)            :: um0(2,0:ntau-1,npp)

sll_real64            :: time
sll_real64            :: xmin
sll_real64            :: xmax
sll_real64            :: ymin
sll_real64            :: ymax
sll_int32             :: istep
sll_int32             :: iplot
sll_int32             :: iargc
sll_int32             :: n,m
sll_int32             :: i
sll_int32             :: j
sll_int32             :: error

sll_real64            :: aux1, aux2
sll_real64            :: t1
sll_real64            :: s, dum

real    :: start_time, stop_time
integer :: dat_file_id, ref_file_id
logical :: file_exists
real(8) :: cost, sint


character(len=272)    :: argv
class(sll_c_poisson_2d_base), pointer :: poisson

call cpu_time(start_time)
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
call sll_s_fft_init_c2c_1d(PlnF,Ntau,temp(1,:), &
  temptilde(1,:),sll_p_fft_forward,optimization=sll_p_fft_measure)
call sll_s_fft_init_c2c_1d(PlnB,Ntau,temp(1,:), &
  temptilde(1,:),sll_p_fft_backward,optimization=sll_p_fft_measure)
m=ntau/2
ltau=(/ (n, n=0,m-1), (n, n=-m,-1 )/)
pl(0)=dt
ql(0)=dt**2/2.0d0
do i=1,Ntau-1
  pl(i)=2.0d0*ep**2*sll_p_i1*(exp(-sll_p_i1*ltau(i)*dt/2.0d0/ep**2)-1.0d0)/ltau(i)
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
print"('ep = ', g15.3)", ep
print"('nbpart = ', g15.3)", nbpart
print"('dt = ', g15.3)", dt
poisson => sll_f_new_poisson_2d_periodic(xmin,xmax,nx,ymin,ymax,ny)

auxpx(1,:)=(p%dpx+p%idx)*dx
auxpx(2,:)=(p%dpy+p%idy)*dy
wp(1,:)=2.0d0*(auxpx(1,:)+ep*p%vpy)
wp(2,:)=2.0d0*(auxpx(2,:)-ep*p%vpx)
wm(1,:)=-2.0d0*ep*p%vpy
wm(2,:)=2.0d0*ep*p%vpx

do n=0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))
  xtemp1(1,n,:)=cost*wp(1,:)-sint*wp(2,:)
  xtemp1(2,n,:)=sint*wp(1,:)+cost*wp(2,:)
  xtemp2(1,n,:)=cost*wm(1,:)+sint*wm(2,:)
  xtemp2(2,n,:)=-sint*wm(1,:)+cost*wm(2,:)
enddo

xtemp1=(xtemp1+xtemp2)/2.0d0

do n=0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))
  xtemp2(1,n,:)=cost*xtemp1(1,n,:)+sint*xtemp1(2,n,:)
  xtemp2(2,n,:)=-sint*xtemp1(1,n,:)+cost*xtemp1(2,n,:)
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
    cost = cos(2.0d0*tau(n))
    sint = sin(2.0d0*tau(n))
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
  cost = cos(tau(n))
  sint = sin(tau(n))
  xtemp1(1,n,:)=cost*up(1,n,:)-sint*up(2,n,:)!v_+ 2scaled
  xtemp1(2,n,:)=cost*up(2,n,:)+sint*up(1,n,:)
  xtemp2(1,n,:)=cost*um(1,n,:)+sint*um(2,n,:)!v_- 2scaled
  xtemp2(2,n,:)=cost*um(2,n,:)-sint*um(1,n,:)
  xtemp1(:,n,:)=(xtemp1(:,n,:)+xtemp2(:,n,:))/2.0d0!z 2scaled
  xtemp2(1,n,:)=cost*xtemp1(1,n,:)+sint*xtemp1(2,n,:)!x 2scaled
  xtemp2(2,n,:)=cost*xtemp1(2,n,:)-sint*xtemp1(1,n,:)

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
  cost = cos(tau(n))
  sint = sin(tau(n))
  xtemp1(1,n,:)=cost*up(1,n,:)-sint*up(2,n,:)
  xtemp1(2,n,:)=sint*up(1,n,:)+cost*up(2,n,:)
  xtemp2(1,n,:)=cost*um(1,n,:)+sint*um(2,n,:)
  xtemp2(2,n,:)=-sint*um(1,n,:)+cost*um(2,n,:)
enddo

xtemp1=(xtemp1+xtemp2)/2.0d0

do n=0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))
  xtemp2(1,n,:)=cost*xtemp1(1,n,:)+sint*xtemp1(2,n,:)
  xtemp2(2,n,:)=-sint*xtemp1(1,n,:)+cost*xtemp1(2,n,:)
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
    cost = cos(2d0*tau(n))
    sint = sin(2d0*tau(n))
    xtemp2(1,n,m)=-2.0d0*(sint*Et(1,n,m)+cost*Et(2,n,m))
    xtemp2(2,n,m)=2.0d0*(-sint*Et(2,n,m)+cost*Et(1,n,m))!g_-
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
  cost = cos(tau(n))
  sint = sin(tau(n))
  xtemp1(1,n,:)=cost*up(1,n,:)-sint*up(2,n,:)!v_+ 2scaled
  xtemp1(2,n,:)=cost*up(2,n,:)+sint*up(1,n,:)
  xtemp2(1,n,:)=cost*um(1,n,:)+sint*um(2,n,:)!v_- 2scaled
  xtemp2(2,n,:)=cost*um(2,n,:)-sint*um(1,n,:)
  xtemp1(:,n,:)=0.5d0*(xtemp1(:,n,:)+xtemp2(:,n,:))!z 2scaled
  xtemp2(1,n,:)=cost*xtemp1(1,n,:)+sint*xtemp1(2,n,:)!x 2scaled
  xtemp2(2,n,:)=cost*xtemp1(2,n,:)-sint*xtemp1(1,n,:)
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
  cost = cos(tau(n))
  sint = sin(tau(n))
  xtemp1(1,n,:) =  cost*up(1,n,:) - sint*up(2,n,:)
  xtemp1(2,n,:) =  sint*up(1,n,:) + cost*up(2,n,:)
  xtemp2(1,n,:) =  cost*um(1,n,:) + sint*um(2,n,:)
  xtemp2(2,n,:) = -sint*um(1,n,:) + cost*um(2,n,:)
enddo
xtemp1=(xtemp1+xtemp2)/2.0d0
do n=0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))
  xtemp2(1,n,:)=cost*xtemp1(1,n,:)+sint*xtemp1(2,n,:)
  xtemp2(2,n,:)=-sint*xtemp1(1,n,:)+cost*xtemp1(2,n,:)
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
    cost = cos(2d0*tau(n))
    sint = sin(2d0*tau(n))
    temp(1,n)=-2.0d0*(sint*Et(1,n,m)+cost*Et(2,n,m))
    temp(2,n)=2.0d0*(-sint*Et(2,n,m)+cost*Et(1,n,m))
  enddo
  call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), temptilde(1,:))
  call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), temptilde(2,:))
  xtemp2(:,:,m)=temptilde/ntau!g_-tilde(t=0)
enddo

do m=1,npp
  call sll_s_fft_exec_c2c_1d(PlnF, up(1,:,m), temptilde(1,:))
  call sll_s_fft_exec_c2c_1d(PlnF, up(2,:,m), temptilde(2,:))
  do n=0,ntau-1
    temp(:,n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/ep**2)*temptilde(:,n)/ntau+pl(n)*xtemp1(:,n,m)!utilde_+^1,predict
  enddo
  call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), up0(1,:,m))!u_+(t1),predict
  call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), up0(2,:,m))
  call sll_s_fft_exec_c2c_1d(PlnF, um(1,:,m), temptilde(1,:))
  call sll_s_fft_exec_c2c_1d(PlnF, um(2,:,m), temptilde(2,:))
  do n=0,ntau-1
    temp(:,n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/ep**2)*temptilde(:,n)/ntau+pl(n)*xtemp2(:,n,m)!utilde_-^1,predict
  enddo
  call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), um0(1,:,m))!u_-(t1),predict
  call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), um0(2,:,m))
enddo

gp=xtemp1
gm=xtemp2
do n=0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))
  xtemp1(1,n,:)=cost*up0(1,n,:)-sint*up0(2,n,:)!v_+ 2scaled
  xtemp1(2,n,:)=cost*up0(2,n,:)+sint*up0(1,n,:)
  xtemp2(1,n,:)=cost*um0(1,n,:)+sint*um0(2,n,:)!v_- 2scaled
  xtemp2(2,n,:)=cost*um0(2,n,:)-sint*um0(1,n,:)
  xtemp1(:,n,:)=(xtemp1(:,n,:)+xtemp2(:,n,:))*0.5d0!z 2scaled
  xtemp2(1,n,:)=cost*xtemp1(1,n,:)+sint*xtemp1(2,n,:)!x 2scaled
  xtemp2(2,n,:)=cost*xtemp1(2,n,:)-sint*xtemp1(1,n,:)
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
  cost = cos(tau(n))
  sint = sin(tau(n))
  xtemp1(1,n,:) =  cost*up0(1,n,:) - sint*up0(2,n,:)
  xtemp1(2,n,:) =  sint*up0(1,n,:) + cost*up0(2,n,:)
  xtemp2(1,n,:) =  cost*um0(1,n,:) + sint*um0(2,n,:)
  xtemp2(2,n,:) = -sint*um0(1,n,:) + cost*um0(2,n,:)
enddo
xtemp1=(xtemp1+xtemp2)/2.0d0
do n=0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))
  xtemp2(1,n,:)=cost*xtemp1(1,n,:)+sint*xtemp1(2,n,:)
  xtemp2(2,n,:)=-sint*xtemp1(1,n,:)+cost*xtemp1(2,n,:)
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
    cost = cos(2d0*tau(n))
    sint = sin(2d0*tau(n))
    temp(1,n)=-2.0d0*(sint*Et(1,n,m)+cost*Et(2,n,m))
    temp(2,n)=2.0d0*(-sint*Et(2,n,m)+cost*Et(1,n,m))
  enddo
  call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), temptilde(1,:))
  call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), temptilde(2,:))
  xtemp2(:,:,m)=temptilde/ntau!g_-tilde(t1) predict
enddo

do m=1,npp
  call sll_s_fft_exec_c2c_1d(PlnF, up(1,:,m), temptilde(1,:))
  call sll_s_fft_exec_c2c_1d(PlnF, up(2,:,m), temptilde(2,:))
  do n=0,ntau-1
    temp(:,n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/ep**2)*temptilde(:,n)/ntau &
             +pl(n)*xtemp1(:,n,m) &
             +ql(n)*(xtemp1(:,n,m)-gp(:,n,m))/dt
  enddo
  call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), up(1,:,m))!u_+(t1)
  call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), up(2,:,m))
  call sll_s_fft_exec_c2c_1d(PlnF, um(1,:,m), temptilde(1,:))
  call sll_s_fft_exec_c2c_1d(PlnF, um(2,:,m), temptilde(2,:))
  do n=0,ntau-1
    temp(:,n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/ep**2)*temptilde(:,n)/ntau &
             +pl(n)*xtemp2(:,n,m) &
             +ql(n)*(xtemp2(:,n,m)-gm(:,n,m))/dt
  enddo
  call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), um(1,:,m))!u_-(t1)
  call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), um(2,:,m))
enddo

do n=0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))
  xtemp1(1,n,:)=cost*up(1,n,:)-sint*up(2,n,:)!v_+ 2scaled
  xtemp1(2,n,:)=cost*up(2,n,:)+sint*up(1,n,:)
  xtemp2(1,n,:)=cost*um(1,n,:)+sint*um(2,n,:)!v_- 2scaled
  xtemp2(2,n,:)=cost*um(2,n,:)-sint*um(1,n,:)
  xtemp1(:,n,:)=0.5_f64*(xtemp1(:,n,:)+xtemp2(:,n,:))!z 2scaled
  xtemp2(1,n,:)=cost*xtemp1(1,n,:)+sint*xtemp1(2,n,:)!x 2scaled
  xtemp2(2,n,:)=cost*xtemp1(2,n,:)-sint*xtemp1(1,n,:)
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

do istep = 2, nstep
  do n=0,ntau-1
    cost = cos(tau(n))
    sint = sin(tau(n))
    xtemp1(1,n,:)= cost*up(1,n,:)-sint*up(2,n,:)
    xtemp1(2,n,:)= sint*up(1,n,:)+cost*up(2,n,:)
    xtemp2(1,n,:)= cost*um(1,n,:)+sint*um(2,n,:)
    xtemp2(2,n,:)=-sint*um(1,n,:)+cost*um(2,n,:)
  enddo
  xtemp1=(xtemp1+xtemp2)/2.0d0
  do n=0,ntau-1
    cost = cos(tau(n))
    sint = sin(tau(n))
    xtemp2(1,n,:)=cost*xtemp1(1,n,:)+sint*xtemp1(2,n,:)
    xtemp2(2,n,:)=-sint*xtemp1(1,n,:)+cost*xtemp1(2,n,:)
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
    Et(1,n,:)=p%epx 
    Et(2,n,:)=p%epy
  enddo
  do m=1,npp
    temp(1,:)=2.0d0*Et(2,:,m)
    temp(2,:)=-2.0d0*Et(1,:,m)
    call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), temptilde(2,:))
    xtemp1(:,:,m)=temptilde/ntau
    !---
    do n=0,ntau-1
      cost = cos(2d0*tau(n))
      sint = sin(2d0*tau(n))
      temp(1,n)=-2.0d0*(sint*Et(1,n,m)+cost*Et(2,n,m))
      temp(2,n)=2.0d0*(-sint*Et(2,n,m)+cost*Et(1,n,m))
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), temptilde(2,:))
    xtemp2(:,:,m)=temptilde/ntau
  enddo

  do m=1,npp
    call sll_s_fft_exec_c2c_1d(PlnF, up(1,:,m), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, up(2,:,m), temptilde(2,:))
    do n=0,ntau-1
      temp(:,n)= exp(-sll_p_i1*ltau(n)*dt*0.5d0/ep**2)*temptilde(:,n)/ntau &
               + pl(n)*xtemp1(:,n,m) &
               + ql(n)*(xtemp1(:,n,m)-gp(:,n,m))/dt
    enddo
    call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), up(1,:,m))
    call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), up(2,:,m))
    call sll_s_fft_exec_c2c_1d(PlnF, um(1,:,m), temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, um(2,:,m), temptilde(2,:))
    do n=0,ntau-1
      temp(:,n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/ep**2)*temptilde(:,n)/ntau &
               +pl(n)*xtemp2(:,n,m) &
               +ql(n)*(xtemp2(:,n,m)-gm(:,n,m))/dt
    enddo
    call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), um(1,:,m))
    call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), um(2,:,m))
  enddo
  gp=xtemp1
  gm=xtemp2
  !--updata E--
  time=dt*istep
  do n=0,ntau-1
    cost = cos(tau(n))
    sint = sin(tau(n))
    xtemp1(1,n,:)=cost*up(1,n,:)-sint*up(2,n,:)!v_+ 2scaled
    xtemp1(2,n,:)=cost*up(2,n,:)+sint*up(1,n,:)
    xtemp2(1,n,:)=cost*um(1,n,:)+sint*um(2,n,:)!v_- 2scaled
    xtemp2(2,n,:)=cost*um(2,n,:)-sint*um(1,n,:)
    xtemp1(:,n,:)=0.5_f64*(xtemp1(:,n,:)+xtemp2(:,n,:))!z 2scaled
    xtemp2(1,n,:)=cost*xtemp1(1,n,:)+sint*xtemp1(2,n,:)!x 2scaled
    xtemp2(2,n,:)=cost*xtemp1(2,n,:)-sint*xtemp1(1,n,:)
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
enddo

do m=1,npp
  call sll_s_fft_exec_c2c_1d(PlnF, up(1,:,m),temptilde(1,:))
  call sll_s_fft_exec_c2c_1d(PlnF, up(2,:,m),temptilde(2,:))
  wp(:,m)=0.0d0
  do n=0,ntau-1
    wp(:,m)=wp(:,m)+temptilde(:,n)/Ntau*cdexp(sll_p_i1*ltau(n)*time*0.5d0/ep**2)
  enddo
  call sll_s_fft_exec_c2c_1d(PlnF, um(1,:,m),temptilde(1,:))
  call sll_s_fft_exec_c2c_1d(PlnF, um(2,:,m),temptilde(2,:))
  wm(:,m)=0.0d0
  do n=0,ntau-1
    wm(:,m)=wm(:,m)+temptilde(:,n)/Ntau*cdexp(sll_p_i1*ltau(n)*time*0.5d0/ep**2)
  enddo
  cost = dcos(0.5_f64*time/ep**2)
  sint = dsin(0.5_f64*time/ep**2)
  vp(1,m)=cost*wp(1,m)-sint*wp(2,m)
  vp(2,m)=cost*wp(2,m)+sint*wp(1,m)
  vm(1,m)=cost*wm(1,m)+sint*wm(2,m)
  vm(2,m)=cost*wm(2,m)-sint*wm(1,m)
  z=0.5_f64*(vp(:,m)+vm(:,m))
  xxt(1)=dreal(cost*z(1)+sint*z(2))
  xxt(2)=dreal(cost*z(2)-sint*z(1))
  call apply_bc()
  p%idx(m) = floor(xxt(1)/dimx*nx)
  p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
  p%idy(m) = floor(xxt(2)/dimy*ny)
  p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
enddo

call calcul_rho_m6( p, f )
call sll_s_fft_free(PlnF)
call sll_s_fft_free(PlnB)

call cpu_time(stop_time)

s = 0.0_f64
inquire( file= 'fh64.ref', exist=file_exists )
if (file_exists) then
  open(newunit=ref_file_id,file='fh64.ref')
  do i=1,nx
    do j=1,ny
      read(ref_file_id,*) dum
      s = s + abs(dum-f%r0(i,j))
    enddo
    read(ref_file_id,*)
  enddo
else
  open(newunit=dat_file_id,file='fh64.dat')
  do i=1,nx
    do j=1,ny
      write(dat_file_id,*)f%r0(i,j)
    enddo
    write(dat_file_id,*)
  enddo
  close(851)
endif
print *, "CPU time:", stop_time - start_time, "seconds, error = ", s

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

  cost = dcos(0.5_f64*time/ep**2)
  sint = dsin(0.5_f64*time/ep**2)
  
  do m=1,npp
    call sll_s_fft_exec_c2c_1d(PlnF, up(1,:,m),temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, up(2,:,m),temptilde(2,:))
    wp(:,m)=0.0d0
    do n=0,ntau-1
      wp(:,m)=wp(:,m)+temptilde(:,n)/Ntau*cdexp(sll_p_i1*ltau(n)*time*0.5d0/ep**2)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, um(1,:,m),temptilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, um(2,:,m),temptilde(2,:))
    wm(:,m)=0.0d0
    do n=0,ntau-1
      wm(:,m)=wm(:,m)+temptilde(:,n)/Ntau*cdexp(sll_p_i1*ltau(n)*time*0.5d0/ep**2)
    enddo
    vp(1,m) = cost*wp(1,m) - sint*wp(2,m)
    vp(2,m) = cost*wp(2,m) + sint*wp(1,m)
    vm(1,m) = cost*wm(1,m) + sint*wm(2,m)
    vm(2,m) = cost*wm(2,m) - sint*wm(1,m)
    z       = 0.5d0*(vp(:,m)+vm(:,m))
    xxt(1)  = dreal(cost*z(1)+sint*z(2))
    xxt(2)  = dreal(cost*z(2)-sint*z(1))
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    p%vpx(m) =  ( sint*vm(1,m)+cost*vm(2,m))*0.5d0/ep
    p%vpy(m) = -(-sint*vm(2,m)+cost*vm(1,m))*0.5d0/ep
  enddo
  
  call calcul_rho_m6( p, f )
  call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)

end subroutine energyuse

end program test_pic2d

