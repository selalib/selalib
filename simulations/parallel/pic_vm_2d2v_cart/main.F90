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
!$ use omp_lib

implicit none

type(tm_mesh_fields)  :: f
type(particle)        :: p

type(sll_t_fft)       :: fw
type(sll_t_fft)       :: bw
integer(4), parameter :: Ntau=32
integer(4), parameter :: npp=204800
real(8)               :: xxt(2)
real(8)               :: ep
real(8)               :: epsq
real(8)               :: dtau

real(8)    , allocatable:: tau(:)
real(8)    , allocatable:: ltau(:)
complex(8) , allocatable:: pl(:)
complex(8) , allocatable:: ql(:)
complex(8) , allocatable:: gp1(:,:)
complex(8) , allocatable:: gp2(:,:)
complex(8) , allocatable:: gm1(:,:)
complex(8) , allocatable:: gm2(:,:)
complex(8) , allocatable:: wp1(:)
complex(8) , allocatable:: wp2(:)
complex(8) , allocatable:: wm1(:)
complex(8) , allocatable:: wm2(:)
complex(8) , allocatable:: xt1(:,:)
complex(8) , allocatable:: xt2(:,:)
complex(8) , allocatable:: Et(:,:,:)
complex(8) , allocatable:: temp1(:)
complex(8) , allocatable:: temp2(:)
complex(8) , allocatable:: z(:)
complex(8) , allocatable:: fex(:,:,:)
complex(8) , allocatable:: fey(:,:,:)

complex(8) , allocatable:: up(:,:,:)
complex(8) , allocatable:: um(:,:,:)
complex(8) , allocatable:: up0(:,:,:)
complex(8) , allocatable:: um0(:,:,:)

complex(8) :: utmp
complex(8) :: vtmp

sll_real64            :: time
sll_real64            :: xmin
sll_real64            :: xmax
sll_real64            :: ymin
sll_real64            :: ymax
sll_int32             :: istep
sll_int32             :: iargc
sll_int32             :: n,m
sll_int32             :: i
sll_int32             :: j
sll_int32             :: error

sll_real64            :: aux1, aux2
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

SLL_ALLOCATE(tau(0:Ntau-1),           error)
SLL_ALLOCATE(ltau(0:Ntau-1),          error)
SLL_ALLOCATE(pl(0:ntau-1),            error)
SLL_ALLOCATE(ql(0:ntau-1),            error)
SLL_ALLOCATE(wp1(npp),                error)
SLL_ALLOCATE(wp2(npp),                error)
SLL_ALLOCATE(wm1(npp),                error)
SLL_ALLOCATE(wm2(npp),                error)
SLL_ALLOCATE(Et(npp,0:ntau-1,2),      error)
SLL_ALLOCATE(temp1(0:Ntau-1),         error)
SLL_ALLOCATE(temp2(0:Ntau-1),         error)

SLL_ALLOCATE(gp1(0:Ntau-1,npp),       error)
SLL_ALLOCATE(gp2(0:Ntau-1,npp),       error)
SLL_ALLOCATE(gm1(0:Ntau-1,npp),       error)
SLL_ALLOCATE(gm2(0:Ntau-1,npp),       error)
SLL_ALLOCATE(up(0:ntau-1,npp,2),      error)
SLL_ALLOCATE(um(0:ntau-1,npp,2),      error)
SLL_ALLOCATE(up0(0:ntau-1,npp,2),     error)
SLL_ALLOCATE(um0(0:ntau-1,npp,2),     error)
SLL_ALLOCATE(xt1(0:ntau-1,2),         error)
SLL_ALLOCATE(xt2(0:ntau-1,2),         error)

SLL_ALLOCATE(z(2),                    error)
SLL_ALLOCATE(fex(0:64,0:32,0:ntau-1), error)
SLL_ALLOCATE(fey(0:64,0:32,0:ntau-1), error)

SLL_CLEAR_ALLOCATE(f%ex(0:nx,0:ny),   error)
SLL_CLEAR_ALLOCATE(f%ey(0:nx,0:ny),   error)
SLL_CLEAR_ALLOCATE(f%bz(0:nx,0:ny),   error)
SLL_CLEAR_ALLOCATE(f%r0(0:nx,0:ny),   error)

time  = 0.d0

ep    = 0.1d0/1.0d0
epsq  = ep * ep
dtau  = 2.0d0*sll_p_pi/ntau
call sll_s_fft_init_c2c_1d(fw,Ntau,temp1,temp1,sll_p_fft_forward)
call sll_s_fft_init_c2c_1d(bw,Ntau,temp1,temp1,sll_p_fft_backward)
m=ntau/2
ltau=(/ (n, n=0,m-1), (n, n=-m,-1 )/)
pl(0)=dt
ql(0)=dt**2/2.0d0
do i=1,Ntau-1
  pl(i)=2.0d0*epsq*sll_p_i1*(exp(-sll_p_i1*ltau(i)*dt/2.0d0/epsq)-1.0d0)/ltau(i)
  ql(i)=2.0d0*epsq*(2.0d0*epsq*(1.0d0-exp(-sll_p_i1*ltau(i)*dt/2.0d0/epsq)) &
                    -sll_p_i1*ltau(i)*dt)/ltau(i)**2
enddo
do i=0,Ntau-1
  tau(i) =i*dtau
enddo

do i=0,nx
  aux1 = alpha/kx * sin(kx*i*dx)
  aux2 = alpha * cos(kx*i*dx)
  do j=0,ny
    f%ex(i,j) = aux1
    f%r0(i,j) = aux2+sin(ky*j*dy)
    f%ey(i,j) = -cos(ky*j*dy)/ky!0.0d0!
  enddo
enddo
      
xmin = 0.0_f64; xmax = dimx
ymin = 0.0_f64; ymax = dimy

call plasma( p )
print"('ep = ', g15.3)", ep
print"('nbpart = ', g15.3)", nbpart
print"('dt = ', g15.3)", dt
poisson => sll_f_new_poisson_2d_periodic(xmin,xmax,nx,ymin,ymax,ny)

wp1(:) =   2.0d0*((p%dpx+p%idx)*dx+ep*p%vpy)
wp2(:) =   2.0d0*((p%dpy+p%idy)*dy-ep*p%vpx)
wm1(:) = - 2.0d0*ep*p%vpy
wm2(:) =   2.0d0*ep*p%vpx

do n=0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))
  do m=1,nbpart
    utmp   = 0.5_f64*(cost*wp1(m)-sint*wp2(m)+cost*wm1(m)+sint*wm2(m))
    vtmp   = 0.5_f64*(sint*wp1(m)+cost*wp2(m)-sint*wm1(m)+cost*wm2(m))
    xxt(1) = real( cost*utmp+sint*vtmp)
    xxt(2) = real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo
  call interpol_eb_m6( f, p )
  Et(:,n,1)=p%epx !g(0,tau,w(0))
  Et(:,n,2)=p%epy
enddo

do m=1,npp

  temp1= 2.0d0*Et(m,:,2)
  temp2=-2.0d0*Et(m,:,1)!g_+
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)

  do n=1,Ntau-1
    temp1(n)=-sll_p_i1*temp1(n)/ltau(n)/Ntau
    temp2(n)=-sll_p_i1*temp2(n)/ltau(n)/Ntau
  enddo

  temp1(0)=0.0d0
  temp2(0)=0.0d0
  call sll_s_fft_exec_c2c_1d(bw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(bw, temp2, temp2)!AF+
  up(:,m,1)=wp1(m)+2.0d0*epsq*(temp1-temp1(0))
  up(:,m,2)=wp2(m)+2.0d0*epsq*(temp2-temp2(0))!1st ini data of U_+
  !---
  do n=0,ntau-1
    cost = cos(2.0d0*tau(n))
    sint = sin(2.0d0*tau(n))
    temp1(n)=-2.0d0*( sin(2.0d0*tau(n))*Et(m,n,1)+cos(2.0d0*tau(n))*Et(m,n,2))
    temp2(n)= 2.0d0*(-sin(2.0d0*tau(n))*Et(m,n,2)+cos(2.0d0*tau(n))*Et(m,n,1))!g_-
  enddo
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
  do n=1,Ntau-1
    temp1(n)=-sll_p_i1*temp1(n)/ltau(n)/Ntau
    temp2(n)=-sll_p_i1*temp2(n)/ltau(n)/Ntau
  enddo
  temp1(0)=0.0d0
  temp2(0)=0.0d0
  call sll_s_fft_exec_c2c_1d(bw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(bw, temp2, temp2)!AF-
  um(:,m,1)=wm1(m)+2.0d0*epsq*(temp1-temp1(0))
  um(:,m,2)=wm2(m)+2.0d0*epsq*(temp2-temp2(0))!1st ini data of U_-
enddo

!--corrected more initial data

do n=0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))

  do m=1,nbpart
    utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2)+cost*um(n,m,1)+sint*um(n,m,2))
    vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2)-sint*um(n,m,1)+cost*um(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
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
  do m=1,nbpart
    utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2)+cost*um(n,m,1)+sint*um(n,m,2))
    vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2)-sint*um(n,m,1)+cost*um(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo
  f%ex=fex(:,:,n)
  f%ey=fey(:,:,n)
  call interpol_eb_m6( f, p )
  Et(:,n,1)=p%epx !g_1st(0,tau,U_1st(0))
  Et(:,n,2)=p%epy
enddo

do m=1,npp
  temp1= 2.0d0*Et(m,:,2)
  temp2=-2.0d0*Et(m,:,1)!g_+
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
  up0(0,m,1)=temp1(0)/ntau!Pi g_+
  up0(0,m,2)=temp2(0)/ntau!Pi g_+
  do n=1,Ntau-1
    temp1(n)=-sll_p_i1*temp1(n)/ltau(n)/Ntau
    temp2(n)=-sll_p_i1*temp2(n)/ltau(n)/Ntau
  enddo
  temp1(0)=0.0d0
  temp2(0)=0.0d0
  call sll_s_fft_exec_c2c_1d(bw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(bw, temp2, temp2)!AF+
  up(:,m,1)=wp1(m)+2.0d0*epsq*(temp1-temp1(0))
  up(:,m,2)=wp2(m)+2.0d0*epsq*(temp2-temp2(0))!3rd ini data of U_+
  !---
  do n=0,ntau-1
    cost = cos(2d0*tau(n))
    sint = sin(2d0*tau(n))
    temp1(n)=-2.0d0*(sint*Et(m,n,1)+cost*Et(m,n,2))
    temp2(n)=2.0d0*(-sint*Et(m,n,2)+cost*Et(m,n,1))!g_-
  enddo
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
  um0(0,m,1)=temp1(0)/ntau!Pi g_-
  um0(0,m,2)=temp2(0)/ntau!Pi g_-
  do n=1,Ntau-1
    temp1(n)=-sll_p_i1*temp1(n)/ltau(n)/Ntau
    temp2(n)=-sll_p_i1*temp2(n)/ltau(n)/Ntau
  enddo
  temp1(0)=0.0d0
  temp2(0)=0.0d0
  call sll_s_fft_exec_c2c_1d(bw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(bw, temp2, temp2)!AF-
  um(:,m,1)=wm1(m)+2.0d0*epsq*(temp1-temp1(0))
  um(:,m,2)=wm2(m)+2.0d0*epsq*(temp2-temp2(0))!3rd ini data of U_-

enddo

do n=0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))
  do m=1,nbpart
    utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2)+cost*um(n,m,1)+sint*um(n,m,2))
    vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2)-sint*um(n,m,1)+cost*um(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
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
  do m=1,nbpart
    utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2)+cost*um(n,m,1)+sint*um(n,m,2))
    vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2)-sint*um(n,m,1)+cost*um(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo
  f%ex=fex(:,:,n)
  f%ey=fey(:,:,n)
  call interpol_eb_m6( f, p )
  Et(:,n,1)=p%epx !g_3rd(0,tau,U_3rd(0))
  Et(:,n,2)=p%epy
enddo

do m=1,npp
  temp1= 2.0d0*Et(m,:,2)
  temp2=-2.0d0*Et(m,:,1)
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
  xt1(:,1)=temp1/ntau!g_+tilde(t=0)
  xt1(:,2)=temp2/ntau!g_+tilde(t=0)
  !---
  do n=0,ntau-1
    cost = cos(2d0*tau(n))
    sint = sin(2d0*tau(n))
    temp1(n)=-2.0d0*(sint*Et(m,n,1)+cost*Et(m,n,2))
    temp2(n)=2.0d0*(-sint*Et(m,n,2)+cost*Et(m,n,1))
  enddo
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
  xt2(:,1)=temp1/ntau!g_-tilde(t=0)
  xt2(:,2)=temp2/ntau!g_-tilde(t=0)

  call sll_s_fft_exec_c2c_1d(fw, up(:,m,1), temp1)
  call sll_s_fft_exec_c2c_1d(fw, up(:,m,2), temp2)
  do n=0,ntau-1
    temp1(n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/epsq)*temp1(n)/ntau &
             +pl(n)*xt1(n,1)!utilde_+^1,predict
    temp2(n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/epsq)*temp2(n)/ntau &
             +pl(n)*xt1(n,2)!utilde_+^1,predict
  enddo
  call sll_s_fft_exec_c2c_1d(bw, temp1, up0(:,m,1))!u_+(t1),predict
  call sll_s_fft_exec_c2c_1d(bw, temp2, up0(:,m,2))
  call sll_s_fft_exec_c2c_1d(fw, um(:,m,1), temp1)
  call sll_s_fft_exec_c2c_1d(fw, um(:,m,2), temp2)
  do n=0,ntau-1
    temp1(n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/epsq)*temp1(n)/ntau &
            +pl(n)*xt2(n,1)!utilde_-^1,predict
    temp2(n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/epsq)*temp2(n)/ntau &
            +pl(n)*xt2(n,2)!utilde_-^1,predict
  enddo
  call sll_s_fft_exec_c2c_1d(bw, temp1, um0(:,m,1))!u_-(t1),predict
  call sll_s_fft_exec_c2c_1d(bw, temp2, um0(:,m,2))

  gp1(:,m) = xt1(:,1)
  gp2(:,m) = xt1(:,2)
  gm1(:,m) = xt2(:,1)
  gm2(:,m) = xt2(:,2)

enddo


do n= 0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))
  do m=1,nbpart
    utmp = 0.5_f64*(cost*up0(n,m,1)-sint*up0(n,m,2)+cost*um0(n,m,1)+sint*um0(n,m,2))
    vtmp = 0.5_f64*(sint*up0(n,m,1)+cost*up0(n,m,2)-sint*um0(n,m,1)+cost*um0(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo
  call calcul_rho_m6( p, f )
  call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
  fex(:,:,n)=cmplx(f%ex,0.0,f64)
  fey(:,:,n)=cmplx(f%ey,0.0,f64)!prediction
enddo
!--correction--
do n=0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))
  do m=1,nbpart
    utmp = 0.5_f64*(cost*up0(n,m,1)-sint*up0(n,m,2)+cost*um0(n,m,1)+sint*um0(n,m,2))
    vtmp = 0.5_f64*(sint*up0(n,m,1)+cost*up0(n,m,2)-sint*um0(n,m,1)+cost*um0(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
    call apply_bc()
    p%idx(m) = floor(xxt(1)/dimx*nx)
    p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
    p%idy(m) = floor(xxt(2)/dimy*ny)
    p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
  enddo
  f%ex=fex(:,:,n)
  f%ey=fey(:,:,n)
  call interpol_eb_m6( f, p )
  Et(:,n,1)=p%epx !g(t1,tau,U(t1))
  Et(:,n,2)=p%epy
enddo

do m=1,npp

  temp1 =  2.0d0*Et(m,:,2)
  temp2 = -2.0d0*Et(m,:,1)
  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
  xt1(:,1)=temp1/ntau!g_+tilde(t1) predict
  xt1(:,2)=temp2/ntau!g_+tilde(t1) predict
  !---
  do n=0,ntau-1
    cost = cos(2d0*tau(n))
    sint = sin(2d0*tau(n))
    temp1(n) = - 2.0d0*( sint*Et(m,n,1)+cost*Et(m,n,2))
    temp2(n) =   2.0d0*(-sint*Et(m,n,2)+cost*Et(m,n,1))
  enddo

  call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
  call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)

  xt2(:,1)=temp1/ntau!g_-tilde(t1) predict
  xt2(:,2)=temp2/ntau!g_-tilde(t1) predict

  call sll_s_fft_exec_c2c_1d(fw, up(:,m,1), temp1)
  call sll_s_fft_exec_c2c_1d(fw, up(:,m,2), temp2)

  do n=0,ntau-1
    temp1(n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/epsq)*temp1(n)/ntau &
             +pl(n)*xt1(n,1)+ql(n)*(xt1(n,1)-gp1(n,m))/dt
    temp2(n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/epsq)*temp2(n)/ntau &
             +pl(n)*xt1(n,2)+ql(n)*(xt1(n,2)-gp2(n,m))/dt
  enddo
  call sll_s_fft_exec_c2c_1d(bw, temp1, up(:,m,1))!u_+(t1)
  call sll_s_fft_exec_c2c_1d(bw, temp2, up(:,m,2))
  call sll_s_fft_exec_c2c_1d(fw, um(:,m,1), temp1)
  call sll_s_fft_exec_c2c_1d(fw, um(:,m,2), temp2)
  do n=0,ntau-1
    temp1(n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/epsq)*temp1(n)/ntau &
             +pl(n)*xt2(n,1)+ql(n)*(xt2(n,1)-gm1(n,m))/dt
    temp2(n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/epsq)*temp2(n)/ntau &
             +pl(n)*xt2(n,2)+ql(n)*(xt2(n,2)-gm2(n,m))/dt
  enddo
  call sll_s_fft_exec_c2c_1d(bw, temp1, um(:,m,1))!u_-(t1)
  call sll_s_fft_exec_c2c_1d(bw, temp2, um(:,m,2))

enddo

do n=0,ntau-1
  cost = cos(tau(n))
  sint = sin(tau(n))
  do m=1,nbpart
    utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2)+cost*um(n,m,1)+sint*um(n,m,2))
    vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2)-sint*um(n,m,1)+cost*um(n,m,2))
    xxt(1)=real( cost*utmp+sint*vtmp)
    xxt(2)=real(-sint*utmp+cost*vtmp)
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

!*** Loop over time ***

do istep = 2, nstep

  do n=0,ntau-1
    cost = cos(tau(n))
    sint = sin(tau(n))
    do m=1,nbpart
      utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2)+cost*um(n,m,1)+sint*um(n,m,2))
      vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2)-sint*um(n,m,1)+cost*um(n,m,2))
      xxt(1)=real( cost*utmp+sint*vtmp)
      xxt(2)=real(-sint*utmp+cost*vtmp)
      call apply_bc()
      p%idx(m) = floor(xxt(1)/dimx*nx)
      p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
      p%idy(m) = floor(xxt(2)/dimy*ny)
      p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
    enddo
    f%ex=fex(:,:,n)
    f%ey=fey(:,:,n)
    call interpol_eb_m6( f, p )
    Et(:,n,1)=p%epx 
    Et(:,n,2)=p%epy
  enddo

  do m=1,npp
    temp1=  2.0d0*Et(m,:,2)
    temp2= -2.0d0*Et(m,:,1)
    call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
    call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
    xt1(:,1)=temp1/ntau
    xt1(:,2)=temp2/ntau
    !---
    do n=0,ntau-1
      cost = cos(2d0*tau(n))
      sint = sin(2d0*tau(n))
      temp1(n) = -2.0d0*( sint*Et(m,n,1)+cost*Et(m,n,2))
      temp2(n) =  2.0d0*(-sint*Et(m,n,2)+cost*Et(m,n,1))
    enddo
    call sll_s_fft_exec_c2c_1d(fw, temp1, temp1)
    call sll_s_fft_exec_c2c_1d(fw, temp2, temp2)
    xt2(:,1)=temp1/ntau
    xt2(:,2)=temp2/ntau

    call sll_s_fft_exec_c2c_1d(fw, up(:,m,1), temp1)
    call sll_s_fft_exec_c2c_1d(fw, up(:,m,2), temp2)
    do n=0,ntau-1
      temp1(n)= exp(-sll_p_i1*ltau(n)*dt*0.5d0/epsq)*temp1(n)/ntau &
               + pl(n)*xt1(n,1)+ql(n)*(xt1(n,1)-gp1(n,m))/dt
      temp2(n)= exp(-sll_p_i1*ltau(n)*dt*0.5d0/epsq)*temp2(n)/ntau &
               + pl(n)*xt1(n,2)+ql(n)*(xt1(n,2)-gp2(n,m))/dt
    enddo
    call sll_s_fft_exec_c2c_1d(bw, temp1, up(:,m,1))
    call sll_s_fft_exec_c2c_1d(bw, temp2, up(:,m,2))
    call sll_s_fft_exec_c2c_1d(fw, um(:,m,1), temp1)
    call sll_s_fft_exec_c2c_1d(fw, um(:,m,2), temp2)
    do n=0,ntau-1
      temp1(n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/epsq)*temp1(n)/ntau &
               +pl(n)*xt2(n,1)+ql(n)*(xt2(n,1)-gm1(n,m))/dt
      temp2(n)=exp(-sll_p_i1*ltau(n)*dt*0.5d0/epsq)*temp2(n)/ntau &
               +pl(n)*xt2(n,2)+ql(n)*(xt2(n,2)-gm2(n,m))/dt
    enddo
    call sll_s_fft_exec_c2c_1d(bw, temp1, um(:,m,1))
    call sll_s_fft_exec_c2c_1d(bw, temp2, um(:,m,2))

    gp1(:,m)=xt1(:,1)
    gp2(:,m)=xt1(:,2)
    gm1(:,m)=xt2(:,1)
    gm2(:,m)=xt2(:,2)

  enddo


  !--updata E--
  time=dt*istep
  do n=0,ntau-1
    cost = cos(tau(n))
    sint = sin(tau(n))
    do m=1,nbpart
      utmp = 0.5_f64*(cost*up(n,m,1)-sint*up(n,m,2)+cost*um(n,m,1)+sint*um(n,m,2))
      vtmp = 0.5_f64*(sint*up(n,m,1)+cost*up(n,m,2)-sint*um(n,m,1)+cost*um(n,m,2))
      xxt(1)=real( cost*utmp+sint*vtmp)
      xxt(2)=real(-sint*utmp+cost*vtmp)
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

!*** Next time step ***

cost = cos(0.5_f64*time/epsq)
sint = sin(0.5_f64*time/epsq)

do m=1,npp
  call sll_s_fft_exec_c2c_1d(fw, up(:,m,1),temp1)
  call sll_s_fft_exec_c2c_1d(fw, up(:,m,2),temp2)
  wp1(m) = sll_p_i0
  wp2(m) = sll_p_i0
  do n=0,ntau-1
    wp1(m)=wp1(m)+temp1(n)/Ntau*exp(sll_p_i1*ltau(n)*time*0.5d0/epsq)
    wp2(m)=wp2(m)+temp2(n)/Ntau*exp(sll_p_i1*ltau(n)*time*0.5d0/epsq)
  enddo
  call sll_s_fft_exec_c2c_1d(fw, um(:,m,1),temp1)
  call sll_s_fft_exec_c2c_1d(fw, um(:,m,2),temp2)
  wm1(m) = sll_p_i0
  wm2(m) = sll_p_i0
  do n=0,ntau-1
    wm1(m)=wm1(m)+temp1(n)/Ntau*exp(sll_p_i1*ltau(n)*time*0.5d0/epsq)
    wm2(m)=wm2(m)+temp2(n)/Ntau*exp(sll_p_i1*ltau(n)*time*0.5d0/epsq)
  enddo
  utmp   = 0.5_f64*(cost*wp1(m)-sint*wp2(m)+cost*wm1(m)+sint*wm2(m))
  vtmp   = 0.5_f64*(cost*wp2(m)+sint*wp1(m)+cost*wm2(m)-sint*wm1(m))
  xxt(1) = real(cost*utmp+sint*vtmp)
  xxt(2) = real(cost*vtmp-sint*utmp)
  call apply_bc()
  p%idx(m) = floor(xxt(1)/dimx*nx)
  p%dpx(m) = real(xxt(1)/dx- p%idx(m), f64)
  p%idy(m) = floor(xxt(2)/dimy*ny)
  p%dpy(m) = real(xxt(2)/dy- p%idy(m), f64)
enddo

call calcul_rho_m6( p, f )
call sll_s_fft_free(fw)
call sll_s_fft_free(bw)

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

end program test_pic2d

