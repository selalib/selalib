!Mar 1st 2016 two-scale solver
!2nd order scheme and initial data
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
type(sll_t_fft) :: PlnF, PlnB
integer(4)  ,parameter :: Ntau=16,npp=120000
real(8)  :: xt(2,0:Ntau-1),xu(2,npp),vu(2,npp),h(2,0:Ntau-1,npp),w(2,0:Ntau-1,npp)
real(8)  :: xxt(2),xuh(2,npp),vuh(2,npp),hh(2,0:Ntau-1,npp),wh(2,0:Ntau-1,npp)
complex(8) :: temp(2,0:Ntau-1),gn(2,0:Ntau-1,npp),gntilde(2,0:Ntau-1,npp)
real(8)  :: ep,dtau,tau(0:Ntau-1),ltau(0:Ntau-1),Et(2,0:Ntau-1,npp),aux(2,npp),auxpx(2,npp)
!real(8),  allocatable :: elec(:)
complex(8) :: g(2,0:Ntau-1),gtilde(2,0:Ntau-1),fn(2,0:Ntau-1),ftilde(2,0:Ntau-1),xtemp(2)
complex(8) :: temp5th(2,0:Ntau-1),temp5th1(2,0:Ntau-1),temp4th(2,0:Ntau-1),temp4th1(2,0:Ntau-1)
complex(8) ::  dtw(2,0:Ntau-1),dth(2,0:Ntau-1),dtx2(2)
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
iplot = 0

!---added---
ep=0.5d0/1.0d0
dtau=2.0d0*sll_p_pi/ntau
call sll_s_fft_init_c2c_1d(PlnF,Ntau,g(1,:),gtilde(1,:),sll_p_fft_FORWARD,optimization=sll_p_FFT_MEASURE)
call sll_s_fft_init_c2c_1d(PlnB,Ntau,g(1,:),gtilde(1,:),sll_p_FFT_BACKWARD,optimization=sll_p_FFT_MEASURE)
!PlnF= fftw_plan_dft_1d(Ntau,g(1,:),gtilde(1,:),FFTW_FORWARD,FFTW_MEASURE+FFTW_UNALIGNED)
!PlnB= fftw_plan_dft_1d(Ntau,g(1,:),gtilde(1,:),FFTW_BACKWARD,FFTW_MEASURE+FFTW_UNALIGNED)
m=ntau/2
ltau=(/ (n, n=0,m-1), (n, n=-m,-1 )/)
do i=0,Ntau-1
    tau(i) =i*dtau
enddo
!--end added--

!if( nstep > nstepmax ) nstep = nstepmax


!********************************************************************

istep = 1

do i=0,nx
  aux1 = alpha/kx * sin(kx*i*dx)
  aux2 = alpha * cos(kx*i*dx)
  do j=0,ny
    f%ex(i,j) = aux1
    f%r0(i,j) = aux2!+dsin(ky*j*dy)+1.0d0
    !f%ey(i,j) = -dcos(ky*j*dy)/ky
  enddo
enddo
      
xmin = 0.0_f64; xmax = dimx
ymin = 0.0_f64; ymax = dimy

call plasma( p ) 
print"('nbpart = ', g15.3)", nbpart
print"('dt = ', g15.3)", dt
call calcul_rho( p, f )

!gnuplot -p rho.gnu (to plot the initial rho)
!call sll_o_gnuplot_2d(xmin, xmax, nx+1, &
!                      ymin, ymax, ny+1, &
!                      f%r0, 'rho', 1, error)

poisson => sll_f_new_poisson_2d_periodic(xmin,xmax,nx,ymin,ymax,ny)
call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)

!gnuplot -p ex.gnu (to plot the initial ex field)
!call sll_o_gnuplot_2d(xmin, xmax, nx+1, &
!                      ymin, ymax, ny+1, &
!                      f%ex, 'ex', 1, error)
call interpol_eb( f, p )
!----------prepared initial data-----
!--3rd order--
!aux(1,:)=p%epx
!aux(2,:)=p%epy
!auxpx(1,:)=p%dpx+xmin+p%idx*dx
!auxpx(2,:)=p%dpy+ymin+p%idy*dy
!xu(1,:)=auxpx(1,:)+ep*p%vpy
!xu(2,:)=auxpx(2,:)-ep*p%vpx
!do n=0,Ntau-1
!    h(1,n,:)=ep*(dsin(tau(n))*p%vpx-dcos(tau(n))*p%vpy)
!    h(2,n,:)=ep*(dsin(tau(n))*p%vpy+dcos(tau(n))*p%vpx)
!enddo
!do n=0,Ntau-1
!    do m=1,nbpart
!        xt=xu(:,m)+h(:,n,m)
!        call apply_bc()
!        p%idx(m) = floor((xt(1)-xmin)/dimx*nx)
!        p%dpx(m) = real(xt(1)-xmin- p%idx(m)*dx, f32)
!        p%idy(m) = floor((xt(2)-ymin)/dimy*ny)
!        p%dpy(m) = real(xt(2)-ymin- p%idy(m)*dy, f32)
!    enddo
!    call interpol_eb( f, p )
!    Et(1,n,:)=p%epx
!    Et(2,n,:)=p%epy
!enddo
!do m=1,nbpart
!    do n=0,Ntau-1
!        g(1,n)=-dsin(tau(n))*Et(2,n,m)+dcos(tau(n))*Et(1,n,m)
!        g(2,n)=dsin(tau(n))*Et(1,n,m)+dcos(tau(n))*Et(2,n,m)
!    enddo
!    call sll_s_fft_exec_c2c_1d(PlnF, g(1,:), gtilde(1,:))
!    call sll_s_fft_exec_c2c_1d(PlnF, g(2,:), gtilde(2,:))
!    do n=1,Ntau-1
!        gtilde(:,n)=-sll_p_i1*gtilde(:,n)/ltau(n)/Ntau
!    enddo
!    gtilde(:,0)=0.0d0
!    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:), g(1,:))
!    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), g(2,:)) !g=L^-1(I-Pi)exp(-J*tau)E(X_1st,0)
!!    vu(1,m)=p%vpx(m)-ep*g(1,0)
!!    vu(2,m)=p%vpy(m)-ep*g(2,0)
!!    w(:,:,m)=ep*g
!    do n=0,Ntau-1
!        h(1,n,m)=h(1,n,m)-ep**2*(dsin(tau(n))*aux(2,m)+dcos(tau(n))*aux(1,m))
!        h(2,n,m)=h(2,n,m)-ep**2*(-dsin(tau(n))*aux(1,m)+dcos(tau(n))*aux(2,m))
!    enddo
!    xu(:,m)=auxpx(:,m)-h(:,0,m)
!    Et(1,:,m)=xu(1,m)+h(1,:,m)
!    Et(2,:,m)=xu(2,m)+h(2,:,m) !X_2nd
!    do n=0,Ntau-1
!        h(1,n,m)=ep*(dsin(tau(n))*p%vpx(m)-dcos(tau(n))*p%vpy(m))-ep**2*(dsin(tau(n))*g(1,0)-dcos(tau(n))*g(2,0))
!        h(2,n,m)=ep*(dsin(tau(n))*p%vpy(m)+dcos(tau(n))*p%vpx(m))-ep**2*(dsin(tau(n))*g(2,0)+dcos(tau(n))*g(1,0))
!    enddo
!    do n=0,Ntau-1
!        fn(1,n)=dsin(tau(n))*g(2,n)+dcos(tau(n))*g(1,n)
!        fn(2,n)=-dsin(tau(n))*g(1,n)+dcos(tau(n))*g(2,n)
!    enddo
!    call sll_s_fft_exec_c2c_1d(PlnF, fn(1,:), gtilde(1,:))
!    call sll_s_fft_exec_c2c_1d(PlnF, fn(2,:), gtilde(2,:))
!    do n=1,Ntau-1
!        gtilde(:,n)=-sll_p_i1*gtilde(:,n)/ltau(n)/Ntau
!    enddo
!    gtilde(:,0)=0.0d0
!    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:), fn(1,:))
!    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), fn(2,:)) !L^-1(I-Pi)exp(J*tau)g(tau)
!    h(:,:,m)=h(:,:,m)+ep**2*fn
!!    do n=0,Ntau-1
!!        g(1,n)=alpha*dcos(kx*auxpx(1,m))*dcos(tau(n))*(dsin(tau(n))*p%vpx(m)-dcos(tau(n))*p%vpy(m))
!!        g(2,n)=alpha*dcos(kx*auxpx(1,m))*dsin(tau(n))*(dsin(tau(n))*p%vpx(m)-dcos(tau(n))*p%vpy(m))!e^{-tauJ}dxE(x,0)Av_0
!!    enddo
!!    call sll_s_fft_exec_c2c_1d(PlnF, g(1,:), gtilde(1,:))
!!    call sll_s_fft_exec_c2c_1d(PlnF, g(2,:), gtilde(2,:))
!    do n=0,Ntau-1
!        h(1,n,m)=h(1,n,m)+ep**3*alpha*dcos(kx*auxpx(1,m))*(dsin(tau(n))*p%vpx(m)-dcos(tau(n))*p%vpy(m))/2.0d0
!        h(2,n,m)=h(2,n,m)+ep**3*alpha*dcos(kx*auxpx(1,m))*(dsin(tau(n))*p%vpy(m)+dcos(tau(n))*p%vpx(m))/2.0d0  ! oringial
!!        h(1,n,m)=h(1,n,m)+ep**3*(alpha*dcos(kx*auxpx(1,m))+dsin(auxpx(2,m)))*(dsin(tau(n))*p%vpx(m)-dcos(tau(n))*p%vpy(m))/2.0d0 !(5.38)
!!        h(2,n,m)=h(2,n,m)+ep**3*(alpha*dcos(kx*auxpx(1,m))+dsin(auxpx(2,m)))*(dsin(tau(n))*p%vpy(m)+dcos(tau(n))*p%vpx(m))/2.0d0
!    enddo
!    xu(:,m)=auxpx(:,m)-h(:,0,m)
!    do n=0,Ntau-1   ! here involves derivatives: dx_E
!        w(1,n,m)=ep**3*(alpha*dcos(tau(n))*dcos(kx*auxpx(1,m))*aux(2,m))
!        w(2,n,m)=ep**3*(alpha*dsin(tau(n))*dcos(kx*auxpx(1,m))*aux(2,m))
!    enddo
!enddo
!do n=0,Ntau-1
!    do m=1,nbpart
!        xt=Et(:,n,m)
!        call apply_bc()
!        p%idx(m) = floor((xt(1)-xmin)/dimx*nx)
!        p%dpx(m) = real(xt(1)-xmin- p%idx(m)*dx, f32)
!        p%idy(m) = floor((xt(2)-ymin)/dimy*ny)
!        p%dpy(m) = real(xt(2)-ymin- p%idy(m)*dy, f32)
!    enddo
!    call interpol_eb( f, p )
!    Et(1,n,:)=p%epx
!    Et(2,n,:)=p%epy
!enddo
!do m=1,nbpart
!    do n=0,Ntau-1
!        g(1,n)=-dsin(tau(n))*Et(2,n,m)+dcos(tau(n))*Et(1,n,m)
!        g(2,n)=dsin(tau(n))*Et(1,n,m)+dcos(tau(n))*Et(2,n,m)
!    enddo
!    call sll_s_fft_exec_c2c_1d(PlnF, g(1,:), gtilde(1,:))
!    call sll_s_fft_exec_c2c_1d(PlnF, g(2,:), gtilde(2,:))
!    do n=1,Ntau-1
!        gtilde(:,n)=-sll_p_i1*gtilde(:,n)/ltau(n)/Ntau
!    enddo
!    gtilde(:,0)=0.0d0
!    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:), g(1,:))
!    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), g(2,:)) !g_2nd
!    w(:,:,m)=w(:,:,m)+ep*g
!    vu(1,m)=p%vpx(m)-w(1,0,m)
!    vu(2,m)=p%vpy(m)-w(2,0,m)
!enddo
!!    !!call dtEsolver()!dtE(x_0,0)==0
!!--end 3rd order--
!--5th order--
aux(1,:)=p%epx
aux(2,:)=p%epy
auxpx(1,:)=p%dpx+p%idx*dx
auxpx(2,:)=p%dpy+p%idy*dy
xu(1,:)=auxpx(1,:)+ep*p%vpy
xu(2,:)=auxpx(2,:)-ep*p%vpx
do n=0,Ntau-1
h(1,n,:)=ep*(dsin(tau(n))*p%vpx-dcos(tau(n))*p%vpy)
h(2,n,:)=ep*(dsin(tau(n))*p%vpy+dcos(tau(n))*p%vpx)
enddo
!do n=0,Ntau-1
!    do m=1,nbpart
!        xt=xu(:,m)+h(:,n,m)
!        call apply_bc()
!        p%idx(m) = floor((xt(1)-xmin)/dimx*nx)
!        p%dpx(m) = real(xt(1)-xmin- p%idx(m)*dx, f32)
!        p%idy(m) = floor((xt(2)-ymin)/dimy*ny)
!        p%dpy(m) = real(xt(2)-ymin- p%idy(m)*dy, f32)
!    enddo
!    call interpol_eb( f, p )
!    Et(1,n,:)=p%epx
!    Et(2,n,:)=p%epy !E(X_1st,0)
!enddo
do m=1,nbpart
    do n=0,Ntau-1
        xt(:,n)=xu(:,m)+h(:,n,m)
        Et(1,n,m)=alpha/kx*dsin(kx*xt(1,n))
        Et(2,n,m)=0.0d0
        g(1,n)=-dsin(tau(n))*Et(2,n,m)+dcos(tau(n))*Et(1,n,m)
        g(2,n)=dsin(tau(n))*Et(1,n,m)+dcos(tau(n))*Et(2,n,m)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, g(1,:), gtilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, g(2,:), gtilde(2,:))
    do n=0,Ntau-1
        temp4th(1,n)=dsin(tau(n))*gtilde(1,0)/Ntau-dcos(tau(n))*gtilde(2,0)/Ntau
        temp4th(2,n)=dcos(tau(n))*gtilde(1,0)/Ntau+dsin(tau(n))*gtilde(2,0)/Ntau!APie^{-Jtau}E(X_1st,0)
    enddo
    do n=1,Ntau-1
        gtilde(:,n)=-sll_p_i1*gtilde(:,n)/ltau(n)/Ntau
    enddo
    gtilde(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:), g(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), g(2,:)) !g=L^-1(I-Pi)exp(-J*tau)E(X_1st,0)
    do n=0,Ntau-1
        temp4th1(1,n)=dsin(tau(n))*g(2,n)+dcos(tau(n))*g(1,n)
        temp4th1(2,n)=-dsin(tau(n))*g(1,n)+dcos(tau(n))*g(2,n)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, temp4th1(1,:), temp(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp4th1(2,:), temp(2,:))
    do n=0,Ntau-1
        temp4th(:,n)=temp4th(:,n)+temp(:,0)/Ntau
    enddo
    do n=0,Ntau-1
        temp4th1(1,n)=alpha*dcos(kx*xt(1,n))*temp4th(1,n)
        temp4th1(2,n)=0.0d0
        temp4th(1,n)=-dsin(tau(n))*temp4th1(2,n)+dcos(tau(n))*temp4th1(1,n)
        temp4th(2,n)=dsin(tau(n))*temp4th1(1,n)+dcos(tau(n))*temp4th1(2,n)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, temp4th(1,:), temp(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp4th(2,:), temp(2,:))
    do n=0,Ntau-1
        dth(1,n)=-ep**4*(dsin(tau(n))*temp(1,0)/Ntau-dcos(tau(n))*temp(2,0)/Ntau)
        dth(2,n)=-ep**4*(dcos(tau(n))*temp(1,0)/Ntau+dsin(tau(n))*temp(2,0)/Ntau)!O(eps^4) in h_5
    enddo
    do n=1,Ntau-1
        temp4th1(:,n)=-sll_p_i1*temp(:,n)/ltau(n)/Ntau
    enddo
    temp4th1(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, temp4th1(1,:), temp5th(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, temp4th1(2,:), temp5th(2,:))
    do n=0,Ntau-1
        temp5th1(1,n)=dsin(tau(n))*temp5th(2,n)+dcos(tau(n))*temp5th(1,n)
        temp5th1(2,n)=-dsin(tau(n))*temp5th(1,n)+dcos(tau(n))*temp5th(2,n)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, temp5th1(1,:), temp5th(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp5th1(2,:), temp5th(2,:))
    do n=1,Ntau-1
        temp5th1(:,n)=-temp5th(:,n)/(ltau(n))**2/Ntau
    enddo
    temp5th1(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, temp5th1(1,:), temp5th(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, temp5th1(2,:), temp5th(2,:))!O(eps^3) in h_5
    dth=dth-ep**4*temp5th
    do n=1,Ntau-1
        temp(:,n)=-temp(:,n)/(ltau(n))**2/Ntau
    enddo
    temp(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, temp(1,:), temp4th(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, temp(2,:), temp4th(2,:))! O(eps^3) in w_4
    do n=0,Ntau-1
        h(1,n,m)=h(1,n,m)-ep**2*(dsin(tau(n))*aux(2,m)+dcos(tau(n))*aux(1,m))
        h(2,n,m)=h(2,n,m)-ep**2*(-dsin(tau(n))*aux(1,m)+dcos(tau(n))*aux(2,m))
    enddo
    xu(:,m)=auxpx(:,m)-h(:,0,m)
    do n=0,Ntau-1
        xt(:,n)=xu(:,m)+h(:,n,m)!X_2nd
        Et(1,n,m)=alpha/kx*dsin(kx*xt(1,n))
        Et(2,n,m)=0.0d0
    enddo
   do n=0,Ntau-1
        h(1,n,m)=ep*(dsin(tau(n))*p%vpx(m)-dcos(tau(n))*p%vpy(m))-ep**2*(dsin(tau(n))*g(1,0)-dcos(tau(n))*g(2,0))
        h(2,n,m)=ep*(dsin(tau(n))*p%vpy(m)+dcos(tau(n))*p%vpx(m))-ep**2*(dsin(tau(n))*g(2,0)+dcos(tau(n))*g(1,0))
    enddo
    do n=0,Ntau-1
        fn(1,n)=dsin(tau(n))*g(2,n)+dcos(tau(n))*g(1,n)
        fn(2,n)=-dsin(tau(n))*g(1,n)+dcos(tau(n))*g(2,n)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, fn(1,:), gtilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, fn(2,:), gtilde(2,:))
    do n=1,Ntau-1
        gtilde(:,n)=-sll_p_i1*gtilde(:,n)/ltau(n)/Ntau
    enddo
    gtilde(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:), fn(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), fn(2,:)) !L^-1(I-Pi)exp(J*tau)g(tau)
    h(:,:,m)=h(:,:,m)+ep**2*fn
    do n=0,Ntau-1
        h(1,n,m)=h(1,n,m)+ep**3*alpha*dcos(kx*auxpx(1,m))*(dsin(tau(n))*p%vpx(m)-dcos(tau(n))*p%vpy(m))/2.0d0
        h(2,n,m)=h(2,n,m)+ep**3*alpha*dcos(kx*auxpx(1,m))*(dsin(tau(n))*p%vpy(m)+dcos(tau(n))*p%vpx(m))/2.0d0  ! h_3
    enddo
    xu(:,m)=auxpx(:,m)-h(:,0,m)
    do n=0,Ntau-1   ! here involves derivatives: dx_E
        w(1,n,m)=ep**3*(alpha*dcos(tau(n))*dcos(kx*auxpx(1,m))*aux(2,m))
        w(2,n,m)=ep**3*(alpha*dsin(tau(n))*dcos(kx*auxpx(1,m))*aux(2,m))
    enddo
!enddo
!do n=0,Ntau-1
!    do m=1,nbpart
!        xt=Et(:,n,m)
!        call apply_bc()
!        p%idx(m) = floor((xt(1)-xmin)/dimx*nx)
!        p%dpx(m) = real(xt(1)-xmin- p%idx(m)*dx, f32)
!        p%idy(m) = floor((xt(2)-ymin)/dimy*ny)
!        p%dpy(m) = real(xt(2)-ymin- p%idy(m)*dy, f32)
!    enddo
!    call interpol_eb( f, p )
!    Et(1,n,:)=p%epx !E(X_2nd,0)
!    Et(2,n,:)=p%epy
!enddo
!do m=1,nbpart
    do n=0,Ntau-1
        g(1,n)=-dsin(tau(n))*Et(2,n,m)+dcos(tau(n))*Et(1,n,m)
        g(2,n)=dsin(tau(n))*Et(1,n,m)+dcos(tau(n))*Et(2,n,m)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, g(1,:), gtilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, g(2,:), gtilde(2,:))
    do n=0,Ntau-1
        dtw(1,n)=(dsin(tau(n))*gtilde(1,0)-dcos(tau(n))*gtilde(2,0))/Ntau
        dtw(2,n)=(dsin(tau(n))*gtilde(2,0)+dcos(tau(n))*gtilde(1,0))/Ntau
        temp4th1(1,n)=(dsin(tau(n))*gtilde(2,0)+dcos(tau(n))*gtilde(1,0))/Ntau
        temp4th1(2,n)=(-dsin(tau(n))*gtilde(1,0)+dcos(tau(n))*gtilde(2,0))/Ntau
    enddo
    do n=1,Ntau-1
        gtilde(:,n)=-sll_p_i1*gtilde(:,n)/ltau(n)/Ntau
    enddo
    gtilde(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:), g(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), g(2,:)) !g_2nd
    w(:,:,m)=w(:,:,m)+ep*g  !w_3
    vu(1,m)=p%vpx(m)-w(1,0,m)
    vu(2,m)=p%vpy(m)-w(2,0,m)
    do n=0,Ntau-1
        temp5th(1,n)=dsin(tau(n))*g(2,n)+dcos(tau(n))*g(1,n)
        temp5th(2,n)=-dsin(tau(n))*g(1,n)+dcos(tau(n))*g(2,n)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, temp5th(1,:), temp5th1(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp5th(2,:), temp5th1(2,:))
    dtx2(1)=temp5th1(1,0)/Ntau! changed here!----
    dtx2(2)=temp5th1(2,0)/Ntau
    do n=0,Ntau-1
        temp5th(1,n)=alpha*dcos(kx*xt(1,n))*(dtw(1,n)+dtx2(1))! changed here!----
        temp5th(2,n)=0.0d0! changed here!----
    enddo
    do n=0,Ntau-1
        temp5th1(1,n)=-dsin(tau(n))*temp5th(2,n)+dcos(tau(n))*temp5th(1,n)
        temp5th1(2,n)=dsin(tau(n))*temp5th(1,n)+dcos(tau(n))*temp5th(2,n)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF,temp5th1(1,:), temp5th(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp5th1(2,:), temp5th(2,:))
    do n=1,Ntau-1
        temp5th(:,n)=-temp5th(:,n)/(ltau(n))**2/Ntau
    enddo
    temp5th(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB,temp5th(1,:), temp5th1(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, temp5th(2,:), temp5th1(2,:))
    dtw=-ep**3*temp5th1
!    dtx2(2)=-(2.0d0*E(2)-E(1)+0.0d0*E(1)/2.0d0)! need change here!----
!    dtx2(1)=-2.0d0*E(1)+E(2)+0.0d0*E(2)/2.0d0! need change here!----
!    temp5th(1,1)=2.0d0*dtx2(1)+dtx2(2)! need change here!----
!    temp5th(2,1)=2.0d0*dtx2(2)+dtx2(1)! need change here!----
    do n=0,Ntau-1
        temp5th1(1,n)=dcos(tau(n))*temp5th(2,1)+dsin(tau(n))*temp5th(1,1)
        temp5th1(2,n)=-dcos(tau(n))*temp5th(1,1)+dsin(tau(n))*temp5th(2,1)
    enddo
    dtw(:,:)=dtw(:,:)!-ep**5*temp5th1  !-ep^2L^-1d_tw
    do n=0,Ntau-1
        xt(:,n)=xu(:,m)+h(:,n,m) !X_3
    enddo
    h(:,:,m)=(ep**2)*temp4th1
    temp4th1=w(:,:,m)
    w(:,:,m)=-(ep**3)*temp4th
!enddo
!do n=0,Ntau-1
!    do m=1,nbpart
!        xt=Et(:,n,m)
!        call apply_bc()
!        p%idx(m) = floor((xt(1)-xmin)/dimx*nx)
!        p%dpx(m) = real(xt(1)-xmin- p%idx(m)*dx, f32)
!        p%idy(m) = floor((xt(2)-ymin)/dimy*ny)
!        p%dpy(m) = real(xt(2)-ymin- p%idy(m)*dy, f32)
!    enddo
!    call interpol_eb( f, p )
!    Et(1,n,:)=p%epx !E(X_3,0)
!    Et(2,n,:)=p%epy
!enddo
!do m=1,nbpart
    do n=0,Ntau-1
        Et(1,n,m)=alpha/kx*dsin(kx*xt(1,n)) !E(X_3,0)
        Et(2,n,m)=0.0d0
    enddo
    do n=0,Ntau-1
        h(1,n,m)=h(1,n,m)+ep*(-dcos(tau(n))*p%vpy(m)+dsin(tau(n))*p%vpx(m))
        h(2,n,m)=h(2,n,m)+ep*(dcos(tau(n))*p%vpx(m)+dsin(tau(n))*p%vpy(m))
    enddo
    do n=0,Ntau-1
        temp(1,n)=-dsin(tau(n))*Et(2,n,m)+dcos(tau(n))*Et(2,n,m)
        temp(2,n)=dsin(tau(n))*Et(1,n,m)+dcos(tau(n))*Et(2,n,m)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), gtilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), gtilde(2,:))
    do n=0,Ntau-1
        dth(1,n)=dth(1,n)+ep**2/Ntau*(dsin(tau(n))*gtilde(2,0)+dcos(tau(n))*gtilde(1,0))
        dth(2,n)=dth(2,n)+ep**2/Ntau*(dcos(tau(n))*gtilde(2,0)-dsin(tau(n))*gtilde(1,0))!O(eps^2) in h_5
    enddo
    do n=1,Ntau-1
        gtilde(:,n)=-sll_p_i1*gtilde(:,n)/ltau(n)/Ntau
    enddo
    gtilde(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:), g(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), g(2,:))
    w(:,:,m)=w(:,:,m)+ep*g
    vu(1,m)=p%vpx(m)-w(1,0,m)
    vu(2,m)=p%vpy(m)-w(2,0,m)
    do n=0,Ntau-1
        temp(1,n)=dsin(tau(n))*temp4th1(2,n)+dcos(tau(n))*temp4th1(1,n)
        temp(2,n)=-dsin(tau(n))*temp4th1(1,n)+dcos(tau(n))*temp4th1(2,n)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), gtilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), gtilde(2,:))
    do n=1,Ntau-1
        gtilde(:,n)=-sll_p_i1*gtilde(:,n)/ltau(n)/Ntau
    enddo
    gtilde(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:), g(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), g(2,:))
    h(:,:,m)=h(:,:,m)+ep*g
    xu(:,m)=auxpx(:,m)-h(:,0,m)
    do n=0,Ntau-1
        xt(:,n)=xu(:,m)+h(:,n,m) !X_4
    enddo
!enddo
!do n=0,Ntau-1
!    do m=1,nbpart
!        xt=Et(:,n,m)
!        call apply_bc()
!        p%idx(m) = floor((xt(1)-xmin)/dimx*nx)
!        p%dpx(m) = real(xt(1)-xmin- p%idx(m)*dx, f32)
!        p%idy(m) = floor((xt(2)-ymin)/dimy*ny)
!        p%dpy(m) = real(xt(2)-ymin- p%idy(m)*dy, f32)
!    enddo
!    call interpol_eb( f, p )
!    Et(1,n,:)=p%epx !E(X_4,0)
!    Et(2,n,:)=p%epy
!enddo
!do m=1,nbpart
    do n=0,Ntau-1
        Et(1,n,m)=alpha/kx*dsin(kx*xt(1,n)) !E(X_4,0)
        Et(2,n,m)=0.0d0
    enddo
    do n=0,Ntau-1
        temp(1,n)=-dsin(tau(n))*Et(2,n,m)+dcos(tau(n))*Et(1,n,m)
        temp(2,n)=dsin(tau(n))*Et(1,n,m)+dcos(tau(n))*Et(2,n,m)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, temp(1,:), gtilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp(2,:), gtilde(2,:))
    do n=1,Ntau-1
        gtilde(:,n)=-sll_p_i1*gtilde(:,n)/ltau(n)/Ntau
    enddo
    gtilde(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:), g(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), g(2,:))
    temp4th=w(:,:,m)   !w_4
    w(:,:,m)=ep*g+dtw !w_5
    do n=0,Ntau-1
        h(1,n,m)=(dsin(tau(n))*p%vpx(m)-dcos(tau(n))*p%vpy(m))*ep
        h(2,n,m)=(dsin(tau(n))*p%vpy(m)+dcos(tau(n))*p%vpx(m))*ep
    enddo
    do n=0,Ntau-1
        temp4th1(1,n)=dsin(tau(n))*temp4th(2,n)+dcos(tau(n))*temp4th(1,n)
        temp4th1(2,n)=-dsin(tau(n))*temp4th(1,n)+dcos(tau(n))*temp4th(2,n)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, temp4th1(1,:), gtilde(1,:))
    call sll_s_fft_exec_c2c_1d(PlnF, temp4th1(2,:), gtilde(2,:))
    do n=1,Ntau-1
        gtilde(:,n)=-sll_p_i1*gtilde(:,n)/ltau(n)/Ntau
    enddo
    gtilde(:,0)=0.0d0
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:), g(1,:))
    call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), g(2,:))
    h(:,:,m)=h(:,:,m)+ep*g+dth
    xu(:,m)=auxpx(:,m)-h(:,0,m)
    vu(1,m)=p%vpx(m)-w(1,0,m)
    vu(2,m)=p%vpy(m)-w(2,0,m)
enddo
!--end 5th order--
!-------end of preparation-----

!allocate (elec(nstep))
do istep = 1, nstep
!    !---ode solver 1st order scheme
!    time = time + dt
!    do m=1,nbpart
!        do n=0,Ntau-1
!            gn(1,n,m)=(dcos(tau(n))*w(1,n,m)+dsin(tau(n))*w(2,n,m))/ep
!            gn(2,n,m)=(dcos(tau(n))*w(2,n,m)-dsin(tau(n))*w(1,n,m))/ep
!        enddo
!        call sll_s_fft_exec_c2c_1d(PlnF, gn(1,:,m), gntilde(1,:,m))
!        call sll_s_fft_exec_c2c_1d(PlnF, gn(2,:,m), gntilde(2,:,m))
!        xu(:,m)=xu(:,m)+dt*dreal(gntilde(:,0,m))/Ntau!update \underline X
!    enddo
!    do n=0,Ntau-1
!        do m=1,nbpart
!            xt=xu(:,m)+h(:,n,m)
!            call apply_bc()
!            p%idx(m) = floor((xt(1)-xmin)/dimx*nx)
!            p%dpx(m) = real(xt(1)-xmin- p%idx(m)*dx, f32)
!            p%idy(m) = floor((xt(2)-ymin)/dimy*ny)
!            p%dpy(m) = real(xt(2)-ymin- p%idy(m)*dy, f32)
!        enddo
!        call interpol_eb( f, p )
!        Et(1,n,:)=p%epx
!        Et(2,n,:)=p%epy
!    enddo
!    do m=1,nbpart
!        do n=0,Ntau-1
!            fn(1,n)=dcos(tau(n))*Et(1,n,m)-dsin(tau(n))*Et(2,n,m)
!            fn(2,n)=dcos(tau(n))*Et(2,n,m)+dsin(tau(n))*Et(1,n,m)
!        enddo
!        call sll_s_fft_exec_c2c_1d(PlnF, fn(1,:), ftilde(1,:))
!        call sll_s_fft_exec_c2c_1d(PlnF, fn(2,:), ftilde(2,:))
!        vu(:,m)=vu(:,m)+dt/ep*dreal(ftilde(:,0))/Ntau
!        do n=0,Ntau-1
!            g(1,n)=gn(1,n,m)-gntilde(1,0,m)/Ntau+(dcos(tau(n))*vu(1,m)+dsin(tau(n))*vu(2,m))/ep
!            g(2,n)=gn(2,n,m)-gntilde(2,0,m)/Ntau+(dcos(tau(n))*vu(2,m)-dsin(tau(n))*vu(1,m))/ep
!            g(:,n)=g(:,n)*dt+h(:,n,m)
!        enddo
!        call sll_s_fft_exec_c2c_1d(PlnF, g(1,:), gtilde(1,:))
!        call sll_s_fft_exec_c2c_1d(PlnF, g(2,:), gtilde(2,:))
!        do n=0,Ntau-1
!            gtilde(:,n)=gtilde(:,n)/(1.0d0+sll_p_i1*ltau(n)*dt/ep**2)/Ntau
!        enddo
!        call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:),temp(1,:))
!        call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), temp(2,:))
!        h(:,:,m)=dreal(temp)
!        do n=0,Ntau-1
!            fn(:,n)=(fn(:,n)-ftilde(:,0)/Ntau)*ep*dt+w(:,n,m)*ep**2
!        enddo
!        call sll_s_fft_exec_c2c_1d(PlnF, fn(1,:), ftilde(1,:))
!        call sll_s_fft_exec_c2c_1d(PlnF, fn(2,:), ftilde(2,:))
!        do n=0,Ntau-1
!            ftilde(:,n)=ftilde(:,n)/(ep**2+sll_p_i1*ltau(n)*dt)/Ntau
!        enddo
!        call sll_s_fft_exec_c2c_1d(PlnB, ftilde(1,:), temp(1,:))
!        call sll_s_fft_exec_c2c_1d(PlnB, ftilde(2,:), temp(2,:))
!        w(:,:,m)=dreal(temp)
!    enddo
!    do m=1,nbpart
!        g(1,:)=xu(1,m)+h(1,:,m)
!        g(2,:)=xu(2,m)+h(2,:,m)
!        call sll_s_fft_exec_c2c_1d(PlnF, g(1,:),gtilde(1,:))
!        call sll_s_fft_exec_c2c_1d(PlnF, g(2,:),gtilde(2,:))
!        xtemp=0.0d0
!        do n=0,ntau-1
!            xtemp=xtemp+gtilde(:,n)/Ntau*cdexp(sll_p_i1*ltau(n)*time/ep**2)
!        enddo
!        xt=dreal(xtemp)
!        call apply_bc()
!        p%idx(m) = floor((xt(1)-xmin)/dimx*nx)
!        p%dpx(m) = real(xt(1) -xmin- p%idx(m)*dx, f32)
!        p%idy(m) = floor((xt(2)-ymin)/dimy*ny)
!        p%dpy(m) = real(xt(2) -ymin- p%idy(m)*dy, f32)
!    enddo
!    !---end 1st solver

    !---ode solver 2nd order scheme
    time = time + dt/2.0d0
    do m=1,nbpart
        do n=0,Ntau-1
            gn(1,n,m)=(dcos(tau(n))*w(1,n,m)+dsin(tau(n))*w(2,n,m))/ep
            gn(2,n,m)=(dcos(tau(n))*w(2,n,m)-dsin(tau(n))*w(1,n,m))/ep
        enddo
        call sll_s_fft_exec_c2c_1d(PlnF, gn(1,:,m), gntilde(1,:,m))
        call sll_s_fft_exec_c2c_1d(PlnF, gn(2,:,m), gntilde(2,:,m))
        xuh(:,m)=xu(:,m)+dt/2.0d0*dreal(gntilde(:,0,m))/Ntau! \underline X^n+1/2
    enddo
    do n=0,Ntau-1
        do m=1,nbpart
            xxt=xuh(:,m)+h(:,n,m)
            call apply_bc()
            p%idx(m) = floor((xxt(1)-xmin)/dimx*nx)
            p%dpx(m) = real(xxt(1)-xmin- p%idx(m)*dx, f32)
            p%idy(m) = floor((xxt(2)-ymin)/dimy*ny)
            p%dpy(m) = real(xxt(2)-ymin- p%idy(m)*dy, f32)
        enddo
        call interpol_eb( f, p )
        Et(1,n,:)=p%epx
        Et(2,n,:)=p%epy
    enddo
    do m=1,nbpart
        do n=0,Ntau-1
            fn(1,n)=dcos(tau(n))*Et(1,n,m)-dsin(tau(n))*Et(2,n,m)
            fn(2,n)=dcos(tau(n))*Et(2,n,m)+dsin(tau(n))*Et(1,n,m)
        enddo
        call sll_s_fft_exec_c2c_1d(PlnF, fn(1,:), ftilde(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, fn(2,:), ftilde(2,:))
        vuh(:,m)=vu(:,m)+dt/2.0d0/ep*dreal(ftilde(:,0))/Ntau ! \underline Y^n+1/2
        do n=0,Ntau-1
            g(1,n)=gn(1,n,m)-gntilde(1,0,m)/Ntau+(dcos(tau(n))*vuh(1,m)+dsin(tau(n))*vuh(2,m))/ep
            g(2,n)=gn(2,n,m)-gntilde(2,0,m)/Ntau+(dcos(tau(n))*vuh(2,m)-dsin(tau(n))*vuh(1,m))/ep
            g(:,n)=g(:,n)*dt/2.0d0+h(:,n,m)
        enddo
        call sll_s_fft_exec_c2c_1d(PlnF, g(1,:), gtilde(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, g(2,:), gtilde(2,:))
        do n=0,Ntau-1
            gtilde(:,n)=gtilde(:,n)/(1.0d0+sll_p_i1*ltau(n)*dt/2.0d0/ep**2)/Ntau
        enddo
        call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:),temp(1,:))
        call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), temp(2,:))
        hh(:,:,m)=dreal(temp)! h^n+1/2
        do n=0,Ntau-1
            fn(:,n)=(fn(:,n)-ftilde(:,0)/Ntau)*ep*dt/2.0d0+w(:,n,m)*ep**2
        enddo
        call sll_s_fft_exec_c2c_1d(PlnF, fn(1,:), ftilde(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, fn(2,:), ftilde(2,:))
        do n=0,Ntau-1
            ftilde(:,n)=ftilde(:,n)/(ep**2+sll_p_i1*ltau(n)*dt/2.0d0)/Ntau
        enddo
        call sll_s_fft_exec_c2c_1d(PlnB, ftilde(1,:), temp(1,:))
        call sll_s_fft_exec_c2c_1d(PlnB, ftilde(2,:), temp(2,:))
        wh(:,:,m)=dreal(temp)! w^n+1/2
    enddo
    do m=1,nbpart
        g(1,:)=xuh(1,m)+hh(1,:,m)
        g(2,:)=xuh(2,m)+hh(2,:,m)
        call sll_s_fft_exec_c2c_1d(PlnF, g(1,:),gtilde(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, g(2,:),gtilde(2,:))
        xtemp=0.0d0
        do n=0,ntau-1
            xtemp=xtemp+gtilde(:,n)/Ntau*cdexp(sll_p_i1*ltau(n)*time/ep**2)
        enddo
        xxt=dreal(xtemp)
        call apply_bc()
        p%idx(m) = floor((xxt(1)-xmin)/dimx*nx)
        p%dpx(m) = real(xxt(1) -xmin- p%idx(m)*dx, f32)
        p%idy(m) = floor((xxt(2)-ymin)/dimy*ny)
        p%dpy(m) = real(xxt(2) -ymin- p%idy(m)*dy, f32)
    enddo
    call calcul_rho( p, f )
    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0) !E^n+1/2
    do m=1,nbpart
        do n=0,Ntau-1
            gn(1,n,m)=(dcos(tau(n))*wh(1,n,m)+dsin(tau(n))*wh(2,n,m))/ep
            gn(2,n,m)=(dcos(tau(n))*wh(2,n,m)-dsin(tau(n))*wh(1,n,m))/ep
        enddo
        call sll_s_fft_exec_c2c_1d(PlnF, gn(1,:,m), gntilde(1,:,m))
        call sll_s_fft_exec_c2c_1d(PlnF, gn(2,:,m), gntilde(2,:,m))
        xuh(:,m)=xu(:,m)+dt*dreal(gntilde(:,0,m))/Ntau! \underline X^n+1
    enddo
    do n=0,Ntau-1
        do m=1,nbpart
            xxt=(xuh(:,m)+xu(:,m))/2.0d0+hh(:,n,m)
            call apply_bc()
            p%idx(m) = floor((xxt(1))/dimx*nx)
            p%dpx(m) = real(xxt(1)- p%idx(m)*dx, f32)
            p%idy(m) = floor((xxt(2))/dimy*ny)
            p%dpy(m) = real(xxt(2)- p%idy(m)*dy, f32)
        enddo
        call interpol_eb( f, p )
        Et(1,n,:)=p%epx
        Et(2,n,:)=p%epy
    enddo
    xu=xuh
    do m=1,nbpart
        do n=0,Ntau-1
            fn(1,n)=dcos(tau(n))*Et(1,n,m)-dsin(tau(n))*Et(2,n,m)
            fn(2,n)=dcos(tau(n))*Et(2,n,m)+dsin(tau(n))*Et(1,n,m)
        enddo
        call sll_s_fft_exec_c2c_1d(PlnF, fn(1,:), ftilde(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, fn(2,:), ftilde(2,:))
        vuh(:,m)=vu(:,m)+dt/ep*dreal(ftilde(:,0))/Ntau ! \underline Y^n+1
        aux(:,m)=vuh(:,m)
        vuh(:,m)=(vuh(:,m)+vu(:,m))/2.0d0
        vu(:,m)=aux(:,m)
        do n=0,Ntau-1
            g(1,n)=gn(1,n,m)-gntilde(1,0,m)/Ntau+(dcos(tau(n))*vuh(1,m)+dsin(tau(n))*vuh(2,m))/ep
            g(2,n)=gn(2,n,m)-gntilde(2,0,m)/Ntau+(dcos(tau(n))*vuh(2,m)-dsin(tau(n))*vuh(1,m))/ep
            g(:,n)=g(:,n)*dt
        enddo
        call sll_s_fft_exec_c2c_1d(PlnF, g(1,:), gtilde(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, g(2,:), gtilde(2,:))
        gn(:,:,1)=h(:,:,m)
        call sll_s_fft_exec_c2c_1d(PlnF, gn(1,:,1), temp(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, gn(2,:,1), temp(2,:))
        do n=0,Ntau-1
            gtilde(:,n)=(gtilde(:,n)+temp(:,n)*(1.0d0-sll_p_i1*ltau(n)*dt/2.0d0/ep**2)) &
            /(1.0d0+sll_p_i1*ltau(n)*dt/2.0d0/ep**2)/Ntau
        enddo
        call sll_s_fft_exec_c2c_1d(PlnB, gtilde(1,:),temp(1,:))
        call sll_s_fft_exec_c2c_1d(PlnB, gtilde(2,:), temp(2,:))
        h(:,:,m)=dreal(temp)! h^n+1
        do n=0,Ntau-1
            fn(:,n)=(fn(:,n)-ftilde(:,0)/Ntau)*ep*dt
        enddo
        call sll_s_fft_exec_c2c_1d(PlnF, fn(1,:), ftilde(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, fn(2,:), ftilde(2,:))
        gn(:,:,1)=w(:,:,m)
        call sll_s_fft_exec_c2c_1d(PlnF, gn(1,:,1), temp(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, gn(2,:,1), temp(2,:))
        do n=0,Ntau-1
            ftilde(:,n)=(ftilde(:,n)+temp(:,n)*(ep**2-sll_p_i1*ltau(n)*dt/2.0d0))   &
            /(ep**2+sll_p_i1*ltau(n)*dt/2.0d0)/Ntau
        enddo
        call sll_s_fft_exec_c2c_1d(PlnB, ftilde(1,:), temp(1,:))
        call sll_s_fft_exec_c2c_1d(PlnB, ftilde(2,:), temp(2,:))
        w(:,:,m)=dreal(temp)! w^n+1
    enddo
    time = time + dt/2.0d0
    do m=1,nbpart
        g(1,:)=xu(1,m)+h(1,:,m)
        g(2,:)=xu(2,m)+h(2,:,m)
        call sll_s_fft_exec_c2c_1d(PlnF, g(1,:),gtilde(1,:))
        call sll_s_fft_exec_c2c_1d(PlnF, g(2,:),gtilde(2,:))
        xtemp=0.0d0
        do n=0,ntau-1
            xtemp=xtemp+gtilde(:,n)/Ntau*cdexp(sll_p_i1*ltau(n)*time/ep**2)
        enddo
        xxt=dreal(xtemp)
        call apply_bc()
        p%idx(m) = floor((xxt(1)-xmin)/dimx*nx)
        p%dpx(m) = real(xxt(1) -xmin- p%idx(m)*dx, f32)
        p%idy(m) = floor((xxt(2)-ymin)/dimy*ny)
        p%dpy(m) = real(xxt(2) -ymin- p%idy(m)*dy, f32)
    enddo
    !---end 2nd solver
    call calcul_rho( p, f )
    call poisson%compute_e_from_rho( f%ex, f%ey, f%r0)
    !call maximumvalue()
    !call modeE( f, istep, time )
    !write(*,"('istep = ', i6, ' time = ',g15.3)", advance='no') istep, time
enddo
call sll_s_fft_free(PlnF)
call sll_s_fft_free(PlnB)
!call cpu_time(t1)
!print"('CPU time = ', g15.3)", t1-t0
!print*,'PASSED'

!open(unit=850,file='fh.dat')
!do i=1,nstep
!write(850,*)elec(i)!fh_fsl(i,j)!  , !eta2feet(i,j),eta1feet(i,j),
!enddo
!close(850)
open(unit=851,file='fh1.dat')
do i=1,nx
do j=1,ny
write(851,*)f%r0(i,j)!fh_fsl(i,j)!  , !eta2feet(i,j),eta1feet(i,j),
enddo
enddo
close(851)

!deallocate (elec)

contains

!subroutine maximumvalue()
!elec(istep)=0.0d0
!do i=1,nx
!    do j=1,ny
!        if (dabs(f%ex(i,j))>elec(istep)) then
!        elec(istep)=dabs(f%ex(i,j))
!        endif
!    enddo
!enddo
!end subroutine maximumvalue

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
