! conservation de la masse
!plot 'diag.dat' u 1:3 w lp title 'BSL', 'diag.dat' u 1:7 w lp title 'BSL NC', 'diag.dat' u 1:11 w lp title 'FSL', 'diag.dat' u 1:15 w lp title 'FSL NC'
! convergence en espace
!plot 'Conv_collela_rot_f3.dat' u 1:2 w lp title 'BSL', 'Conv_collela_rot_f3.dat' u 1:3 w lp title 'BSL NC', 'Conv_collela_rot_f3.dat' u 1:4 w lp title 'FSL', 'Conv_collela_rot_f3.dat' u 1:5 w lp title 'FSL NC'
! in rotating framework
program test_deposit_cubic_splines
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use sll_m_cubic_splines
use sll_m_constants
use sll_m_boundary_condition_descriptors
!use sll_fft
use, intrinsic :: iso_c_binding
implicit none

include "fftw3.f03"

type(C_PTR) :: PlnFwd,PlnBwd !,PlnFwdc
type(sll_cubic_spline_2D), pointer :: spl_fsl
sll_int32  :: Neta1,Neta2,step,nb_step
sll_int32  :: i,j,bc1_type,bc2_type,err
sll_real64 :: eta1,delta_eta1,eta1_min,eta1_max,eta2,delta_eta2,eta2_min,eta2_max
sll_real64 :: xi1_0,xi2_0
sll_real64 :: T,eps
sll_real64, dimension(:,:), pointer :: eta1feet,eta2feet
sll_int32, parameter ::  ntau =64,Nn=128!, Mn=514

sll_real64 :: x1_array(1:Nn+1),x2_array(1:Nn+1),r_array(1:Nn),v_array(1:Nn)
sll_real64 :: fh_fsl(1:Nn+1,1:Nn+1),f0(1:Nn,1:Nn)
sll_comp64 :: temp1(0:ntau-1),temp2(0:ntau-1), AF1(0:ntau-1), AF2(0:ntau-1),w1c(0:ntau-1),w2c(0:ntau-1)
sll_real64 :: k,h,gn,F1s(0:ntau-1),F2s(0:ntau-1),En(0:Ntau-1,1:Nn),Enr(0:Ntau-1,1:Nn),Ent(0:Ntau-1,1:Nn)
!sll_real64 :: f,tmpp, tmpm,nd,ff2,tmp
sll_comp64 :: sumup1,sumup2,dtgn
sll_real64 :: tau(0:ntau-1), ltau(0:ntau-1),lx(1:Nn),Etemp(1:Nn),w1_0(0:ntau-1),w2_0(0:ntau-1),fvr(1:Nn,1:Nn)!,ltaus(0:mn-1)
sll_comp64 :: F1(0:ntau-1),F2(0:ntau-1),Ftilde1(0:ntau-1),Ftilde2(0:ntau-1)
!sll_comp64 :: erftilde(0:mn-1),erf(1:Mn),tmp1(1:mn+1),tmp2(1:mn+1)
sll_comp64 :: Term2(0:ntau-1),Term1(0:ntau-1),dtF1_F(0:ntau-1),dtF2_F(0:ntau-1)
sll_comp64 :: temp1_F(0:ntau-1),temp2_F(0:ntau-1),dtF1(0:ntau-1),dtF2(0:ntau-1)
sll_int32  :: n,m

Neta1 = Nn
Neta2 = Nn

! mesh type : cartesian
! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
! BC        : periodic-periodic
eta1_min = -4._f64
eta1_max = 4._f64
eta2_min = -4._f64
eta2_max = 4._f64

PlnFwd = fftw_plan_dft_1d(ntau,Ftilde2,Ftilde1,FFTW_FORWARD, FFTW_MEASURE+FFTW_UNALIGNED)
PlnBwd = fftw_plan_dft_1d(ntau,Ftilde1,Ftilde2,FFTW_BACKWARD,FFTW_MEASURE+FFTW_UNALIGNED)

! ---- * Parameters * ----

! --- Space and time parameters --

! Final time
T =0.4_f64

! ---- * Construction of the mesh * ----
bc1_type = SLL_PERIODIC
bc2_type = SLL_PERIODIC

! ---- * Time and space steps * ----

! space steps
delta_eta1 = (eta1_max-eta1_min)/real(Neta1,f64)
delta_eta2 = (eta2_max-eta2_min)/real(Neta2,f64)

! time step and number of steps
k = 0.05_f64  !T*delta_eta1
nb_step = floor(T/k)
eps=0.001_f64
h=2.0_f64*sll_pi/ntau

! ---- * Messages * ----

print *,'# N=',Nn
print *,'# k=',k
print *,'# T=',T
print *,'# eps=',eps

! ---- * Allocation and creation of the splines * ----

! allocations of the arrays
SLL_ALLOCATE(eta1feet(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(eta2feet(Neta1+1,Neta2+1), err)

spl_fsl => new_cubic_spline_2D(Neta1+1, Neta2+1, &
eta1_min, eta1_max, &
eta2_min, eta2_max, &
bc1_type, bc2_type)

! ---- * Initializations * ----

! Analytic distribution function and data for the mesh
do i=1,Neta1+1
    eta1 = eta1_min + (i-1)*delta_eta1
    x1_array(i) = eta1
    do j=1,Neta2+1
        eta2 = eta2_min + (j-1)*delta_eta2
        x2_array(j) = eta2
        fh_fsl(i,j) = exp(-2.0_f64*eta2**2)*exp(-2.0_f64*eta1**2)
    enddo
enddo
r_array=x1_array(1:Nn)
v_array=x2_array(1:Nn)
do m=0,ntau-1
    tau(m)=dble(m)*h
enddo
m=ntau/2
ltau=(/ (n, n=0,m-1), (n, n=-m,-1 )/)
m=Nn/2
lx=(/ (n, n=0,m-1), (n, n=-m,-1 )/)*2.0_f64*sll_pi/(eta1_max-eta1_min)
t=0.0_f64
do step=1,nb_step ! ---- * Evolution in time * ----
    f0=fh_fsl(1:Nn,1:Nn)
    fh_fsl(:,Neta2+1)   = fh_fsl(:,1)
    call compute_cubic_spline_2D(fh_fsl,spl_fsl)
    call poissonsolver(f0,r_array,v_array,tau+t,lx,Nn,Ntau,En,Enr,Ent)
    do i=1,Neta1+1
    do j=1,Neta2+1
        xi1_0=x1_array(i)
        xi2_0=x2_array(j)
        do m=0,ntau-1
        Etemp=En(m,:)
        gn=g(tau(m)+t,xi1_0,xi2_0,Etemp,lx,Nn)
        F1(m)=-cos(2.0_f64*tau(m)+2.0_f64*t)**2*(0.5_f64*sin(2.0_f64*tau(m)+2.0_f64*t)*xi1_0+sin(tau(m)+t)**2*xi2_0)-sin(tau(m)+t)*gn  !g(tau,xi_1,xi_2) should be a function with the inputs
        F2(m)=cos(2.0_f64*tau(m)+2.0_f64*t)**2*(cos(tau(m)+t)**2*xi1_0+0.5_f64*sin(2.0_f64*tau(m)+2.0_f64*t)*xi2_0)+cos(tau(m)+t)*gn
        enddo
        call fftw_execute_dft(PlnFwd, F1, Ftilde1)
        call fftw_execute_dft(PlnFwd, F2, Ftilde2)
        Ftilde1=Ftilde1/dble(ntau)
        Ftilde2=Ftilde2/dble(ntau)
        do m=0,ntau-1
        sumup1=cmplx(0.0_f64,0.0_f64,kind=f64)
        sumup2=cmplx(0.0_f64,0.0_f64,kind=f64)
            do n=1,ntau-1
            sumup1=sumup1-sll_i1*Ftilde1(n)/ltau(n)*(exp(sll_i1*ltau(n)*tau(m))-1.0_f64)
            sumup2=sumup2-sll_i1*Ftilde2(n)/ltau(n)*(exp(sll_i1*ltau(n)*tau(m))-1.0_f64)
            enddo
        w1_0(m)=xi1_0+eps*dreal(sumup1)
        w2_0(m)=xi2_0+eps*dreal(sumup2)
        enddo
!------------for 2nd order correction------------------
        do m=0,ntau-1
        Etemp=En(m,:)
        gn=g(tau(m)+t,w1_0(m),w2_0(m),Etemp,lx,Nn)
        F1s(m)=-cos(2.0_f64*tau(m)+2.0_f64*t)**2*(0.5_f64*sin(2.0_f64*tau(m)+2.0_f64*t)*w1_0(m)+sin(tau(m)+t)**2*w2_0(m))-sin(tau(m)+t)*gn
        F2s(m)=cos(2.0_f64*tau(m)+2.0_f64*t)**2*(cos(tau(m)+t)**2*w1_0(m)+0.5_f64*sin(2.0_f64*tau(m)+2.0_f64*t)*w2_0(m))+cos(tau(m)+t)*gn
        Etemp=Ent(m,:)
        dtgn=g(tau(m)+t,w1_0(m),w2_0(m),Etemp,lx,Nn)
        Etemp=Enr(m,:)
        dtgn=dtgn+g(tau(m)+t,w1_0(m),w2_0(m),Etemp,lx,Nn)*(cos(tau(m)+t)*Ftilde1(0)+sin(tau(m)+t)*Ftilde2(0))
        dtF1(m)=-cos(2.0_f64*tau(m)+2.0_f64*t)**2*(0.5_f64*sin(2.0_f64*tau(m)+2.0_f64*t)*Ftilde1(0)+sin(tau(m)+t)**2*Ftilde2(0))-sin(tau(m)+t)*dtgn
        dtF2(m)=cos(2.0_f64*tau(m)+2.0_f64*t)**2*(cos(tau(m)+t)**2*Ftilde1(0)+0.5_f64*sin(2.0_f64*tau(m)+2.0_f64*t)*Ftilde2(0))+cos(tau(m)+t)*dtgn
        enddo
        call fftw_execute_dft(PlnFwd,dtF1,dtF1_F)
        call fftw_execute_dft(PlnFwd,dtF2,dtF2_F)
        dtF1_F=dtF1_F/dble(ntau)
        dtF2_F=dtF2_F/dble(ntau)
        do m=0,ntau-1
        temp1(m)=cmplx(0.0_f64,0.0_f64,kind=f64)
        temp2(m)=cmplx(0.0_f64,0.0_f64,kind=f64)
            do n=1,ntau-1
            temp1(m)=(exp(sll_i1*tau(m)*ltau(n))-1.0_f64)*(dtF1_F(n)/sll_i1/ltau(n))+temp1(m)
            temp2(m)=(exp(sll_i1*tau(m)*ltau(n))-1.0_f64)*(dtF2_F(n)/sll_i1/ltau(n))+temp2(m)
            enddo
        enddo
        call fftw_execute_dft(PlnFwd,temp1,temp1_F)
        call fftw_execute_dft(PlnFwd,temp2,temp2_F)
        temp1_F=temp1_F/dble(ntau)
        temp2_F=temp2_F/dble(ntau)
        AF1=(temp1-temp1_F(0))*eps
        AF2=(temp2-temp2_F(0))*eps
        Term1=F1s-AF1
        Term2=F2s-AF2
        call fftw_execute_dft(PlnFwd,Term1, temp1_F)
        call fftw_execute_dft(PlnFwd,Term2, temp2_F)
        temp1_F=temp1_F/dble(ntau)
        temp2_F=temp2_F/dble(ntau)
        do m=0,ntau-1
        temp1(m)=cmplx(0.0_f64,0.0_f64,kind=f64)
        temp2(m)=cmplx(0.0_f64,0.0_f64,kind=f64)
            do n=1,ntau-1
            temp1(m)=(exp(sll_i1*tau(m)*ltau(n))-1.0_f64)*(temp1_F(n)/sll_i1/ltau(n))+temp1(m)
            temp2(m)=(exp(sll_i1*tau(m)*ltau(n))-1.0_f64)*(temp2_F(n)/sll_i1/ltau(n))+temp2(m)
            enddo
        enddo
        w1_0=xi1_0+eps*dreal(temp1)
        w2_0=xi2_0+eps*dreal(temp2)
!-----------------time solver----------------
!---1st UA
!        do m=0,ntau-1
!        Etemp=En(m,:)
!        gn=g(tau(m)+t,w1_0(m),w2_0(m),Etemp,lx,Nn)
!        F1(m)=-cos(2.0_f64*tau(m)+2.0_f64*t)**2*(0.5_f64*sin(2.0_f64*tau(m)+2.0_f64*t)*w1_0(m)+sin(tau(m)+t)**2*w2_0(m))-sin(tau(m)+t)*gn
!        F2(m)=cos(2.0_f64*tau(m)+2.0_f64*t)**2*(cos(tau(m)+t)**2*w1_0(m)+0.5_f64*sin(2.0_f64*tau(m)+2.0_f64*t)*w2_0(m))+cos(tau(m)+t)*gn
!        enddo
!        temp1=w1_0+k*F1
!        temp2=w2_0+k*F2
!        call fftw_execute_dft(PlnFwd, temp1, Ftilde1)
!        call fftw_execute_dft(PlnFwd, temp2, Ftilde2)
!        do m=0,ntau-1
!        Ftilde1(m)=Ftilde1(m)/(1.0_f64+sll_i1*k*ltau(m)/eps)/dble(ntau)
!        Ftilde2(m)=Ftilde2(m)/(1.0_f64+sll_i1*k*ltau(m)/eps)/dble(ntau)
!        enddo
!        sumup1=cmplx(0.0_f64,0.0_f64,kind=f64)
!        sumup2=cmplx(0.0_f64,0.0_f64,kind=f64)
!        do n=0,ntau-1
!            sumup1=sumup1+Ftilde1(n)*exp(sll_i1*ltau(n)*k/eps)
!            sumup2=sumup2+Ftilde2(n)*exp(sll_i1*ltau(n)*k/eps)
!        enddo
!---2nd UA
        do m=0,ntau-1
            Etemp=En(m,:)
            gn=g(tau(m)+t,w1_0(m),w2_0(m),Etemp,lx,Nn)
            F1(m)=-cos(2.0_f64*tau(m)+2.0_f64*t)**2*(0.5_f64*sin(2.0_f64*tau(m)+2.0_f64*t)*w1_0(m)+sin(tau(m)+t)**2*w2_0(m))-sin(tau(m)+t)*gn
            F2(m)=cos(2.0_f64*tau(m)+2.0_f64*t)**2*(cos(tau(m)+t)**2*w1_0(m)+0.5_f64*sin(2.0_f64*tau(m)+2.0_f64*t)*w2_0(m))+cos(tau(m)+t)*gn
        enddo
        temp1=w1_0+k/2.0_f64*F1
        temp2=w2_0+k/2.0_f64*F2
        call fftw_execute_dft(PlnFwd, temp1, Ftilde1)
        call fftw_execute_dft(PlnFwd, temp2, Ftilde2)
        do m=0,ntau-1
            Ftilde1(m)=Ftilde1(m)/(1.0_f64+sll_i1*k/2.0_f64*ltau(m)/eps)
            Ftilde2(m)=Ftilde2(m)/(1.0_f64+sll_i1*k/2.0_f64*ltau(m)/eps)
        enddo
        call fftw_execute_dft(PlnBwd, Ftilde1,AF1)
        call fftw_execute_dft(PlnBwd, Ftilde2,AF2)
        AF1=AF1/dble(ntau) !w1_h
        AF2=AF2/dble(ntau) !w2_h
        do m=0,ntau-1
        Etemp=En(m,:)
        gn=g(tau(m)+t,dreal(AF1(m)),dreal(AF2(m)),Etemp,lx,Nn)
        F1(m)=-cos(2.0_f64*tau(m)+2.0_f64*t)**2*(0.5_f64*sin(2.0_f64*tau(m)+2.0_f64*t)*AF1(m)+sin(tau(m)+t)**2*AF2(m))-sin(tau(m)+t)*gn
        F2(m)=cos(2.0_f64*tau(m)+2.0_f64*t)**2*(cos(tau(m)+t)**2*AF1(m)+0.5_f64*sin(2.0_f64*tau(m)+2.0_f64*t)*AF2(m))+cos(tau(m)+t)*gn
        enddo
        call fftw_execute_dft(PlnFwd, F1,  Ftilde1)
        call fftw_execute_dft(PlnFwd, F2,  Ftilde2)
        w1c=w1_0
        w2c=w2_0
        call fftw_execute_dft(PlnFwd, w1c, F1)
        call fftw_execute_dft(PlnFwd, w2c, F2)
        do m=0,ntau-1
        temp1(m)=(F1(m)*(1.0_f64-sll_i1*k/eps/2.0_f64*ltau(m))+k*Ftilde1(m))/(1.0_f64+sll_i1*k/2.0_f64*ltau(m)/eps)
        temp2(m)=(F2(m)*(1.0_f64-sll_i1*k/eps/2.0_f64*ltau(m))+k*Ftilde2(m))/(1.0_f64+sll_i1*k/2.0_f64*ltau(m)/eps)
        enddo
        temp1=temp1/dble(ntau)
        temp2=temp2/dble(ntau)
        sumup1=cmplx(0.0_f64,0.0_f64,kind=f64)
        sumup2=cmplx(0.0_f64,0.0_f64,kind=f64)
        do n=0,ntau-1
        sumup1=sumup1+temp1(n)*exp(sll_i1*ltau(n)*k/eps)
        sumup2=sumup2+temp2(n)*exp(sll_i1*ltau(n)*k/eps)
        enddo
!---------------end time solve-------------------------
        eta1=dreal(sumup1)
        eta2=dreal(sumup2)
        call apply_bc()
        eta1feet(i,j)=eta1
        eta2feet(i,j)=eta2
    enddo
    enddo
    call deposit_value_2D(eta1feet,eta2feet,spl_fsl,fh_fsl)
    t=dble(step)*k/eps
enddo
f0=fh_fsl(1:Nn,1:Nn)
call fvrinterp(f0,t,r_array,v_array,Nn,lx,fvr)
call fftw_destroy_plan(PlnFwd)
call fftw_destroy_plan(PlnBwd)
print *,'# t=',t
open(unit=850,file='fh.dat')
do i=1,Neta1
do j=1,Neta2
write(850,*)  fvr(i,j)!fh_fsl(i,j)!, !eta2feet(i,j),eta1feet(i,j),,
enddo
enddo
close(850)

contains

subroutine apply_bc()
! --- Corrections on the BC ---
do while (eta1>eta1_max)
eta1 = eta1-(eta1_max-eta1_min)
enddo
do while (eta1<eta1_min)
eta1 = eta1+(eta1_max-eta1_min)
enddo
do while (eta2>eta2_max)
eta2 = eta2-(eta2_max-eta2_min)
enddo
do while (eta2<eta2_min)
eta2 = eta2+(eta2_max-eta2_min)
enddo
if (abs(eta1)<1.0e-12_f64) then
eta1=0.0_f64
endif
if (abs(eta2)<1.0e-12_f64) then
eta2=0.0_f64
endif
end subroutine apply_bc

end program


!subroutine fvrinterp(fh_fsl,t,r,v,Nn,lx,fvr)
!! ---interpolation-------
!use sll_m_constants
!use, intrinsic :: iso_c_binding
!implicit none
!include "fftw3.f03"
!integer(4), intent(in)   :: Nn
!real(8), intent(in)      :: fh_fsl(1:Nn,1:Nn),t,r(1:Nn),v(1:Nn),lx(1:Nn)
!real(8), intent(inout)   :: fvr(1:Nn,1:Nn)
!complex(8) :: ftemp(1:Nn,1:Nn),ftilde(1:Nn,1:Nn), sum0
!type(C_PTR) :: PlnF2
!real(8)    :: x,y,L
!integer(4) :: n,m,p,q
!L=4.0_f64
!PlnF2= fftw_plan_dft_2d(Nn,Nn,ftemp,ftilde,FFTW_FORWARD,FFTW_ESTIMATE+FFTW_UNALIGNED)
!ftemp=fh_fsl
!call  fftw_execute_dft(PlnF2, ftemp, ftilde)
!ftilde=ftilde/dble(Nn**2)
!do n=1,Nn
!    do m=1,Nn
!    x=cos(t)*r(n)-sin(t)*v(m)
!    y=sin(t)*r(n)+cos(t)*v(m)
!    sum0=cmplx(0.0_f64,0.0_f64,kind=f64)
!        if (abs(x)<L .and. abs(y)<L) then
!            do p=1,Nn
!                do q=1,Nn
!                sum0=sum0+exp(sll_i1*lx(p)*(x+L)+sll_i1*lx(q)*(y+L))*ftilde(p,q)
!                enddo
!            enddo
!        endif
!    fvr(n,m)=dreal(sum0)
!    enddo
!enddo
!call fftw_destroy_plan(PlnF2)
!end subroutine fvrinterp

subroutine fvrinterp(fh_fsl,t,r,v,Nn,lx,fvr)
!------interpolation by NUFFT-------
use sll_m_constants
use, intrinsic :: iso_c_binding
implicit none
include "fftw3.f03"

sll_int32,  intent(in)    :: Nn
sll_real64, intent(in)    :: fh_fsl(1:Nn,1:Nn),t,r(1:Nn),v(1:Nn),lx(1:Nn)
sll_real64, intent(inout) :: fvr(1:Nn,1:Nn)

sll_comp64 :: ftilde(1:Nn,1:Nn), f(1:Nn)
sll_comp64, allocatable :: ftemp(:,:)
sll_comp64, allocatable :: ftemp_1d(:)
type(C_PTR) :: PlnF2
sll_real64  :: x,y,x_array(Nn),y_array(Nn),L,epsnufft
sll_int32   :: n,m,p,ier,ntrace(Nn),mtrace(Nn)

epsnufft=1.0e-9_f64
L=4.0_f64
allocate(ftemp(1:Nn,1:Nn))
PlnF2= fftw_plan_dft_2d(Nn,Nn,ftemp,ftilde,FFTW_FORWARD,FFTW_ESTIMATE+FFTW_UNALIGNED)
ftemp=fh_fsl
call  fftw_execute_dft(PlnF2, ftemp, ftilde)
ftilde=ftilde/dble(Nn**2)
m=Nn/2

do n=1,Nn
  f=ftilde(n,:)
  ftilde(n,:)=f((/ (/(p,p=m+1,Nn)/), (/(p,p=1,m)/) /) )
enddo

do n=1,Nn
  f=ftilde(:,n)
  ftilde(:,n)=f((/ (/(p,p=m+1,Nn)/), (/(p,p=1,m)/) /) )
enddo

p=0
do n=1,Nn
    do m=1,Nn
    x=cos(t)*r(n)-sin(t)*v(m)
    y=sin(t)*r(n)+cos(t)*v(m)
        if (abs(x)<L .and. abs(y)<L) then
            p=p+1
            ntrace(p)=n
            mtrace(p)=m
            x_array(p)=x
            y_array(p)=y
        else
            fvr(n,m)=0.0_f64
        endif
    enddo
enddo

allocate( ftemp_1d(1:p) )
call nufft2d2f90(p,x_array(1:p),y_array(1:p),ftemp_1d,1,epsnufft, Nn,Nn,ftilde,ier)

do n=1,p
  fvr(ntrace(n),mtrace(n))=dreal(ftemp_1d(n))
enddo

call fftw_destroy_plan(PlnF2)
deallocate( ftemp_1d )

end subroutine fvrinterp


subroutine poissonsolver(fh_fsl,r,v,tau,lx,Nn,Ntau,En,Enr,Ent)
! ---PoissonSolver-------
use sll_m_constants
use, intrinsic :: iso_c_binding
implicit none
include "fftw3.f03"

sll_int32,  intent(in)    :: Nn,Ntau
sll_real64, intent(in)    :: fh_fsl(1:Nn,1:Nn),r(1:Nn),v(1:Nn),lx(1:Nn),tau(0:Ntau-1)
sll_real64, intent(inout) :: En(0:Ntau-1,1:Nn),Enr(0:Ntau-1,1:Nn),Ent(0:Ntau-1,1:Nn)

sll_real64  :: x(1:Nn),L,fvr(1:Nn,1:Nn),ftv(1:Nn,1:Nn),ftr(1:Nn,1:Nn)
sll_real64  :: ftemp1(1:Nn,1:Nn),ftemp2(1:Nn,1:Nn),xi1,xi2
sll_int32   :: n,m,i
sll_comp64  :: fvptilde(1:Nn),fvptilde0(1:Nn),temp(1:Nn),sum0(1:Nn)
sll_comp64  :: vctmp(Nn),uctmp(Nn)
type(C_PTR) :: PlnF,PlnB

L=4.0_f64
x=r
x(Nn/2+1)=1.0_f64
PlnF= fftw_plan_dft_1d(Nn,sum0,fvptilde,FFTW_FORWARD,FFTW_ESTIMATE+FFTW_UNALIGNED)
PlnB= fftw_plan_dft_1d(Nn,fvptilde,sum0,FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_UNALIGNED)

do i=0,Ntau-1
    call fvrinterp(fh_fsl,tau(i),r,v,Nn,lx,fvr)
    do n=1,Nn
        do m=1,Nn
            vctmp(m) = fvr(n,m)
        enddo
        call fftw_execute_dft(PlnF, vctmp, fvptilde)
        sum0(n)=fvptilde(1)/dble(Nn)*(2.0_f64*L)*r(n) !r*int_R fdv
    enddo
    call fftw_execute_dft(PlnF,sum0, fvptilde)
    do n=2,Nn
        fvptilde(n)=fvptilde(n)/sll_i1/lx(n)/dble(Nn)
    enddo
    fvptilde(1)=cmplx(0.0_f64,0.0_f64,kind=f64)
    call fftw_execute_dft(PlnB, fvptilde,temp)
    do n=1,Nn
        En(i,n)=dreal(temp(n)-temp(Nn/2+1))/x(n) !g(tau,r)
        Enr(i,n)=dreal(sum0(n)-En(i,n))/x(n)
    enddo
enddo

do n=1,Nn
    do m = 1,Nn
        vctmp(m) = fh_fsl(n,m)
        uctmp(m) = fh_fsl(m,n)
    enddo
    call fftw_execute_dft(PlnF, vctmp, fvptilde)
    call fftw_execute_dft(PlnF, uctmp, fvptilde0)
    do m=1,Nn
    fvptilde0(m)=fvptilde0(m)/dble(Nn)*sll_i1*lx(m)
    fvptilde(m)=fvptilde(m)/dble(Nn)*sll_i1*lx(m)
    enddo
    call fftw_execute_dft(PlnB, fvptilde, temp)
    ftv(n,:)=dreal(temp)  !\partial_\xi1 f_filde(\xi1,\xi2)
    call fftw_execute_dft(PlnB, fvptilde0, temp)
    ftr(:,n)=dreal(temp)  !\partial_\xi2 f_filde(\xi1,\xi2)
enddo

do i=0,Ntau-1
    call fvrinterp(ftv,tau(i),r,v,Nn,lx,ftemp1)
    call fvrinterp(ftr,tau(i),r,v,Nn,lx,ftemp2)
    do n=1,Nn
        do m=1,Nn
        uctmp=En(i,:)
        xi1=r(n)*cos(tau(i))-v(m)*sin(tau(i))
        xi2=r(n)*sin(tau(i))+v(m)*cos(tau(i))
        vctmp(m)=(cos(2.0_f64*tau(i))**2*(xi1*cos(tau(i))+xi2*sin(tau(i)))+g(tau(i),xi1,xi2,dreal(uctmp),lx,Nn))*(-sin(tau(i))*ftemp1(n,m)+cos(tau(i))*ftemp2(n,m))!partial_t f_tilde(xi1,xi2)
        enddo
        call fftw_execute_dft(PlnF, vctmp, fvptilde)
        sum0(n)=fvptilde(1)/dble(Nn)*(2.0_f64*L)*r(n)
    enddo
    call fftw_execute_dft(PlnF,sum0, fvptilde)
    do n=2,Nn
        fvptilde(n)=fvptilde(n)/sll_i1/lx(n)/dble(Nn)
    enddo
    fvptilde(1)=cmplx(0.0_f64,0.0_f64,kind=f64)
    call fftw_execute_dft(PlnB, fvptilde,temp)
    do n=1,Nn
        Ent(i,n)=dreal(temp(n)-temp(Nn/2+1))/x(n)!E_tilde(tau,r)
    enddo
enddo

call fftw_destroy_plan(PlnF)
call fftw_destroy_plan(PlnB)

contains

  !-----function g(tau,xi1,xi2,E,lx)---
double precision function g(ta,x1,x2,E,lx,Nn)

use sll_m_constants
use, intrinsic :: iso_c_binding
implicit none
include "fftw3.f03"

sll_int32,  intent(in) :: Nn
sll_real64, intent(in) :: x1,x2,ta,E(1:Nn),lx(1:Nn)

type(C_PTR) :: PlnF
sll_real64  :: x,L
sll_int32   :: n
sll_comp64  :: Etilde1(1:Nn),Etilde2(1:Nn),sum0, E0(1:Nn)

L=4.0_f64
PlnF = fftw_plan_dft_1d(Nn,Etilde2,Etilde1,FFTW_FORWARD, FFTW_ESTIMATE+FFTW_UNALIGNED)
E0=E
call fftw_execute_dft(PlnF,E0,Etilde1)
Etilde1=Etilde1/dble(Nn)
x=cos(ta)*x1+sin(ta)*x2
sum0=cmplx(0.0_f64,0.0_f64,kind=f64)

if (abs(x)<L) then
    do n=1,Nn
    sum0=sum0+exp(sll_i1*lx(n)*(x+L))*Etilde1(n)
    enddo
endif

g = dreal(sum0)
call fftw_destroy_plan(PlnF)

end function g


end subroutine poissonsolver


