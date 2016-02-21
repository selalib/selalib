! Mouton example
! N_tau optimised
! taut optimised; 1d nufft corrected
! in rotating framework
program test_deposit_cubic_splines
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use sll_m_cubic_splines
use sll_m_constants
use sll_m_boundary_condition_descriptors
use sll_m_fft
implicit none

type(sll_t_fft) :: PlnFwd,PlnBwd,PlnF,PlnB

type(sll_t_cubic_spline_2D), pointer :: spl_fsl
sll_int32  :: step,nb_step
sll_int32  :: i,j,bc1_type,bc2_type,err
sll_real64 :: eta1,delta_eta1,eta1_min,eta1_max,eta2
sll_real64 :: xi1_0,xi2_0,x
sll_real64 :: T,eps!,L
sll_real64,dimension(:,:), pointer :: eta1feet,eta2feet
integer(4),parameter ::  ntau =32,Nn=128

sll_real64 :: x1_array(Nn+1),r_array(Nn),t1,t2
sll_real64 :: fh_fsl(Nn+1,Nn+1),f0(Nn,Nn)
complex(8) :: temp1(0:ntau-1),temp2(0:ntau-1), AF1(0:ntau-1), AF2(0:ntau-1),w1c(0:ntau-1),w2c(0:ntau-1)
real(8)    :: k,h,En(0:Ntau-1,1:Nn),Enr(0:Ntau-1,1:Nn),Ent(0:Ntau-1,1:Nn)
real(8)    :: gn(0:Ntau-1,Nn+1,Nn+1),gnr(0:Ntau-1,Nn+1,Nn+1),gnt(0:Ntau-1,Nn+1,Nn+1)
complex(8) :: sumup1,sumup2,dtgn
real(8)    :: tau(0:ntau-1), ltau(0:ntau-1),lx(1:Nn),fvr(1:Nn,1:Nn),taut(0:ntau-1)
real(8)    :: w1_0(0:ntau-1,Nn+1,Nn+1),w2_0(0:ntau-1,Nn+1,Nn+1)
complex(8) :: F1(0:ntau-1),F2(0:ntau-1),Ftilde1(0:ntau-1,Nn+1,Nn+1),Ftilde2(0:ntau-1,Nn+1,Nn+1)
complex(8) :: Term2(0:ntau-1),Term1(0:ntau-1)
complex(8) :: temp1_F(0:ntau-1),temp2_F(0:ntau-1),dtF1(0:ntau-1),dtF2(0:ntau-1)
integer(4) :: n,m
complex(8),  allocatable :: fgen1(:), fgen2(:)
!real(8),external :: seconds


! mesh type : cartesian
! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
! BC        : periodic-periodic
eta1_min = -4.0_f64
eta1_max =  4.0_f64
!L=4.0d0
allocate (fgen1(Nn))
allocate (fgen2(Nn))
call sll_s_fft_init_c2c_1d(PlnFwd ,ntau, AF2,   AF1,   sll_p_fft_forward)
call sll_s_fft_init_c2c_1d(PlnBwd ,ntau, AF1,   AF2,   sll_p_fft_forward)
call sll_s_fft_init_c2c_1d(PlnF   ,Nn,   Fgen1, Fgen2, sll_p_fft_forward)
call sll_s_fft_init_c2c_1d(PlnB   ,Nn,   Fgen2, Fgen1, sll_p_fft_forward)
deallocate (fgen1)
deallocate (fgen2)
! ---- * Parameters * ----

! --- Space and time parameters --

! Final time
T =0.4_f64

! ---- * Construction of the mesh * ----
bc1_type = SLL_P_PERIODIC
bc2_type = SLL_P_PERIODIC

! ---- * Time and space steps * ----

! space steps
delta_eta1 = (eta1_max-eta1_min)/real(Nn,f64)

! time step and number of steps
k =0.05d0  !T*delta_eta1
nb_step = floor(T/k)
eps=1.00d0
h=2.0d0*sll_p_pi/ntau

! ---- * Messages * ----

print *,'# N=',Nn
print *,'# k=',k
print *,'# T=',T
print *,'# eps=',eps

! ---- * Allocation and creation of the splines * ----

! allocations of the arrays
SLL_ALLOCATE(eta1feet(Nn+1,Nn+1), err)
SLL_ALLOCATE(eta2feet(Nn+1,Nn+1), err)
spl_fsl => sll_f_new_cubic_spline_2D(Nn+1, Nn+1, &
eta1_min, eta1_max, &
eta1_min, eta1_max, &
bc1_type, bc2_type)

! ---- * Initializations * ----

! Analytic distribution function and data for the mesh
do i=1,Nn+1
    eta1 = eta1_min + (i-1)*delta_eta1
    x1_array(i) = eta1
!    if (dabs(eta1)<=1.85) then
!        x=1.0d0
!    else
!        x=0.0d0
!    endif
    do j=1,Nn+1
        eta2 = eta1_min + (j-1)*delta_eta1
        fh_fsl(i,j) = dexp(-2.0d0*eta2**2)*dexp(-2.0d0*eta1**2)
!        fh_fsl(i,j) =4.0d0/dsqrt(0.4d0*sll_p_pi)*dexp(-eta2**2/0.4d0)*0.5d0*(derf((eta1+1.2d0)/0.3d0)-derf((eta1-1.2d0)/0.3d0))
        !fh_fsl(i,j) =4.0d0/dsqrt(2.0d0*sll_p_pi)/0.1d0*dexp(-eta2**2/2.0d0/(0.1d0)**2)*x
    enddo
enddo
do m=0,ntau-1
    tau(m)=dble(m)*h
enddo
m=ntau/2
ltau=[ real(8) :: (n, n=0,m-1), (n, n=-m,-1 ) ]
m=Nn/2
lx=(/ (n, n=0,m-1), (n, n=-m,-1 )/)*2.0d0*sll_p_pi/(eta1_max-eta1_min)
t=0.0d0
r_array=x1_array(1:Nn)
!t1= second();
!-------- * Evolution in time * ---------
do step=1,nb_step
    taut=tau+t
    f0=fh_fsl(1:Nn,1:Nn)
    call sll_s_compute_cubic_spline_2D(fh_fsl,spl_fsl)
    call poissonsolver(f0,r_array,taut,lx,Nn,Ntau,En,Enr,Ent,PlnF,PlnB)
    call ge0(Nn,Ntau,taut,x1_array,En,Ent,Enr,gn,gnt,gnr)  !gn=gn(tau,xi1,xi2) is a 3d array output

    do i=1,Nn+1
    do j=1,Nn+1
        xi1_0=x1_array(i)
        xi2_0=x1_array(j)
    !------------no correction------------------
!        do m=0,ntau-1
!            w1_0(m)=xi1_0
!            w2_0(m)=xi2_0
!        enddo
    !------------for 1st order correction------------------
        do m=0,ntau-1
        F1(m)=(-dcos(2.0d0*taut(m))**2*(0.5d0*dsin(2.0d0*taut(m))*xi1_0+dsin(taut(m))**2*xi2_0)-dsin(taut(m))*gn(m,i,j))/dble(ntau)
        F2(m)=(dcos(2.0d0*taut(m))**2*(dcos(taut(m))**2*xi1_0+0.5d0*dsin(2.0d0*taut(m))*xi2_0)+dcos(taut(m))*gn(m,i,j))/dble(ntau)
        enddo

        call sll_s_fft_exec_c2c_1d(PlnFwd, F1, AF1)
        call sll_s_fft_exec_c2c_1d(PlnFwd, F2, AF2)
        Ftilde1(:,i,j)=AF1
        Ftilde2(:,i,j)=AF2
        do m=1,ntau-1
            temp1(m)=-sll_p_i1*Ftilde1(m,i,j)/ltau(m)
            temp2(m)=-sll_p_i1*Ftilde2(m,i,j)/ltau(m)
        enddo
        temp1(0)=0.0d0
        temp2(0)=0.0d0
        call sll_s_fft_exec_c2c_1d(PlnBwd, temp1,F1)
        call sll_s_fft_exec_c2c_1d(PlnBwd, temp2,F2)
        F1=F1-sum(temp1)
        F2=F2-sum(temp2)
        w1_0(:,i,j)=xi1_0+eps*dreal(F1)
        w2_0(:,i,j)=xi2_0+eps*dreal(F2)
    enddo
    enddo
    call ge1(Nn,Ntau,taut,w1_0,w2_0,En,Ent,Enr,gn,gnt,gnr) ! gn,gnt,gnr 3d array output
    do i=1,Nn+1
    do j=1,Nn+1
    !------------for 2nd order correction------------------
        do m=0,ntau-1
        F1(m)=-dcos(2.0d0*taut(m))**2*(0.5d0*dsin(2.0d0*taut(m))*w1_0(m,i,j)+dsin(taut(m))**2*w2_0(m,i,j))-dsin(taut(m))*gn(m,i,j)
        F2(m)=dcos(2.0d0*taut(m))**2*(dcos(taut(m))**2*w1_0(m,i,j)+0.5d0*dsin(2.0d0*taut(m))*w2_0(m,i,j))+dcos(taut(m))*gn(m,i,j)
        dtgn=gnt(m,i,j)+gnr(m,i,j)*(dcos(taut(m))*Ftilde1(0,i,j)+dsin(taut(m))*Ftilde2(0,i,j))
        dtF1(m)=(-dcos(2.0d0*taut(m))**2*(0.5d0*dsin(2.0d0*taut(m))*Ftilde1(0,i,j)+dsin(taut(m))**2*Ftilde2(0,i,j))-dsin(taut(m))*dtgn)/dble(ntau)
        dtF2(m)=(dcos(2.0d0*taut(m))**2*(dcos(taut(m))**2*Ftilde1(0,i,j)+0.5d0*dsin(2.0d0*taut(m))*Ftilde2(0,i,j))+dcos(taut(m))*dtgn)/dble(ntau)
        enddo
        call sll_s_fft_exec_c2c_1d(PlnFwd,dtF1,AF1)
        call sll_s_fft_exec_c2c_1d(PlnFwd,dtF2,AF2)
        do m=1,ntau-1
            AF1(m)=-sll_p_i1*AF1(m)/ltau(m)
            AF2(m)=-sll_p_i1*AF2(m)/ltau(m)
        enddo
        AF1(0)=0.0d0
        AF2(0)=0.0d0
        call sll_s_fft_exec_c2c_1d(PlnBwd,AF1,dtF1)
        call sll_s_fft_exec_c2c_1d(PlnBwd,AF2,dtF2)
        temp1=(dtF1-sum(AF1))/dble(ntau)
        temp2=(dtF2-sum(AF2))/dble(ntau)
        call sll_s_fft_exec_c2c_1d(PlnFwd,temp1,temp1_F)
        call sll_s_fft_exec_c2c_1d(PlnFwd,temp2,temp2_F)
        Term1=(F1-(temp1-temp1_F(0))*eps)/dble(ntau)
        Term2=(F2-(temp2-temp2_F(0))*eps)/dble(ntau)
        call sll_s_fft_exec_c2c_1d(PlnFwd,Term1, temp1_F)
        call sll_s_fft_exec_c2c_1d(PlnFwd,Term2, temp2_F)
        do m=1,ntau-1
            temp1_F(m)=-sll_p_i1*temp1_F(m)/ltau(m)
            temp2_F(m)=-sll_p_i1*temp2_F(m)/ltau(m)
        enddo
        temp1_F(0)=0.0d0
        temp2_F(0)=0.0d0
        call sll_s_fft_exec_c2c_1d(PlnBwd, temp1_F,temp1)
        call sll_s_fft_exec_c2c_1d(PlnBwd, temp2_F,temp2)
        temp1=temp1-sum(temp1_F)
        temp2=temp2-sum(temp2_F)
        w1_0(:,i,j)=x1_array(i)+eps*dreal(temp1)
        w2_0(:,i,j)=x1_array(j)+eps*dreal(temp2)
    enddo
    enddo
    call ge2(Nn,Ntau,taut,w1_0,w2_0,En,gn)

    do i=1,Nn+1
    do j=1,Nn+1
!-----------------time solver----------------
!---1st UA
!        do m=0,ntau-1
!        Etemp=En(m,:)
!        gn=g(taut(m),w1_0(m),w2_0(m),Etemp,lx,Nn,PlnF)
!        F1(m)=-dcos(2.0d0*taut(m))**2*(0.5d0*dsin(2.0d0*taut(m))*w1_0(m)+dsin(taut(m))**2*w2_0(m))-dsin(taut(m))*gn
!        F2(m)=dcos(2.0d0*taut(m))**2*(dcos(taut(m))**2*w1_0(m)+0.5d0*dsin(2.0d0*taut(m))*w2_0(m))+dcos(taut(m))*gn
!        enddo
!        temp1=w1_0+k*F1
!        temp2=w2_0+k*F2
!        call fftw_execute_dft(PlnFwd, temp1, Ftilde1)
!        call fftw_execute_dft(PlnFwd, temp2, Ftilde2)
!        do m=0,ntau-1
!        Ftilde1(m)=Ftilde1(m)/(1.0d0+sll_p_i1*k*ltau(m)/eps)/dble(ntau)
!        Ftilde2(m)=Ftilde2(m)/(1.0d0+sll_p_i1*k*ltau(m)/eps)/dble(ntau)
!        enddo
!        sumup1=cmplx(0.0d0,0.0d0,kind=f64)
!        sumup2=cmplx(0.0d0,0.0d0,kind=f64)
!        do n=0,ntau-1
!            sumup1=sumup1+Ftilde1(n)*cdexp(sll_p_i1*ltau(n)*k/eps)
!            sumup2=sumup2+Ftilde2(n)*cdexp(sll_p_i1*ltau(n)*k/eps)
!        enddo
!---2nd UA
        do m=0,ntau-1
            F1(m)=-dcos(2.0d0*taut(m))**2*(0.5d0*dsin(2.0d0*taut(m))*w1_0(m,i,j)+dsin(taut(m))**2*w2_0(m,i,j))-dsin(taut(m))*gn(m,i,j)
            F2(m)=dcos(2.0d0*taut(m))**2*(dcos(taut(m))**2*w1_0(m,i,j)+0.5d0*dsin(2.0d0*taut(m))*w2_0(m,i,j))+dcos(taut(m))*gn(m,i,j)
        enddo
        temp1=w1_0(:,i,j)+k/2.0d0*F1
        temp2=w2_0(:,i,j)+k/2.0d0*F2
        call sll_s_fft_exec_c2c_1d(PlnFwd, temp1, AF1)
        call sll_s_fft_exec_c2c_1d(PlnFwd, temp2, AF2)
        do m=0,ntau-1
            AF1(m)=AF1(m)/(1.0d0+sll_p_i1*k/2.0d0*ltau(m)/eps)/dble(ntau)
            AF2(m)=AF2(m)/(1.0d0+sll_p_i1*k/2.0d0*ltau(m)/eps)/dble(ntau)
        enddo
        call sll_s_fft_exec_c2c_1d(PlnBwd, AF1,temp1_F)
        call sll_s_fft_exec_c2c_1d(PlnBwd, AF2,temp2_F)
        Ftilde1(:,i,j)=temp1_F
        Ftilde2(:,i,j)=temp2_F
        !----------Insert half step evaluation--------
        sumup1=cmplx(0.0d0,0.0d0,kind=f64)
        sumup2=cmplx(0.0d0,0.0d0,kind=f64)
        do n=0,ntau-1
            sumup1=sumup1+AF1(n)*cdexp(sll_p_i1*ltau(n)*k/2.0d0/eps)
            sumup2=sumup2+AF2(n)*cdexp(sll_p_i1*ltau(n)*k/2.0d0/eps)
        enddo
        eta1=dreal(sumup1)
        eta2=dreal(sumup2)
        call apply_bc()
        eta1feet(i,j)=eta1
        eta2feet(i,j)=eta2
    enddo
    enddo
    call sll_s_deposit_value_2D(eta1feet,eta2feet,spl_fsl,fh_fsl) !function value at the half time
    f0=fh_fsl(1:Nn,1:Nn)
    call poissonsolver2(f0,r_array,taut,lx,Nn,Ntau,En,PlnF,PlnB)
    !------End evaluation and continue 2nd solver-------------
    call ge2(Nn,Ntau,taut,dreal(Ftilde1),dreal(Ftilde2),En,gn)
    do i=1,Nn+1
    do j=1,Nn+1
        do m=0,ntau-1
        F1(m)=-dcos(2.0d0*taut(m))**2*(0.5d0*dsin(2.0d0*taut(m))*Ftilde1(m,i,j)+dsin(taut(m))**2*Ftilde2(m,i,j))-dsin(taut(m))*gn(m,i,j)
        F2(m)=dcos(2.0d0*taut(m))**2*(dcos(taut(m))**2*Ftilde1(m,i,j)+0.5d0*dsin(2.0d0*taut(m))*Ftilde2(m,i,j))+dcos(taut(m))*gn(m,i,j)
        enddo
        call sll_s_fft_exec_c2c_1d(PlnFwd, F1,  AF1)
        call sll_s_fft_exec_c2c_1d(PlnFwd, F2,  AF2)
        w1c=w1_0(:,i,j)
        w2c=w2_0(:,i,j)
        call sll_s_fft_exec_c2c_1d(PlnFwd, w1c, F1)
        call sll_s_fft_exec_c2c_1d(PlnFwd, w2c, F2)
        do m=0,ntau-1
        temp1(m)=(F1(m)*(1.0d0-sll_p_i1*k/eps/2.0d0*ltau(m))+k*AF1(m))/(1.0d0+sll_p_i1*k/2.0d0*ltau(m)/eps)/dble(ntau)
        temp2(m)=(F2(m)*(1.0d0-sll_p_i1*k/eps/2.0d0*ltau(m))+k*AF2(m))/(1.0d0+sll_p_i1*k/2.0d0*ltau(m)/eps)/dble(ntau)
        enddo
        sumup1=cmplx(0.0d0,0.0d0,kind=f64)
        sumup2=cmplx(0.0d0,0.0d0,kind=f64)
        do n=0,ntau-1
        sumup1=sumup1+temp1(n)*cdexp(sll_p_i1*ltau(n)*k/eps)
        sumup2=sumup2+temp2(n)*cdexp(sll_p_i1*ltau(n)*k/eps)
        enddo
!---------------end time solve-------------------------
        eta1=dreal(sumup1)
        eta2=dreal(sumup2)
        call apply_bc()
        eta1feet(i,j)=eta1
        eta2feet(i,j)=eta2
    enddo
    enddo
    call sll_s_deposit_value_2D(eta1feet,eta2feet,spl_fsl,fh_fsl)
    t=dble(step)*k/eps
enddo
!t2= second();
!print *, 'whole time = ' , t2 - t1
f0=fh_fsl(1:Nn,1:Nn)
call fvrinterp(f0,t,r_array,Nn,fvr)
call sll_s_fft_free(PlnFwd)
call sll_s_fft_free(PlnBwd)
call sll_s_fft_free(PlnF)
call sll_s_fft_free(PlnB)
open(unit=850,file='fh.dat')
do i=1,Nn
do j=1,Nn
write(850,*)i,j,sngl(fvr(i,j))!fh_fsl(i,j)!  , !eta2feet(i,j),eta1feet(i,j),
enddo
write(850,*)
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
do while (eta2>eta1_max)
eta2 = eta2-(eta1_max-eta1_min)
enddo
do while (eta2<eta1_min)
eta2 = eta2+(eta1_max-eta1_min)
enddo
if (dabs(eta1)<1.0d-12) then
eta1=0.0d0
endif
if (dabs(eta2)<1.0d-12) then
eta2=0.0d0
endif
end subroutine apply_bc

end program

subroutine fvrinterp(fh_fsl,t,r,Nn,fvr)
!------interpolation by cubic spline-------
use sll_m_boundary_condition_descriptors
use sll_m_cubic_splines
use sll_m_constants
implicit none
type(sll_t_cubic_spline_2D), pointer :: spl_fsl1
integer(4), intent(in)   :: Nn
real(8), intent(in)      :: fh_fsl(Nn,Nn),t,r(Nn)!,lx(1:Nn)
real(8), intent(inout)   :: fvr(1:Nn,1:Nn)
real(8)    :: x,y,f0(Nn+1,Nn+1),L
integer(4) :: n,m
L=4.0d0
spl_fsl1 => sll_f_new_cubic_spline_2D(Nn+1, Nn+1, &
-L, L, &
-L, L, &
SLL_P_PERIODIC, SLL_P_PERIODIC)
f0(1:Nn,1:Nn)=fh_fsl
do n=1,Nn
f0(n,Nn+1)=0.0d0
f0(Nn+1,n)=0.0d0
enddo

call sll_s_compute_cubic_spline_2D(f0,spl_fsl1)

do n=1,Nn
    do m=1,Nn
        x=dcos(t)*r(n)-dsin(t)*r(m)
        y=dsin(t)*r(n)+dcos(t)*r(m)
        if (dabs(x)<L .and. dabs(y)<L) then
            fvr(n,m)=sll_f_interpolate_value_2D(x,y,spl_fsl1)
        else
        fvr(n,m)=0.0d0
        endif
    enddo
enddo
call sll_o_delete(spl_fsl1)
end subroutine fvrinterp


subroutine poissonsolver(fh_fsl,r,tau,lx,Nn,Ntau,En,Enr,Ent,PlnF,PlnB)
! ---PoissonSolver-------
use sll_m_constants
use sll_m_fft
use sll_m_working_precision
implicit none
type(sll_t_fft), intent(in) :: PlnF,PlnB
integer(4), intent(in)   :: Nn,Ntau
real(8), intent(in)      :: fh_fsl(Nn,Nn),r(Nn),lx(1:Nn),tau(0:Ntau-1)
real(8), intent(inout)   :: En(0:Ntau-1,1:Nn),Enr(0:Ntau-1,1:Nn),Ent(0:Ntau-1,1:Nn)
real(8)    :: x(1:Nn),L,fvr(Nn,Nn),ftv(Nn,Nn),ftr(Nn,Nn)
real(8)    :: ftemp1(Nn,Nn),ftemp2(Nn,Nn),xi1(0:Ntau-1,Nn+1,Nn+1),xi2(0:Ntau-1,Nn+1,Nn+1),v(Nn+1)
integer(4) :: n,m,i
complex(8) :: fvptilde(1:Nn),fvptilde0(Nn),temp(Nn),sum0(Nn)
complex(8) :: vctmp(Nn),uctmp(Nn)
real(8) :: gn(0:Ntau-1,Nn+1,Nn+1)
L=4.0d0
x=r
x(Nn/2+1)=1.0d0
do i=0,Ntau-1
    call fvrinterp(fh_fsl,tau(i),r,Nn,fvr)
    do n=1,Nn
        vctmp = fvr(n,:)
        call sll_s_fft_exec_c2c_1d(PlnF, vctmp, fvptilde)
        sum0(n)=fvptilde(1)/dble(Nn)*(2.0d0*L)*r(n) !r*int_R fdv
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF,sum0, fvptilde)
    do n=2,Nn
        fvptilde(n)=fvptilde(n)/sll_p_i1/lx(n)/dble(Nn)
    enddo
    fvptilde(1)=cmplx(0.0d0,0.0d0,kind=f64)
    call sll_s_fft_exec_c2c_1d(PlnB, fvptilde,temp)
    do n=1,Nn
        En(i,n)=dreal(temp(n)-temp(Nn/2+1))/x(n) !g(tau,r)
        Enr(i,n)=dreal(sum0(n)-En(i,n))/x(n)
    enddo
enddo

do n=1,Nn
    vctmp = fh_fsl(n,:)
    uctmp = fh_fsl(:,n)
    call sll_s_fft_exec_c2c_1d(PlnF, vctmp, fvptilde)
    call sll_s_fft_exec_c2c_1d(PlnF, uctmp, fvptilde0)
    do m=1,Nn
    fvptilde0(m)=fvptilde0(m)/dble(Nn)*sll_p_i1*lx(m)
    fvptilde(m)=fvptilde(m)/dble(Nn)*sll_p_i1*lx(m)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnB, fvptilde, temp)
    ftv(n,:)=dreal(temp)  !\partial_\xi1 f_filde(\xi1,\xi2)
    call sll_s_fft_exec_c2c_1d(PlnB, fvptilde0, temp)
    ftr(:,n)=dreal(temp)  !\partial_\xi2 f_filde(\xi1,\xi2)
enddo

v(1:Nn)=r
v(Nn+1)=L
do i=0,Ntau-1
    do n=1,Nn+1
    do m=1,Nn+1
        xi1(i,n,m)=v(n)*dcos(tau(i))-v(m)*dsin(tau(i))
        xi2(i,n,m)=v(n)*dsin(tau(i))+v(m)*dcos(tau(i))
    enddo
    enddo
enddo

call ge2(Nn,Ntau,tau,xi1,xi2,En,gn)

do i=0,Ntau-1
    call fvrinterp(ftv,tau(i),r,Nn,ftemp1)
    call fvrinterp(ftr,tau(i),r,Nn,ftemp2)
    do n=1,Nn
        do m=1,Nn
        vctmp(m)=(dcos(2.0d0*tau(i))**2*(xi1(i,n,m)*dcos(tau(i))+xi2(i,n,m)*dsin(tau(i)))+gn(i,n,m))*(-dsin(tau(i))*ftemp1(n,m)+dcos(tau(i))*ftemp2(n,m))!partial_t f_tilde(xi1,xi2)
        enddo
        call sll_s_fft_exec_c2c_1d(PlnF, vctmp, fvptilde)
        sum0(n)=fvptilde(1)/dble(Nn)*(2.0d0*L)*r(n)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF,sum0, fvptilde)
    do n=2,Nn
        fvptilde(n)=fvptilde(n)/sll_p_i1/lx(n)/dble(Nn)
    enddo
    fvptilde(1)=cmplx(0.0d0,0.0d0,kind=f64)
    call sll_s_fft_exec_c2c_1d(PlnB, fvptilde,temp)
    do n=1,Nn
        Ent(i,n)=dreal(temp(n)-temp(Nn/2+1))/x(n)!E_tilde(tau,r)
    enddo
enddo
!call fftw_destroy_plan(PlnF)
!call fftw_destroy_plan(PlnB)
end subroutine poissonsolver

subroutine poissonsolver2(fh_fsl,r,tau,lx,Nn,Ntau,En,PlnF,PlnB)
! ---PoissonSolver-------
use sll_m_constants
use sll_m_working_precision
use sll_m_fft
implicit none
type(sll_t_fft), intent(in) :: PlnF,PlnB
integer(4), intent(in)   :: Nn,Ntau
real(8), intent(in)      :: fh_fsl(Nn,Nn),r(Nn),lx(Nn),tau(0:Ntau-1)
real(8), intent(inout)   :: En(0:Ntau-1,Nn)
real(8)    :: x(Nn),L,fvr(Nn,Nn)
integer(4) :: n,m,i
complex(8) :: fvptilde(Nn),temp(Nn),sum0(Nn)
complex(8) :: vctmp(Nn)
L=4.0d0
x=r
x(Nn/2+1)=1.0d0
do i=0,Ntau-1
    call fvrinterp(fh_fsl,tau(i),r,Nn,fvr)
    do n=1,Nn
    do m=1,Nn
        vctmp(m) = fvr(n,m)
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF, vctmp, fvptilde)
    sum0(n)=fvptilde(1)/dble(Nn)*(2.0d0*L)*r(n) !r*int_R fdv
    enddo
    call sll_s_fft_exec_c2c_1d(PlnF,sum0, fvptilde)
    do n=2,Nn
    fvptilde(n)=fvptilde(n)/sll_p_i1/lx(n)/dble(Nn)
    enddo
    fvptilde(1)=cmplx(0.0d0,0.0d0,kind=f64)
    call sll_s_fft_exec_c2c_1d(PlnB, fvptilde,temp)
    do n=1,Nn
    En(i,n)=dreal(temp(n)-temp(Nn/2+1))/x(n) !g(tau,r)
    enddo
enddo
!call fftw_destroy_plan(PlnF)
!call fftw_destroy_plan(PlnB)
end subroutine poissonsolver2

!----- gglobal(tau,xi1,xi2,E,lx)---
subroutine ge0(Nn,Ntau,tau,w0,En,Ent,Enr,gn,gnt,gnr)
use sll_m_boundary_condition_descriptors
use sll_m_cubic_splines
use sll_m_constants
implicit none
type(sll_t_cubic_spline_1D), pointer :: spl_fsl0, spl_fsl1,spl_fsl2
integer(4),intent(in)    :: Nn,Ntau
real(8),intent(in)       :: w0(Nn+1),tau(0:Ntau-1),En(0:Ntau-1,Nn),Enr(0:Ntau-1,Nn),Ent(0:Ntau-1,Nn)
real(8), intent(inout)   :: gn(0:Ntau-1,Nn+1,Nn+1),gnt(0:Ntau-1,Nn+1,Nn+1),gnr(0:Ntau-1,Nn+1,Nn+1)
real(8) :: x,L
integer(4) :: n,i,j
real(8) :: E0(Nn+1),E1(Nn+1),E2(Nn+1)
L=4.0d0
spl_fsl1 => sll_f_new_cubic_spline_1D(Nn+1, &
-L, L, &
SLL_P_PERIODIC)
spl_fsl0=>sll_f_new_cubic_spline_1D(Nn+1, &
-L, L, &
SLL_P_PERIODIC)
spl_fsl2=>sll_f_new_cubic_spline_1D(Nn+1, &
-L, L, &
SLL_P_PERIODIC)

do n=0,Ntau-1
    do j=1,Nn
    E0(j)=En(n,j)
    E1(j)=Enr(n,j)
    E2(j)=Ent(n,j)
    enddo
    E0(Nn+1)=0.0d0
    E1(Nn+1)=0.0d0
    E2(Nn+1)=0.0d0
    call sll_s_compute_cubic_spline_1D(E0,spl_fsl0)
    call sll_s_compute_cubic_spline_1D(E1,spl_fsl1)
    call sll_s_compute_cubic_spline_1D(E2,spl_fsl2)
    do j=1,Nn+1
    do i=1,Nn+1
        x=dcos(tau(n))*w0(i)+dsin(tau(n))*w0(j)
        if (dabs(x)<L) then
        gn(n,i,j)=sll_f_interpolate_from_interpolant_value(x,spl_fsl0)
        gnr(n,i,j)=sll_f_interpolate_from_interpolant_value(x,spl_fsl1)
        gnt(n,i,j)=sll_f_interpolate_from_interpolant_value(x,spl_fsl2)
        else
        gn(n,i,j)=0.0d0
        gnr(n,i,j)=0.0d0
        gnt(n,i,j)=0.0d0
        endif
    enddo
    enddo
enddo
call sll_o_delete(spl_fsl1)
call sll_o_delete(spl_fsl0)
call sll_o_delete(spl_fsl2)
end subroutine ge0

subroutine ge1(Nn,Ntau,tau,w1,w2,En,Ent,Enr,gn,gnt,gnr)
use sll_m_boundary_condition_descriptors
use sll_m_cubic_splines
use sll_m_constants
implicit none
type(sll_t_cubic_spline_1D), pointer :: spl_fsl0, spl_fsl1,spl_fsl2
integer(4),intent(in)    :: Nn,Ntau
real(8),intent(in)       :: tau(0:Ntau-1),En(0:Ntau-1,Nn),Enr(0:Ntau-1,Nn),Ent(0:Ntau-1,Nn)
real(8),intent(in)       :: w1(0:Ntau-1,Nn+1,Nn+1),w2(0:Ntau-1,Nn+1,Nn+1)
real(8), intent(inout)   :: gn(0:Ntau-1,Nn+1,Nn+1),gnt(0:Ntau-1,Nn+1,Nn+1),gnr(0:Ntau-1,Nn+1,Nn+1)
real(8) :: x,L
integer(4) :: n,i,j
real(8) :: E0(Nn+1),E1(Nn+1),E2(Nn+1)
L=4.0d0
spl_fsl1 => sll_f_new_cubic_spline_1D(Nn+1, &
-L, L, &
SLL_P_PERIODIC)
spl_fsl0=>sll_f_new_cubic_spline_1D(Nn+1, &
-L, L, &
SLL_P_PERIODIC)
spl_fsl2=>sll_f_new_cubic_spline_1D(Nn+1, &
-L, L, &
SLL_P_PERIODIC)

do n=0,Ntau-1
    do j=1,Nn
        E0(j)=En(n,j)
        E1(j)=Enr(n,j)
        E2(j)=Ent(n,j)
    enddo
    E0(Nn+1)=0.0d0
    E1(Nn+1)=0.0d0
    E2(Nn+1)=0.0d0
    call sll_s_compute_cubic_spline_1D(E0,spl_fsl0)
    call sll_s_compute_cubic_spline_1D(E1,spl_fsl1)
    call sll_s_compute_cubic_spline_1D(E2,spl_fsl2)
    do j=1,Nn+1
    do i=1,Nn+1
        x=dcos(tau(n))*w1(n,i,j)+dsin(tau(n))*w2(n,i,j)
        if (dabs(x)<L) then
            gn(n,i,j)=sll_f_interpolate_from_interpolant_value(x,spl_fsl0)
            gnr(n,i,j)=sll_f_interpolate_from_interpolant_value(x,spl_fsl1)
            gnt(n,i,j)=sll_f_interpolate_from_interpolant_value(x,spl_fsl2)
        else
            gn(n,i,j)=0.0d0
            gnr(n,i,j)=0.0d0
            gnt(n,i,j)=0.0d0
        endif
    enddo
    enddo
enddo
call sll_o_delete(spl_fsl1)
call sll_o_delete(spl_fsl0)
call sll_o_delete(spl_fsl2)
end subroutine ge1


subroutine ge2(Nn,Ntau,tau,w1,w2,En,gn)
use sll_m_boundary_condition_descriptors
use sll_m_cubic_splines
use sll_m_constants
implicit none
type(sll_t_cubic_spline_1D), pointer :: spl_fsl0
integer(4),intent(in)    :: Nn,Ntau
real(8),intent(in)       :: tau(0:Ntau-1),En(0:Ntau-1,Nn)
real(8),intent(in)       :: w1(0:Ntau-1,Nn+1,Nn+1),w2(0:Ntau-1,Nn+1,Nn+1)
real(8), intent(inout)   :: gn(0:Ntau-1,Nn+1,Nn+1)
real(8) :: x,L
integer(4) :: n,i,j
real(8) :: E0(Nn+1)
L=4.0d0
spl_fsl0=>sll_f_new_cubic_spline_1D(Nn+1, &
-L, L, &
SLL_P_PERIODIC)
do n=0,Ntau-1
    do j=1,Nn
        E0(j)=En(n,j)
    enddo
    E0(Nn+1)=0.0d0
    call sll_s_compute_cubic_spline_1D(E0,spl_fsl0)
    do j=1,Nn+1
    do i=1,Nn+1
        x=dcos(tau(n))*w1(n,i,j)+dsin(tau(n))*w2(n,i,j)
        if (dabs(x)<L) then
        gn(n,i,j)=sll_f_interpolate_from_interpolant_value(x,spl_fsl0)
        else
        gn(n,i,j)=0.0d0
        endif
    enddo
    enddo
enddo
call sll_o_delete(spl_fsl0)
end subroutine ge2

!double precision function seconds( )
!
!INTEGER ICOUNT, ICOUNT_RATE, ICOUNT_MAX
!
!call system_clock(icount,icount_rate,icount_max)
!seconds = icount/dble(icount_rate)
!
!
!RETURN
!
!end

