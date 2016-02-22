! Mouton example
! N_tau optimised
! taut optimised; 1d nufft corrected
! in rotating framework
program test_deposit_cubic_splines
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

use deposit_cubic_splines
use sll_m_cubic_splines
use sll_m_constants
use sll_m_boundary_condition_descriptors
use sll_m_fft

implicit none

sll_int32, parameter :: Nn   = 128
sll_int32, parameter :: ntau = 32

type(sll_t_fft) :: PlnFwd
type(sll_t_fft) :: PlnBwd
type(sll_t_fft) :: PlnF
type(sll_t_fft) :: PlnB

type(sll_t_cubic_spline_2D), pointer :: spl_fsl

sll_int32  :: step,nb_step
sll_int32  :: i,j,bc1_type,bc2_type,err
sll_real64 :: eta1,delta_eta1,eta1_min,eta1_max,eta2
sll_real64 :: xi1_0,xi2_0,x
sll_real64 :: T,eps

sll_real64, dimension(:,:), allocatable :: eta1feet
sll_real64, dimension(:,:), allocatable :: eta2feet

sll_real64 :: x1_array(Nn+1)
sll_real64 :: r_array(Nn)
sll_real64 :: t1
sll_real64 :: t2
sll_real64 :: fh_fsl(Nn+1,Nn+1)
sll_real64 :: f0(Nn,Nn)
sll_comp64 :: temp1(0:ntau-1)
sll_comp64 :: temp2(0:ntau-1)
sll_comp64 :: AF1(0:ntau-1)
sll_comp64 :: AF2(0:ntau-1)
sll_comp64 :: w1c(0:ntau-1)
sll_comp64 :: w2c(0:ntau-1)
sll_real64 :: k,h
sll_real64 :: En(0:Ntau-1,1:Nn)
sll_real64 :: Enr(0:Ntau-1,1:Nn)
sll_real64 :: Ent(0:Ntau-1,1:Nn)
sll_real64 :: gn(0:Ntau-1,Nn+1,Nn+1)
sll_real64 :: gnr(0:Ntau-1,Nn+1,Nn+1)
sll_real64 :: gnt(0:Ntau-1,Nn+1,Nn+1)
sll_comp64 :: sumup1
sll_comp64 :: sumup2
sll_comp64 :: dtgn
sll_real64 :: tau(0:ntau-1)
sll_real64 :: ltau(0:ntau-1)
sll_real64 :: lx(1:Nn)
sll_real64 :: fvr(1:Nn,1:Nn)
sll_real64 :: taut(0:ntau-1)
sll_real64 :: w1_0(0:ntau-1,Nn+1,Nn+1)
sll_real64 :: w2_0(0:ntau-1,Nn+1,Nn+1)
sll_comp64 :: F1(0:ntau-1)
sll_comp64 :: F2(0:ntau-1)
sll_comp64 :: Ftilde1(0:ntau-1,Nn+1,Nn+1)
sll_comp64 :: Ftilde2(0:ntau-1,Nn+1,Nn+1)
sll_comp64 :: Term1(0:ntau-1)
sll_comp64 :: Term2(0:ntau-1)
sll_comp64 :: temp1_F(0:ntau-1)
sll_comp64 :: temp2_F(0:ntau-1)
sll_comp64 :: dtF1(0:ntau-1)
sll_comp64 :: dtF2(0:ntau-1)
sll_int32  :: n,m

sll_comp64,  allocatable :: fgen(:)

! ---- * Parameters * ----

! mesh type : cartesian
! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
! BC        : periodic-periodic
eta1_min = -4.0_f64
eta1_max =  4.0_f64

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
k       = 0.05d0  !T*delta_eta1
nb_step = floor(T/k)
eps     = 1.00d0
h       = 2.0d0 * sll_p_pi/ntau

! ---- * Messages * ----

print *,'# N=',Nn
print *,'# k=',k
print *,'# T=',T
print *,'# eps=',eps

! ---- * Allocation and creation of the splines * ----
!L=4.0d0
allocate (fgen(Nn))
call sll_s_fft_init_c2c_1d(PlnFwd ,ntau, AF2,   AF1, sll_p_fft_forward)
call sll_s_fft_init_c2c_1d(PlnBwd ,ntau, AF1,   AF2, sll_p_fft_forward)
call sll_s_fft_init_c2c_1d(PlnF   ,Nn,   Fgen, Fgen, sll_p_fft_forward)
call sll_s_fft_init_c2c_1d(PlnB   ,Nn,   Fgen, Fgen, sll_p_fft_forward)
deallocate (fgen)

! allocations of the arrays
SLL_ALLOCATE(eta1feet(Nn+1,Nn+1), err)
SLL_ALLOCATE(eta2feet(Nn+1,Nn+1), err)
spl_fsl => sll_f_new_cubic_spline_2D( Nn+1,     &
                                      Nn+1,     &
                                      eta1_min, &
                                      eta1_max, &
                                      eta1_min, &
                                      eta1_max, &
                                      bc1_type, &
                                      bc2_type)

! ---- * Initializations * ----

! Analytic distribution function and data for the mesh
do i=1,Nn+1
  eta1 = eta1_min + (i-1)*delta_eta1
  x1_array(i) = eta1
  do j=1,Nn+1
    eta2 = eta1_min + (j-1)*delta_eta1
    fh_fsl(i,j) = dexp(-2.0d0*(eta1**2+eta2**2))
  enddo
enddo
do m=0,ntau-1
  tau(m)=dble(m)*h
enddo
m       = ntau/2
ltau    = [ (real(n,f64), n=0,m-1), (real(n,f64), n=-m,-1 ) ]
m       = Nn/2
lx      = [ (real(n,f64), n=0,m-1), (real(n,f64), n=-m,-1 ) ]*2.0d0*sll_p_pi/(eta1_max-eta1_min)
t       = 0.0d0
r_array = x1_array(1:Nn)

!t1= second();
!-------- * Evolution in time * ---------
do step=1,nb_step

  taut=tau+t
  f0=fh_fsl(1:Nn,1:Nn)
  call sll_s_compute_cubic_spline_2D(fh_fsl,spl_fsl)
  call poissonsolver(f0,r_array,taut,lx,Nn,Ntau,En,Enr,Ent,PlnF,PlnB)
  call ge0(Nn,Ntau,taut,x1_array,En,Ent,Enr,gn,gnt,gnr)

  do i=1,Nn+1
  do j=1,Nn+1
    xi1_0=x1_array(i)
    xi2_0=x1_array(j)
    !------------for 1st order correction------------------
    do m=0,ntau-1
      F1(m)=(       -cos(2.0d0*taut(m))**2
             *( 0.5d0*sin(2.0d0*taut(m))*xi1_0
               + sin(taut(m))**2*xi2_0)
               - sin(taut(m))*gn(m,i,j))/real(ntau,f64)

      F2(m)=(        cos(2.0d0*taut(m))**2
             *( cos(taut(m))**2*xi1_0
               + 0.5d0 * sin(2.0d0*taut(m))*xi2_0)
               + cos(taut(m))*gn(m,i,j))/real(ntau, f64)
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
      F1(m)=-dcos(2.0d0*taut(m))**2*(0.5d0*dsin(2.0d0*taut(m))*w1_0(m,i,j) &
            +dsin(taut(m))**2*w2_0(m,i,j))-dsin(taut(m))*gn(m,i,j)
      F2(m)=dcos(2.0d0*taut(m))**2*(dcos(taut(m))**2*w1_0(m,i,j) &
           +0.5d0*dsin(2.0d0*taut(m))*w2_0(m,i,j))+dcos(taut(m))*gn(m,i,j)
      dtgn=gnt(m,i,j)+gnr(m,i,j)*(dcos(taut(m))*Ftilde1(0,i,j) &
           +dsin(taut(m))*Ftilde2(0,i,j))
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

f0=fh_fsl(1:Nn,1:Nn)
call fvrinterp(f0,t,r_array,Nn,fvr)
call sll_s_fft_free(PlnFwd)
call sll_s_fft_free(PlnBwd)
call sll_s_fft_free(PlnF)
call sll_s_fft_free(PlnB)

open(unit=850,file='fh.dat')
do i=1,Nn
  do j=1,Nn
    write(850,*) i, j, sngl(fvr(i,j))
  enddo
  write(850,*)
enddo
close(850)

contains

!> Corrections on the BC 
subroutine apply_bc()
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
if (dabs(eta1)<1.0d-12) eta1=0.0d0
if (dabs(eta2)<1.0d-12) eta2=0.0d0
end subroutine apply_bc

end program test_deposit_cubic_splines

