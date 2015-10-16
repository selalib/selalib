! conservation de la masse
!plot 'diag.dat' u 1:3 w lp title 'BSL', 'diag.dat' u 1:7 w lp title 'BSL NC', 'diag.dat' u 1:11 w lp title 'FSL', 'diag.dat' u 1:15 w lp title 'FSL NC'
! convergence en espace
!plot 'Conv_collela_rot_f3.dat' u 1:2 w lp title 'BSL', 'Conv_collela_rot_f3.dat' u 1:3 w lp title 'BSL NC', 'Conv_collela_rot_f3.dat' u 1:4 w lp title 'FSL', 'Conv_collela_rot_f3.dat' u 1:5 w lp title 'FSL NC'

program test_deposit_cubic_splines
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use sll_m_cubic_splines
use sll_m_constants
use sll_m_boundary_condition_descriptors
use sll_m_fft
use iso_c_binding
implicit none

#ifndef __INTEL_COMPILER

include "fftw3.f03"

type(C_PTR) :: PlnFwd,PlnBwd
type(sll_cubic_spline_2D), pointer :: spl_fsl
sll_int32  :: N,Neta1,Neta2,step,nb_step
sll_int32  :: i,j,bc1_type,bc2_type,err
sll_real64 :: eta1,delta_eta1,eta1_min,eta1_max,eta2,delta_eta2,eta2_min,eta2_max
sll_real64 :: x1_min,x2_min,x1_max,x2_max
sll_real64 :: T
sll_real64 :: eps
sll_real64 :: eta1c,eta2c
sll_real64,dimension(:,:), pointer :: f0
sll_real64,dimension(:,:), pointer :: fh_fsl
sll_real64,dimension(:,:), pointer :: x1_array,x2_array,eta1feet,eta2feet
sll_int32                         :: nc_eta1
sll_real64, dimension(:), pointer :: d_dx1
type(sll_fft_plan),       pointer :: fwx1
type(sll_fft_plan),       pointer :: bwx1
sll_comp64, dimension(:), pointer :: fk1
sll_int32     :: error

integer(4),parameter ::  ntau =256
complex(8) :: temp1(0:ntau-1),temp2(0:ntau-1),Ftilde2(0:ntau-1), AF1(0:ntau-1), AF2(0:ntau-1)
real(8)    :: k,h,r0,v0
complex(8) :: w10, w20,f,sumup1,sumup2,ff2
real(8)    :: tau(0:ntau-1), ltau(0:ntau-1)
complex(8) :: bw10(0:ntau-1),bw20(0:ntau-1),F1(0:ntau-1),F2(0:ntau-1),Ftilde1(0:ntau-1)
complex(8) :: Term2(0:ntau-1),Term1(0:ntau-1)
integer(4) :: m

N=128
Neta1 = N
Neta2 = N

! mesh type : cartesian
! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
! BC        : periodic-periodic
eta1_min = -4._f64
eta1_max = 4._f64
eta2_min = -4._f64
eta2_max = 4._f64

nc_eta1    = ntau


PlnFwd = fftw_plan_dft_1d(ntau,Ftilde2,Ftilde1,FFTW_FORWARD, FFTW_ESTIMATE+FFTW_UNALIGNED)
PlnBwd = fftw_plan_dft_1d(ntau,Ftilde1,Ftilde2,FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_UNALIGNED)
!SLL_CLEAR_ALLOCATE(d_dx1(1:nc_eta1+1), error)
!SLL_CLEAR_ALLOCATE(fk1(1:nc_eta1/2+1), error)
!fwx1 => fft_new_plan(nc_eta1, d_dx1, fk1)
!bwx1 => fft_new_plan(nc_eta1, fk1, d_dx1)
!
!SLL_CLEAR_ALLOCATE(kx1(1:nc_eta1/2+1), error)
!kx10 = 2._f64*sll_pi/(eta1_max-eta1_min)
!kx1(1) = 1.0_f64
!do i=2,nc_eta1/2+1
!   kx1(i) = (i-1)*kx10
!end do

! ---- * Parameters * ----

! --- Space and time parameters --

! Final time
T =1._f64

! -- mesh type --
! 1 : cartesian
!mesh_case = 1

! -- distribution function --
! 4 : centered sll_m_gaussian in eta1 and eta2
!test_case = 4

! -- advecton field --
! 1 : translation of vector (a1,a2)
! 2 : rotation
! 3 : non homogeneous rotation
! 4 : divergence free complex symmetric field (polar mesh only)
!field_case = 1

! -- visualization parameters --
!visu_step = 1
  
! ---- * Construction of the mesh * ----
  
bc1_type = SLL_PERIODIC
bc2_type = SLL_PERIODIC
  
eta1c = 0.5_f64*(eta1_max+eta1_min)
eta2c = 0.5_f64*(eta2_max+eta2_min)

! ---- * Time and space steps * ----

! space steps
delta_eta1 = (eta1_max-eta1_min)/real(Neta1,f64)
delta_eta2 = (eta2_max-eta2_min)/real(Neta2,f64)

! time step and number of steps
k = 0.05d0
nb_step = floor(T/k)
eps=0.1d0
h=2.0d0*sll_pi/ntau

! ---- * Messages * ----
  
print *,'# N=',N
print *,'# T=',T
print *,'# eps=',eps
!print *,'# mesh_case=',mesh_case
!print *,'# test_case=',test_case
!print *,'# field_case=',field_case
    
! ---- * Allocation and creation of the splines * ----
	
! allocations of the arrays
SLL_ALLOCATE(f0(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(fh_fsl(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(x1_array(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(x2_array(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(eta1feet(Neta1+1,Neta2+1), err)
SLL_ALLOCATE(eta2feet(Neta1+1,Neta2+1), err)
!SLL_ALLOCATE(diag(10,0:nb_step), err)

spl_fsl => new_cubic_spline_2D(Neta1+1, Neta2+1, &
  eta1_min, eta1_max, &
  eta2_min, eta2_max, &
  bc1_type, bc2_type)
  
! ---- * Initializations * ----
  
! Analytic distribution function and data for the mesh
open(unit=900,file='f0.dat')
do i=1,Neta1+1
  eta1 = eta1_min + (i-1)*delta_eta1
  do j=1,Neta2+1
    eta2 = eta2_min + (j-1)*delta_eta2
    
    x1_min = eta1_min
    x2_min = eta2_min
    x1_max = eta1_max
    x2_max = eta2_max
    
    x1_array(i,j) = eta1
    x2_array(i,j) = eta2
    
    f0(i,j) = exp(-2_f64*(eta1-eta1c)**2)*exp(-2_f64*(eta2-eta2c)**2)

    write(900,*) eta1, eta2, f0(i,j)
  enddo
  write(900,*)
enddo
close(900)

fh_fsl   = f0
fh_fsl(:,Neta2+1)   = fh_fsl(:,1)
call compute_cubic_spline_2D(fh_fsl,spl_fsl)


do m=0,ntau-1
    tau(m)=dble(m)*h
enddo
m=ntau/2
ltau=(/ (n, n=0,m-1), (n, n=-m,-1 )/)

do i=1,Neta1+1
do j=1,Neta2+1
    r0=x1_array(i,j)
    v0=x2_array(i,j)
    w10=(r0-sll_i1*v0)/2.0d0
    w20=-0.5d0*(r0+sll_i1*v0)
    do m=0,ntau-1
        f=(cdexp(sll_i1*tau(m))*w10-cdexp(-sll_i1*tau(m))*w20)*(dcos(2.0d0*tau(m)))**2.0d0
        F1(m)=-sll_i1*0.5d0*f*cdexp(-sll_i1*tau(m))
        F2(m)=-sll_i1*0.5d0*f*cdexp(sll_i1*tau(m))
    enddo
    call fft_apply_plan(fwx1, F1, Ftilde1)
    call fft_apply_plan(fwx1, F2, Ftilde2)
!    call fftw_execute_dft(PlnFwd, F1, Ftilde1)
!    call fftw_execute_dft(PlnFwd, F2, Ftilde2)
    Ftilde1=Ftilde1/dble(ntau)
    Ftilde2=Ftilde2/dble(ntau)
    do m=0,ntau-1
        sumup1=tau(m)*Ftilde1(0)
        sumup2=tau(m)*Ftilde2(0)
        do n=1,ntau-1
            sumup1=sumup1-sll_i1*Ftilde1(n)/ltau(n)*(cdexp(sll_i1*ltau(n)*tau(m))-1.0d0)
            sumup2=sumup2-sll_i1*Ftilde2(n)/ltau(n)*(cdexp(sll_i1*ltau(n)*tau(m))-1.0d0)
        enddo
        bw10(m)=w10-eps*tau(m)*Ftilde1(0)+eps*sumup1
        bw20(m)=w20-eps*tau(m)*Ftilde2(0)+eps*sumup2
    enddo
        !!! For 2nd order initial correction
    f=0.0d0
    ff2=0.0d0
    do n=1,ntau-1
        f=f+Ftilde1(n)/sll_i1/ltau(n)  !AF1_0
        ff2=ff2+Ftilde2(n)/sll_i1/ltau(n) !AF2_0
    enddo
    do m=0,ntau-1
        temp1(m)=tau(m)*Ftilde1(0)
        temp2(m)=tau(m)*Ftilde2(0)
        do n=1,ntau-1
            temp1(m)=(cdexp(sll_i1*tau(m)*ltau(n))-1.0d0)*(Ftilde1(n)/sll_i1/ltau(n))+temp1(m)
            temp2(m)=(cdexp(sll_i1*tau(m)*ltau(n))-1.0d0)*(Ftilde2(n)/sll_i1/ltau(n))+temp2(m)
        enddo
    enddo
    AF1=temp1-tau*Ftilde1(0)+f
    AF2=temp2-tau*Ftilde2(0)+ff2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do m=0,ntau-1
        sumup1=(cdexp(sll_i1*tau(m))*f-cdexp(-sll_i1*tau(m))*ff2)*(dcos(2*tau(m)))**2.0d0
        temp1(m)=-0.5d0*sll_i1*cdexp(-sll_i1*tau(m))*sumup1 !dwFh10_1
        temp2(m)=-0.5d0*sll_i1*cdexp(sll_i1*tau(m))*sumup1  !dwFh10_2
    enddo
    call fft_apply_plan(fwx1,temp1, Term1)
    call fft_apply_plan(fwx1,temp2, Term2)
!    call fftw_execute_dft(PlnFwd,temp1, Term1)
!    call fftw_execute_dft(PlnFwd,temp2, Term2)
    Term1=Term1/dble(ntau)
    Term2=Term2/dble(ntau)
    do m=0,ntau-1
        temp1(m)=tau(m)*Term1(0)
        temp2(m)=tau(m)*Term2(0)
        do n=1,ntau-1
            temp1(m)=(cdexp(sll_i1*tau(m)*ltau(n))-1.0d0)*(Term1(n)/sll_i1/ltau(n))+temp1(m)
            temp2(m)=(cdexp(sll_i1*tau(m)*ltau(n))-1.0d0)*(Term2(n)/sll_i1/ltau(n))+temp2(m)
        enddo
    enddo
    bw10=bw10-eps**2.0d0*(temp1-tau*Term1(0))
    bw20=bw20-eps**2.0d0*(temp2-tau*Term2(0))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do m=0,ntau-1
        f=(cdexp(sll_i1*tau(m))*AF1(m)-cdexp(-sll_i1*tau(m))*AF2(m))*(dcos(2*tau(m)))**2
        temp1(m)=-0.5d0*sll_i1*cdexp(-sll_i1*tau(m))*f !dwFAF1
        temp2(m)=-0.5d0*sll_i1*cdexp(sll_i1*tau(m))*f  !dwFAF2
        f=(cdexp(sll_i1*tau(m))*Ftilde1(0)-cdexp(-sll_i1*tau(m))*Ftilde2(0))*(dcos(2*tau(m)))**2
        AF1(m)=-0.5d0*sll_i1*cdexp(-sll_i1*tau(m))*f   !dwFPiF1
        AF2(m)=-0.5d0*sll_i1*cdexp(sll_i1*tau(m))*f    !dwFPiF2
    enddo
    call fft_apply_plan(fwx1, AF1, F1)
    call fft_apply_plan(fwx1, AF2, F2)
!    call fftw_execute_dft(PlnFwd, AF1, F1)
!    call fftw_execute_dft(PlnFwd, AF2, F2)
    F1=F1/dble(ntau)
    F2=F2/dble(ntau)
    f=0.0d0
    ff2=0.0d0
    do n=1,ntau-1
        f=f+F1(n)/sll_i1/ltau(n)  !AdwFPiF1_0
        ff2=ff2+F2(n)/sll_i1/ltau(n) !AdwFPiF2_0
    enddo
    do m=0,ntau-1
        Ftilde1(m)=tau(m)*F1(0)
        Ftilde2(m)=tau(m)*F2(0)
        do n=1,ntau-1
            Ftilde1(m)=(cdexp(sll_i1*tau(m)*ltau(n))-1.0d0)*(F1(n)/sll_i1/ltau(n))+Ftilde1(m)
            Ftilde2(m)=(cdexp(sll_i1*tau(m)*ltau(n))-1.0d0)*(F2(n)/sll_i1/ltau(n))+Ftilde2(m)
        enddo
    enddo
    AF1=Ftilde1-tau*F1(0)+f    !AdwFPiF1
    AF2=Ftilde2-tau*F2(0)+ff2  !AdwFPiF2
    F1=temp1-AF1 !z1
    F2=temp2-AF2
    call fft_apply_plan(fwx1, F1, Ftilde1)
    call fft_apply_plan(fwx1, F2, Ftilde2)
!    call fftw_execute_dft(PlnFwd, F1, Ftilde1)
!    call fftw_execute_dft(PlnFwd, F2, Ftilde2)
    Ftilde1=Ftilde1/dble(ntau)
    Ftilde2=Ftilde2/dble(ntau)
    do m=0,ntau-1
        temp1(m)=tau(m)*Ftilde1(0)
        temp2(m)=tau(m)*Ftilde2(0)
        do n=1,ntau-1
            temp1(m)=(cdexp(sll_i1*tau(m)*ltau(n))-1.0d0)*(Ftilde1(n)/sll_i1/ltau(n))+temp1(m)
            temp2(m)=(cdexp(sll_i1*tau(m)*ltau(n))-1.0d0)*(Ftilde2(n)/sll_i1/ltau(n))+temp2(m)
        enddo
    enddo
    AF1=temp1-tau*Ftilde1(0) !h2_1
    AF2=temp2-tau*Ftilde2(0) !h2_2
    bw10=bw10+eps**2.0d0*AF1
    bw20=bw20+eps**2.0d0*AF2
    do step=1,nb_step ! ---- * Evolution in time * ----
        do m=0,ntau-1
            f=(cdexp(sll_i1*tau(m))*bw10(m)-cdexp(-sll_i1*tau(m))*bw20(m))*(dcos(2.0d0*tau(m)))**2.0d0
            F1(m)=-f*cdexp(-sll_i1*tau(m))*0.5d0*sll_i1
            F2(m)=-0.5d0*sll_i1*f*cdexp(sll_i1*tau(m))
        enddo
        temp1=bw10+k/2.0d0*F1
        temp2=bw20+k/2.0d0*F2
        call fft_apply_plan(fwx1, temp1, Ftilde1)
        call fft_apply_plan(fwx1, temp2, Ftilde2)
!        call fftw_execute_dft(PlnFwd, temp1, Ftilde1)
!        call fftw_execute_dft(PlnFwd, temp2, Ftilde2)
        do m=0,ntau-1
            Ftilde1(m)=Ftilde1(m)/(1.0d0+sll_i1*k/2.0d0*ltau(m)/eps)
            Ftilde2(m)=Ftilde2(m)/(1.0d0+sll_i1*k/2.0d0*ltau(m)/eps)
        enddo
        call fft_apply_plan(bwx1, Ftilde1,AF1)
        call fft_apply_plan(bwx1, Ftilde2,AF2)
!        call fftw_execute_dft(PlnBwd, Ftilde1,AF1)
!        call fftw_execute_dft(PlnBwd, Ftilde2,AF2)
        AF1=AF1/dble(ntau) !w1_h
        AF2=AF2/dble(ntau) !w2_h
        do m=0,ntau-1
            f=(cdexp(sll_i1*tau(m))*AF1(m)-cdexp(-sll_i1*tau(m))*AF2(m))*(dcos(2.0d0*tau(m)))**2.0d0
            F1(m)=-f*cdexp(-sll_i1*tau(m))*0.5d0*sll_i1
            F2(m)=-0.5d0*sll_i1*f*cdexp(sll_i1*tau(m))
        enddo
        call fft_apply_plan(fwx1, F1,  Ftilde1)
        call fft_apply_plan(fwx1, F2,  Ftilde2)
        call fft_apply_plan(fwx1, bw10, F1)
        call fft_apply_plan(fwx1, bw20, F2)
!        call fftw_execute_dft(PlnFwd, F1,  Ftilde1)
!        call fftw_execute_dft(PlnFwd, F2,  Ftilde2)
!        call fftw_execute_dft(PlnFwd, bw10, F1)
!        call fftw_execute_dft(PlnFwd, bw20, F2)
        do m=0,ntau-1
            temp1(m)=(F1(m)*(1-sll_i1*k/eps/2.0d0*ltau(m))+k*Ftilde1(m))/(1.0d0+sll_i1*k/2.0d0*ltau(m)/eps)
            temp2(m)=(F2(m)*(1-sll_i1*k/eps/2.0d0*ltau(m))+k*Ftilde2(m))/(1.0d0+sll_i1*k/2.0d0*ltau(m)/eps)
        enddo
        call fft_apply_plan(fwx1, temp1,bw10)
        call fft_apply_plan(fwx1, temp2,bw20)
!        call fftw_execute_dft(PlnBwd, temp1,bw10)
!        call fftw_execute_dft(PlnBwd, temp2,bw20)
        bw10=bw10/dble(ntau)
        bw20=bw20/dble(ntau)
    enddo
    call fft_apply_plan(fwx1,bw10,temp1)
    call fft_apply_plan(fwx1,bw20,temp2)
!    call fftw_execute_dft(PlnFwd,bw10,temp1)
!    call fftw_execute_dft(PlnFwd,bw20,temp2)
    temp1=temp1/dble(ntau)
    temp2=temp2/dble(ntau)
    sumup1=0.0d0
    sumup2=0.0d0
    do n=0,ntau-1
        sumup1=sumup1+temp1(n)*cdexp(sll_i1*ltau(n)*T/eps)
        sumup2=sumup2+temp2(n)*cdexp(sll_i1*ltau(n)*T/eps)
    enddo
    eta1=dreal(cdexp(sll_i1*T/eps)*sumup1-cdexp(-sll_i1*T/eps)*sumup2)
    eta2=dreal(sll_i1*cdexp(sll_i1*T/eps)*sumup1+sll_i1*cdexp(-sll_i1*T/eps)*sumup2)
    call apply_bc()
    eta1feet(i,j)=eta1
    eta2feet(i,j)=eta2
enddo
enddo

call fft_delete_plan(fwx1)
call fft_delete_plan(bwx1)
!call fftw_destroy_plan(PlnFwd)
!call fftw_destroy_plan(PlnBwd)
  ! --- Deposition FSL ---
    
call deposit_value_2D(eta1feet,eta2feet,spl_fsl,fh_fsl)

! File name

!mesh_name = "crt"
!
!SELECT CASE (field_case)
!  CASE (1)
!    field_name = "trs"
!  CASE (2)
!    field_name = "rot"
!  CASE (3)
!    field_name = "rnh"
!END SELECT

!i1 = int(T)/100
!i2 =(int(T)-100*i1)/10
!i3 = int(T)-100*i1-10*i2
!time_name = char(i1+48)//char(i2+48)//char(i3+48)
!
!conv_name = 'Conv_'//mesh_name//'_'//field_name//'_'//time_name//'.dat'
!mass_name = 'Mass_'//mesh_name//'_'//field_name//'_'//time_name//'.dat'
!  
!val_bsl    = maxval(abs(f-fh_bsl))
!val_fsl    = maxval(abs(f-fh_fsl))
!val_spe    = maxval(abs(f-fh_spe))
!val        = maxval(abs(fh_fsl-fh_bsl))

!write(*,*) N,'fsl:',val_fsl
open(unit=850,file='fh.dat')  
do i=1,Neta1+1
  do j=1,Neta2+1
    write(850,*) x1_array(i,j),x2_array(i,j),fh_fsl(i,j)
  enddo
  write(850,*) ' '
enddo
close(850)

#endif /* __INTEL_COMPILER */

end program

subroutine apply_bc()

  ! --- Corrections on the BC ---
  if (bc1_type.eq.SLL_HERMITE) eta1 = min(max(eta1,eta1_min),eta1_max)
  if (bc2_type.eq.SLL_HERMITE) eta2 = min(max(eta2,eta2_min),eta2_max)

  if (bc1_type==SLL_PERIODIC) then
    do while (eta1>eta1_max)
      eta1 = eta1-(eta1_max-eta1_min)
    enddo
    do while (eta1<eta1_min)
      eta1 = eta1+(eta1_max-eta1_min)
    enddo
  endif

  if (bc2_type==SLL_PERIODIC) then
    do while (eta2>eta2_max)
      eta2 = eta2-(eta2_max-eta2_min)
    enddo
    do while (eta2<eta2_min)
      eta2 = eta2+(eta2_max-eta2_min)
    enddo
  endif

end subroutine apply_bc



