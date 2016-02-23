! Mouton example
! N_tau optimised
! taut optimised; 1d nufft corrected
! in rotating framework
program test_deposit_cubic_splines
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

use deposit_cubic_splines
use poisson_solver
use sll_m_cubic_splines
use sll_m_constants
use sll_m_boundary_condition_descriptors
use sll_m_fft
use sll_m_gnuplot
use sll_m_xdmf

implicit none

sll_comp64, parameter :: sll_p_i0 = (0.0_f64, 0.0_f64)
sll_int32,  parameter :: n        = 128
sll_int32,  parameter :: ntau     = 32

type(sll_t_fft) :: fw_fft
type(sll_t_fft) :: bw_fft

type(sll_t_cubic_spline_2d), pointer :: spl_2d

sll_int32  :: step,nb_step
sll_int32  :: i,j,l
sll_int32  :: bc1_type,bc2_type,err
sll_real64 :: delta_eta1
sll_real64 :: eta1_min
sll_real64 :: eta1_max
sll_real64 :: xi1_0,xi2_0
sll_real64 :: T,eps
sll_real64 :: eta1,  eta2
sll_real64 :: x1(n+1), x2(n+1)

sll_real64, dimension(:,:), allocatable :: eta1feet
sll_real64, dimension(:,:), allocatable :: eta2feet

sll_comp64 :: F1(0:ntau-1)
sll_comp64 :: F2(0:ntau-1)
sll_comp64 :: tmp1(0:ntau-1)
sll_comp64 :: tmp2(0:ntau-1)
sll_comp64 :: tmp1_f(0:ntau-1)
sll_comp64 :: tmp2_f(0:ntau-1)
sll_comp64 :: w1c(0:ntau-1)
sll_comp64 :: w2c(0:ntau-1)

sll_real64 :: tau(0:ntau-1)
sll_real64 :: ltau(0:ntau-1)

sll_real64 :: fh_fsl(n+1,n+1)
sll_real64 :: f0(n,n)
sll_real64 :: k, h
sll_real64 :: En(0:ntau-1,1:n)
sll_real64 :: Enr(0:ntau-1,1:n)
sll_real64 :: Ent(0:ntau-1,1:n)
sll_real64 :: gn(0:ntau-1,n+1,n+1)
sll_real64 :: gnr(0:ntau-1,n+1,n+1)
sll_real64 :: gnt(0:ntau-1,n+1,n+1)
sll_comp64 :: sumup1
sll_comp64 :: sumup2
sll_comp64 :: dtgn
sll_real64 :: lx(1:n)
sll_real64 :: fvr(1:n,1:n)
sll_real64 :: taut(0:ntau-1)
sll_real64 :: w1_0(0:ntau-1,n+1,n+1)
sll_real64 :: w2_0(0:ntau-1,n+1,n+1)
sll_comp64 :: Ftilde1(0:ntau-1,n+1,n+1)
sll_comp64 :: Ftilde2(0:ntau-1,n+1,n+1)
sll_int32  :: m, ref_id
sll_real32 :: fdum, error
sll_real64 :: tstart, tend
sll_int32  :: ierr


type(poisson) :: solver

! ---- * Parameters * ----

! mesh type : cartesian
! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
! BC        : periodic-periodic
eta1_min = -4.0_f64
eta1_max =  4.0_f64

! --- Space and time parameters --

! Final time
T = 0.4_f64

! ---- * Construction of the mesh * ----
bc1_type = SLL_P_PERIODIC
bc2_type = SLL_P_PERIODIC

! ---- * Time and space steps * ----

! space steps
delta_eta1 = (eta1_max-eta1_min)/real(n,f64)

! time step and number of steps
k       = 0.05d0  !T*delta_eta1
nb_step = floor(T/k)
eps     = 1.00d0
h       = 2.0d0 * sll_p_pi/ real(ntau,f64)

! ---- * Messages * ----

print *,'# N=',n
print *,'# k=',k
print *,'# T=',T
print *,'# eps=',eps

call cpu_time(tstart)

call solver%init( eta1_min, eta1_max, n)

call sll_s_fft_init_c2c_1d(fw_fft ,ntau, tmp1, tmp1_f, sll_p_fft_forward)
call sll_s_fft_init_c2c_1d(bw_fft ,ntau, tmp1_f, tmp1, sll_p_fft_backward)

! allocations of the arrays
SLL_ALLOCATE(eta1feet(n+1,n+1), err)
SLL_ALLOCATE(eta2feet(n+1,n+1), err)
spl_2d => sll_f_new_cubic_spline_2d( n+1,     &
                                     n+1,     &
                                     eta1_min, &
                                     eta1_max, &
                                     eta1_min, &
                                     eta1_max, &
                                     bc1_type, &
                                     bc2_type)

! ---- * Initializations * ----

! Analytic distribution function and data for the mesh
do i=1,n+1
  x1(i) = eta1_min + (i-1)*delta_eta1
enddo
do j=1,n+1
  x2(j) = eta1_min + (j-1)*delta_eta1
enddo

do i=1,n+1
  do j=1,n+1
    fh_fsl(i,j) = exp(-2.0d0*(x1(i)**2+x2(j)**2))
  end do
end do
do m=0,ntau-1
  tau(m)=real(m,f64)*h
enddo

m    = ntau/2
ltau = [(real(l,f64),l=0,m-1),(real(l,f64),l=-m,-1)]
m    = n/2
lx   = [(real(l,f64),l=0,m-1),(real(l,f64),l=-m,-1)]*2.0d0*sll_p_pi/(eta1_max-eta1_min)
t    = 0.0d0

!t1= second();
!-------- * Evolution in time * ---------
do step=1,nb_step

  taut = tau + t
  f0   = fh_fsl(1:n,1:n)
  call sll_s_compute_cubic_spline_2d(fh_fsl,spl_2d)
  call solver%solve1(f0,taut,n,ntau,En,Enr,Ent)
  call ge0(n,ntau,taut,x1,En,Ent,Enr,gn,gnt,gnr)

  do i=1,n+1
    do j=1,n+1

      xi1_0 = x1(i)
      xi2_0 = x2(j)

      !------------for 1st order correction------------------
      do m=0,ntau-1
        F1(m) = rfct1(taut(m), ntau, xi1_0, xi2_0, gn(m,i,j))
        F2(m) = rfct2(taut(m), ntau, xi1_0, xi2_0, gn(m,i,j))
      enddo

      call sll_s_fft_exec_c2c_1d(fw_fft, F1, F1)
      call sll_s_fft_exec_c2c_1d(fw_fft, F2, F2)

      Ftilde1(:,i,j) = F1
      Ftilde2(:,i,j) = F2

      tmp1(1:ntau-1) = - sll_p_i1 * F1(1:ntau-1) / ltau(1:ntau-1)
      tmp2(1:ntau-1) = - sll_p_i1 * F2(1:ntau-1) / ltau(1:ntau-1)

      tmp1(0) = sll_p_i0
      tmp2(0) = sll_p_i0

      call sll_s_fft_exec_c2c_1d(bw_fft, tmp1, F1)
      call sll_s_fft_exec_c2c_1d(bw_fft, tmp2, F2)

      F1 = F1 - sum(tmp1)
      F2 = F2 - sum(tmp2)

      w1_0(:,i,j) = xi1_0 + eps*real(F1)
      w2_0(:,i,j) = xi2_0 + eps*real(F2)

    enddo
  enddo

  call ge1(n,ntau,taut,w1_0,w2_0,En,Ent,Enr,gn,gnt,gnr) ! gn,gnt,gnr 3d array output

  do i=1,n+1
    do j=1,n+1
      !------------for 2nd order correction------------------
      do m=0,ntau-1

        F1(m)= rfct1(taut(m), 1, w1_0(m,i,j), w2_0(m,i,j), gn(m,i,j))
        F2(m)= rfct2(taut(m), 1, w1_0(m,i,j), w2_0(m,i,j), gn(m,i,j))

        dtgn=  cmplx(gnt(m,i,j),0.0,f64) &
             + cmplx(gnr(m,i,j),0.0,f64) &
             *(cmplx(cos(taut(m)),0.0,f64)*Ftilde1(0,i,j) &
             + cmplx(sin(taut(m)),0.0,f64)*Ftilde2(0,i,j))

        tmp1(m)= fct1(taut(m), ntau, Ftilde1(0,i,j), Ftilde2(0,i,j), dtgn)
        tmp2(m)= fct2(taut(m), ntau, Ftilde1(0,i,j), Ftilde2(0,i,j), dtgn)

      enddo

      call sll_s_fft_exec_c2c_1d(fw_fft, tmp1, tmp1_f)
      call sll_s_fft_exec_c2c_1d(fw_fft, tmp2, tmp2_f)

      tmp1_f(1:ntau-1) = -sll_p_i1*tmp1_f(1:ntau-1)/ltau(1:ntau-1)
      tmp2_f(1:ntau-1) = -sll_p_i1*tmp2_f(1:ntau-1)/ltau(1:ntau-1)

      tmp1_f(0) = sll_p_i0
      tmp2_f(0) = sll_p_i0

      call sll_s_fft_exec_c2c_1d(bw_fft, tmp1_f, tmp1)
      call sll_s_fft_exec_c2c_1d(bw_fft, tmp2_f, tmp2)

      tmp1 = (tmp1-sum(tmp1_f))/cmplx(ntau,0.0,f64)
      tmp2 = (tmp2-sum(tmp2_f))/cmplx(ntau,0.0,f64)

      call sll_s_fft_exec_c2c_1d(fw_fft,tmp1,tmp1_f)
      call sll_s_fft_exec_c2c_1d(fw_fft,tmp2,tmp2_f)

      tmp1_f = (F1-(tmp1-tmp1_f(0))*eps)/cmplx(ntau,0.0,f64)
      tmp2_f = (F2-(tmp2-tmp2_f(0))*eps)/cmplx(ntau,0.0,f64)

      call sll_s_fft_exec_c2c_1d(fw_fft, tmp1_f, tmp1_f)
      call sll_s_fft_exec_c2c_1d(fw_fft, tmp2_f, tmp2_f)

      tmp1_f(1:ntau-1) = -sll_p_i1*tmp1_f(1:ntau-1)/ltau(1:ntau-1)
      tmp2_f(1:ntau-1) = -sll_p_i1*tmp2_f(1:ntau-1)/ltau(1:ntau-1)

      tmp1_f(0) = sll_p_i0
      tmp2_f(0) = sll_p_i0

      call sll_s_fft_exec_c2c_1d(bw_fft, tmp1_f, tmp1)
      call sll_s_fft_exec_c2c_1d(bw_fft, tmp2_f, tmp2)

      tmp1 = tmp1 - sum(tmp1_f)
      tmp2 = tmp2 - sum(tmp2_f)

      w1_0(:,i,j) = x1(i) + eps*real(tmp1)
      w2_0(:,i,j) = x2(j) + eps*real(tmp2)

    enddo
  enddo

  call ge2(n,ntau,taut,w1_0,w2_0,En,gn)

  do i=1,n+1
    do j=1,n+1

      do m=0,ntau-1
        F1(m)=rfct1(taut(m),1,w1_0(m,i,j),w2_0(m,i,j),gn(m,i,j))
        F2(m)=rfct2(taut(m),1,w1_0(m,i,j),w2_0(m,i,j),gn(m,i,j))
      enddo

      tmp1 = w1_0(:,i,j) + 0.5_f64*k*F1
      tmp2 = w2_0(:,i,j) + 0.5_f64*k*F2

      call sll_s_fft_exec_c2c_1d(fw_fft, tmp1, tmp1_f)
      call sll_s_fft_exec_c2c_1d(fw_fft, tmp2, tmp2_f)

      tmp1_f = tmp1_f/(1.0d0+sll_p_i1*k/2.0d0*ltau/eps)/cmplx(ntau,0.0,f64)
      tmp2_f = tmp2_f/(1.0d0+sll_p_i1*k/2.0d0*ltau/eps)/cmplx(ntau,0.0,f64)

      call sll_s_fft_exec_c2c_1d(bw_fft, tmp1_f, tmp1)
      call sll_s_fft_exec_c2c_1d(bw_fft, tmp2_f, tmp2)

      Ftilde1(:,i,j) = tmp1
      Ftilde2(:,i,j) = tmp2

      !----------Insert half step evaluation--------

      sumup1 = sum(tmp1_f*exp(0.5*sll_p_i1*ltau*k/eps))
      sumup2 = sum(tmp2_f*exp(0.5*sll_p_i1*ltau*k/eps))

      eta1=dreal(sumup1)
      eta2=dreal(sumup2)

      call apply_bc()

      eta1feet(i,j)=eta1
      eta2feet(i,j)=eta2

    enddo
  enddo

  call sll_s_deposit_value_2d(eta1feet,eta2feet,spl_2d,fh_fsl) !function value at the half time

  f0=fh_fsl(1:n,1:n)

  call solver%solve2(f0,taut,n,ntau,En)

  !------End evaluation and continue 2nd solver-------------

  call ge2(n,ntau,taut,dreal(Ftilde1),dreal(Ftilde2),En,gn)

  do i=1,n+1
    do j=1,n+1

      do m=0,ntau-1
        F1(m)=fct1(taut(m), 1, Ftilde1(m,i,j), Ftilde2(m,i,j), cmplx(gn(m,i,j),0.0,f64))
        F2(m)=fct2(taut(m), 1, Ftilde1(m,i,j), Ftilde2(m,i,j), cmplx(gn(m,i,j),0.0,f64))
      enddo

      call sll_s_fft_exec_c2c_1d(fw_fft, F1, tmp1_f)
      call sll_s_fft_exec_c2c_1d(fw_fft, F2, tmp2_f)

      w1c = cmplx(w1_0(:,i,j),0.0,f64)
      w2c = cmplx(w2_0(:,i,j),0.0,f64)

      call sll_s_fft_exec_c2c_1d(fw_fft, w1c, F1)
      call sll_s_fft_exec_c2c_1d(fw_fft, w2c, F2)

      tmp1=(F1*(1.0_f64-sll_p_i1*0.5_f64*k/eps*ltau)+k*tmp1_f) &
              /(1.0_f64+sll_p_i1*0.5_f64*k*ltau/eps)
      tmp2=(F2*(1.0_f64-sll_p_i1*0.5_f64*k/eps*ltau)+k*tmp2_f) &
              /(1.0_f64+sll_p_i1*0.5_f64*k*ltau/eps)

      sumup1 = sum(tmp1*exp(sll_p_i1*ltau*k/eps))/cmplx(ntau,0.0,f64)
      sumup2 = sum(tmp2*exp(sll_p_i1*ltau*k/eps))/cmplx(ntau,0.0,f64)

!---------------end time solve-------------------------

      eta1 = real(sumup1)
      eta2 = real(sumup2)

      call apply_bc()

      eta1feet(i,j)=eta1
      eta2feet(i,j)=eta2

    enddo
  enddo

  call sll_s_deposit_value_2d(eta1feet,eta2feet,spl_2d,fh_fsl)

  t = real(step,f64)*k/eps

  print"('Step =', i6, ' Time = ', g15.3)", step, t
  f0=fh_fsl(1:n,1:n)
  call solver%interp(f0,t,n,fvr)
  call sll_o_gnuplot_2d(n, x1, n, x2, fvr, 'fh', step, ierr)
  call sll_s_xdmf_rect2d_nodes( 'fh', fvr, 'fh', x1(1:n), x2(1:n), &
                                'HDF5', step) 
enddo

call sll_s_fft_free(fw_fft)
call sll_s_fft_free(bw_fft)
call solver%free()

call cpu_time(tend)
print"('CPU time = ', g15.3)", tend - tstart

error = 0.0
open(newunit = ref_id, file='fh.ref')
do i=1,n
  do j=1,n
    read(ref_id,*) m, l, fdum
    error = error + abs(fdum-sngl(fvr(i,j)))
  enddo
  read(ref_id,*)
enddo
close(ref_id)
print"('Error = ', g15.3)", error / real(n*n,f32)

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
  if (abs(eta1) < 1.0d-12) eta1 = 0.0_f64
  if (abs(eta2) < 1.0d-12) eta2 = 0.0_f64

end subroutine apply_bc

function fct1( tau, ntau, xi1, xi2, gn )

  sll_comp64 :: fct1
  sll_real64 :: tau
  sll_int32  :: ntau
  sll_comp64 :: xi1
  sll_comp64 :: xi2
  sll_comp64 :: gn

  sll_comp64 :: ctau

  ctau = cmplx(tau,0.0,f64)

  fct1 = ( - cos(2.0*ctau)**2 * ( cmplx(0.5,0.0,f64)*sin(2.0*ctau)*xi1 &
           + sin(ctau)**2*xi2) - sin(ctau)*gn)/cmplx(ntau,0.0,f64)

end function fct1

function fct2( tau, ntau, xi1, xi2, gn )

  sll_comp64 :: fct2
  sll_real64 :: tau
  sll_int32  :: ntau
  sll_comp64 :: xi1
  sll_comp64 :: xi2
  sll_comp64 :: gn

  sll_comp64 :: ctau

  ctau = cmplx(tau,0.0,f64)

  fct2 = (  cos(2.0*ctau)**2 * ( cos(ctau)**2*xi1 &
          + cmplx(0.5,0.0,f64)*sin(2.0*ctau)*xi2)  &
          + cos(ctau)*gn )/cmplx(ntau,0.0,f64)

end function fct2

function rfct1( tau, ntau, xi1, xi2, gn )

  sll_comp64 :: rfct1
  sll_real64 :: tau
  sll_int32  :: ntau
  sll_real64 :: xi1
  sll_real64 :: xi2
  sll_real64 :: gn
  sll_real64 :: tmp

  tmp = ( - cos(2.0*tau)**2 * ( 0.5_f64*sin(2.0*tau)*xi1 + sin(tau)**2*xi2) &
          - sin(tau)*gn)/real(ntau,f64)

  rfct1 = cmplx(tmp,0.0,f64)

end function rfct1

function rfct2( tau, ntau, xi1, xi2, gn )

  sll_comp64 :: rfct2
  sll_real64 :: tau
  sll_int32  :: ntau
  sll_real64 :: xi1
  sll_real64 :: xi2
  sll_real64 :: gn
  sll_real64 :: tmp

  tmp = (  cos(2.0*tau)**2 * ( cos(tau)**2*xi1 + 0.5_f64*sin(2.0*tau)*xi2)  &
         + cos(tau)*gn )/real(ntau,f64)

  rfct2 = cmplx(tmp,0.0,f64)

end function rfct2

end program test_deposit_cubic_splines

