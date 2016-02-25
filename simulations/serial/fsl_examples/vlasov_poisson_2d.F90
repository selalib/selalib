
program test_deposit_cubic_splines
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

use sll_m_cubic_splines
use sll_m_constants
use sll_m_boundary_condition_descriptors
use sll_m_fft
use sll_m_gnuplot
use sll_m_xdmf
use sll_m_cubic_splines

!$ use omp_lib

implicit none

sll_comp64, parameter :: sll_p_i0   = (0.0_f64, 0.0_f64)
sll_int32,  parameter :: n          = 128
sll_int32,  parameter :: ntau       = 32
sll_real64, parameter :: final_time = 0.4_f64

type(sll_t_fft) :: fw_fft
type(sll_t_fft) :: bw_fft
type(sll_t_cubic_spline_1d), pointer :: spl_1d
type(sll_t_cubic_spline_2d), pointer :: spl_2d

sll_int32  :: step,nb_step
sll_int32  :: i,j,l
sll_int32  :: bc1_type,bc2_type,err
sll_real64 :: delta_eta1
sll_real64 :: eta1_min
sll_real64 :: eta1_max
sll_real64 :: eps
sll_real64 :: eta1,  eta2

sll_real64, dimension(:,:), allocatable :: eta1feet
sll_real64, dimension(:,:), allocatable :: eta2feet

sll_comp64, allocatable :: ftilde1(:,:,:)
sll_comp64, allocatable :: ftilde2(:,:,:)
sll_comp64, allocatable :: F1     (:)
sll_comp64, allocatable :: F2     (:)
sll_comp64, allocatable :: tmp1   (:)
sll_comp64, allocatable :: tmp2   (:)
sll_comp64, allocatable :: tmp1_f (:)
sll_comp64, allocatable :: tmp2_f (:)
sll_comp64, allocatable :: ltau   (:)
sll_comp64, allocatable :: tmp    (:)
sll_comp64, allocatable :: sum0   (:)
sll_comp64, allocatable :: vctmp  (:)
sll_comp64, allocatable :: uctmp  (:)

sll_real64, allocatable :: gn     (:,:,:)
sll_real64, allocatable :: gnr    (:,:,:)
sll_real64, allocatable :: gnt    (:,:,:)
sll_real64, allocatable :: xi1    (:,:,:)
sll_real64, allocatable :: xi2    (:,:,:)
sll_real64, allocatable :: ftv    (:,:)
sll_real64, allocatable :: ftr    (:,:)
sll_real64, allocatable :: ft1    (:,:)
sll_real64, allocatable :: ft2    (:,:)
sll_real64, allocatable :: fh_fsl (:,:)
sll_real64, allocatable :: Ens    (:,:)
sll_real64, allocatable :: Enr    (:,:)
sll_real64, allocatable :: Ent    (:,:)
sll_real64, allocatable :: fvr    (:,:)
sll_real64, allocatable :: x      (:,:)
sll_real64, allocatable :: x1     (:)
sll_real64, allocatable :: x2     (:)
sll_real64, allocatable :: tau    (:)
sll_real64, allocatable :: r      (:)
sll_real64, allocatable :: v      (:)

sll_comp64 :: dtgn
sll_real64 :: k, h
sll_comp64 :: sumup1
sll_comp64 :: sumup2
sll_real64 :: s1, s2, s3
sll_int32  :: ref_id
sll_real32 :: fdum, error
sll_real64 :: tstart, tend
sll_real64 :: ctau, stau, csq

type :: fsl_solver

  sll_comp64, allocatable              :: lx(:)
  type(sll_t_fft)                      :: fw
  type(sll_t_fft)                      :: bw
  sll_real64                           :: L
  sll_real64, allocatable              :: r(:)
  sll_int32                            :: num_cells
  type(sll_t_cubic_spline_2d), pointer :: spl_2d

end type fsl_solver

type(fsl_solver) :: solver


! ---- * Parameters * ----

! mesh type : cartesian
! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
! BC        : periodic-periodic

eta1_min = -4.0_f64
eta1_max =  4.0_f64

! --- Space and time parameters --


! ---- * Construction of the mesh * ----
bc1_type = SLL_P_PERIODIC
bc2_type = SLL_P_PERIODIC

! ---- * Time and space steps * ----

! space steps
delta_eta1 = (eta1_max-eta1_min)/real(n,f64)

! time step and number of steps
k       = 0.05d0  !T*delta_eta1
nb_step = floor(final_time/k)
eps     = 1.00d0
h       = 2.0d0 * sll_p_pi/ real(ntau,f64)

! ---- * Messages * ----

print *,'# N   = ',n
print *,'# k   = ',k
print *,'# T   = ',final_time
print *,'# eps = ',eps

call cpu_time(tstart)


! allocations of the arrays
SLL_ALLOCATE(ftilde1(0:ntau-1,n+1,n+1),err)
SLL_ALLOCATE(ftilde2(0:ntau-1,n+1,n+1),err)
SLL_ALLOCATE(x1(n+1),err)
SLL_ALLOCATE(x2(n+1),err)
SLL_ALLOCATE(fh_fsl(n+1,n+1),err)
SLL_ALLOCATE(Ens(n+1,0:ntau-1),err)
SLL_ALLOCATE(Enr(n+1,0:ntau-1),err)
SLL_ALLOCATE(Ent(n+1,0:ntau-1),err)
SLL_ALLOCATE(gn(0:ntau-1,n+1,n+1),err)
SLL_ALLOCATE(gnr(0:ntau-1,n+1,n+1),err)
SLL_ALLOCATE(gnt(0:ntau-1,n+1,n+1),err)
SLL_ALLOCATE(F1(0:ntau-1),err)
SLL_ALLOCATE(F2(0:ntau-1),err)
SLL_ALLOCATE(tmp1(0:ntau-1),err)
SLL_ALLOCATE(tmp2(0:ntau-1),err)
SLL_ALLOCATE(tmp1_f(0:ntau-1),err)
SLL_ALLOCATE(tmp2_f(0:ntau-1),err)
SLL_ALLOCATE(tau(0:ntau-1),err)
SLL_ALLOCATE(ltau(0:ntau-1),err)
SLL_ALLOCATE(fvr(1:n,1:n),err)
SLL_ALLOCATE(x(n+1,n+1),err)
SLL_ALLOCATE(tmp(n),err)
SLL_ALLOCATE(sum0(n),err)
SLL_ALLOCATE(r(n),err)
SLL_ALLOCATE(xi1(0:ntau-1,n+1,n+1),err)
SLL_ALLOCATE(xi2(0:ntau-1,n+1,n+1),err)
SLL_ALLOCATE(ftv(n+1,n+1),err)
SLL_ALLOCATE(ftr(n+1,n+1),err)
SLL_ALLOCATE(vctmp(n),err)
SLL_ALLOCATE(uctmp(n),err)
SLL_ALLOCATE(v(n+1),err)
SLL_ALLOCATE(ft1(n,n),err)
SLL_ALLOCATE(ft2(n,n),err)
SLL_ALLOCATE(eta1feet(n+1,n+1), err)
SLL_ALLOCATE(eta2feet(n+1,n+1), err)

!$OMP PARALLEL NUM_THREADS(1)                              & 
!$OMP DEFAULT(SHARED)                                      &
!$OMP PRIVATE(spl_1d, spl_2d, solver, bw_fft, fw_fft)
 
!$OMP CRITICAL
!$ PRINT*, OMP_GET_NUM_THREADS(), OMP_GET_THREAD_NUM()
!$OMP END CRITICAL

call init_fsl_solver( solver, eta1_min, eta1_max, n)
call sll_s_fft_init_c2c_1d(fw_fft ,ntau, tmp1, tmp1_f, sll_p_fft_forward)
call sll_s_fft_init_c2c_1d(bw_fft ,ntau, tmp1_f, tmp1, sll_p_fft_backward)

spl_1d => sll_f_new_cubic_spline_1d( n+1, eta1_min, eta1_max, sll_p_periodic )

spl_2d => sll_f_new_cubic_spline_2d( n+1,      &
                                     n+1,      &
                                     eta1_min, &
                                     eta1_max, &
                                     eta1_min, &
                                     eta1_max, &
                                     bc1_type, &
                                     bc2_type)

! ---- * Initializations * ----

Ens = 0.0_f64
Enr = 0.0_f64
Ent = 0.0_f64
ftv = 0.0_f64
ftr = 0.0_f64

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

do l=0,ntau-1
  tau(l)=real(l,f64)*h
enddo

ltau = [(cmplx(0.,-l,f64),l=0,ntau/2-1),(cmplx(0,-l,f64),l=-ntau/2,-1)]

do step=1,nb_step !-------- * Evolution in time * ---------


  call sll_s_compute_cubic_spline_2d(fh_fsl,spl_2d)

  r        = solver%r
  r(n/2+1) = 1.0_f64
  do l=0,ntau-1
  
    call fsl_interp_2d(solver,fh_fsl,tau(l),n,fvr)
    do i=1,n
      tmp = cmplx(fvr(i,:),0.,f64)
      call sll_s_fft_exec_c2c_1d(solver%fw, tmp, tmp)
      sum0(i)=tmp(1)*cmplx(2.0*solver%L*solver%r(i)/n,0.0,f64) 
    enddo
    call sll_s_fft_exec_c2c_1d(solver%fw, sum0, tmp)
    
    tmp(2:n)=tmp(2:n)/solver%lx(2:n)
    tmp(1)=cmplx(0.0d0,0.0d0,kind=f64)
    call sll_s_fft_exec_c2c_1d(solver%bw, tmp, tmp)
  
    Ens(1:n,l)=real(tmp-tmp(n/2+1))/r
    Enr(1:n,l)=real(sum0-Ens(1:n,l))/r
  
  enddo
  
  do j=1,n
  
    vctmp = cmplx(fh_fsl(j,1:n),0.,f64)
    uctmp = cmplx(fh_fsl(1:n,j),0.,f64)
  
    call sll_s_fft_exec_c2c_1d(solver%fw, vctmp, vctmp)
    call sll_s_fft_exec_c2c_1d(solver%fw, uctmp, uctmp)
    
    uctmp = uctmp/cmplx(n**2,0.,f64)*solver%lx
    vctmp = vctmp/cmplx(n**2,0.,f64)*solver%lx
   
    call sll_s_fft_exec_c2c_1d(solver%bw, vctmp, tmp)
    ftv(j,1:n)=real(tmp)  !\partial_\x1 f_tilde(\xi1,\xi2)
    call sll_s_fft_exec_c2c_1d(solver%bw, uctmp, tmp)
    ftr(1:n,j)=real(tmp)  !\partial_\x2 f_tilde(\xi1,\xi2)
  
  enddo
  
  v(1:n)=solver%r
  v(n+1)=solver%L
  do l=0,ntau-1
    ctau = cos(tau(l))
    stau = sin(tau(l))
    do j=1,n+1
      do i=1,n+1
        xi1(l,i,j)= v(i)*ctau-v(j)*stau
        xi2(l,i,j)= v(i)*stau+v(j)*ctau
        x(i,j)    = v(i) 
      end do
    end do

    call fsl_interp_1d(spl_1d, n, x, ens(:,l), gn(l,:,:))

  enddo

  do l=0,ntau-1
  
    call fsl_interp_2d(solver,ftv,tau(l),n,ft1)
    call fsl_interp_2d(solver,ftr,tau(l),n,ft2)
  
    ctau = cos(tau(l))
    stau = sin(tau(l))
    csq  = cos(2.0_f64*tau(l))**2

    do i=1,n

      do j=1,n
        s1 =  xi1(l,i,j)*ctau+xi2(l,i,j)*stau
        s2 = -stau*ft1(i,j)+ctau*ft2(i,j)
        s3 = (csq * s1 + gn(l,i,j)) * s2
        vctmp(j) = cmplx(s3,0.,f64)
      enddo
  
      call sll_s_fft_exec_c2c_1d(solver%fw, vctmp, tmp)
      sum0(i)=tmp(1)*cmplx(2.0d0*solver%L*solver%r(i)/n,0.,f64)
  
    enddo
  
    call sll_s_fft_exec_c2c_1d(solver%fw, sum0, tmp)
  
    tmp(2:n) = tmp(2:n)/solver%lx(2:n)
    tmp(1)   = cmplx(0.0d0,0.0d0,kind=f64)
  
    call sll_s_fft_exec_c2c_1d(solver%bw, tmp, tmp)
  
    Ent(1:n,l)=real(tmp-tmp(n/2+1))/r
  
  enddo


  do l=0,Ntau-1
  
    ctau = cos(tau(l))
    stau = sin(tau(l))
    do j=1,n+1
      do i=1,n+1
        x(i,j)=ctau*x1(i)+stau*x1(j)
      end do
    end do

    call fsl_interp_1d(spl_1d, n, x, ens(:,l), gn(l,:,:))
    call fsl_interp_1d(spl_1d, n, x, enr(:,l), gnr(l,:,:))
    call fsl_interp_1d(spl_1d, n, x, ent(:,l), gnt(l,:,:))
  
  enddo

  do j=1,n+1
    do i=1,n+1

      !------------for 1st order correction------------------
      do l=0,ntau-1
        F1(l) = rfct1(tau(l), ntau, x1(i), x2(j), gn(l,i,j))
        F2(l) = rfct2(tau(l), ntau, x1(i), x2(j), gn(l,i,j))
      enddo

      call sll_s_fft_exec_c2c_1d(fw_fft, F1, tmp1_f)
      call sll_s_fft_exec_c2c_1d(fw_fft, F2, tmp2_f)

      ftilde1(:,i,j) = tmp1_f
      ftilde2(:,i,j) = tmp2_f

      tmp1_f(1:ntau-1) = - tmp1_f(1:ntau-1) / ltau(1:ntau-1)
      tmp2_f(1:ntau-1) = - tmp2_f(1:ntau-1) / ltau(1:ntau-1)

      tmp1_f(0) = sll_p_i0
      tmp2_f(0) = sll_p_i0

      call sll_s_fft_exec_c2c_1d(bw_fft, tmp1_f, F1)
      call sll_s_fft_exec_c2c_1d(bw_fft, tmp2_f, F2)

      F1 = F1 - sum(tmp1_f)
      F2 = F2 - sum(tmp2_f)

      xi1(:,i,j) = x1(i) + eps*real(F1)
      xi2(:,i,j) = x2(j) + eps*real(F2)

    enddo
  enddo

  do l=0,Ntau-1
  
    ctau = cos(tau(l))
    stau = sin(tau(l))
    do j=1,n+1
      do i=1,n+1
        x(i,j)=ctau*xi1(l,i,j)+stau*xi2(l,i,j)
      end do
    end do

    call fsl_interp_1d(spl_1d, n, x, ens(:,l), gn(l,:,:))
    call fsl_interp_1d(spl_1d, n, x, enr(:,l), gnr(l,:,:))
    call fsl_interp_1d(spl_1d, n, x, ent(:,l), gnt(l,:,:))

  enddo
  
  do j=1,n+1
    do i=1,n+1
      !------------for 2nd order correction------------------
      do l=0,ntau-1

        F1(l)= rfct1(tau(l), 1, xi1(l,i,j), xi2(l,i,j), gn(l,i,j))
        F2(l)= rfct2(tau(l), 1, xi1(l,i,j), xi2(l,i,j), gn(l,i,j))

        dtgn=  cmplx(gnt(l,i,j),0.0,f64) &
             + cmplx(gnr(l,i,j),0.0,f64) &
             *(cmplx(cos(tau(l)),0.0,f64)*ftilde1(0,i,j) &
             + cmplx(sin(tau(l)),0.0,f64)*ftilde2(0,i,j))

        tmp1(l)= cfct1(tau(l), ntau, ftilde1(0,i,j), ftilde2(0,i,j), dtgn)
        tmp2(l)= cfct2(tau(l), ntau, ftilde1(0,i,j), ftilde2(0,i,j), dtgn)

      enddo

      call sll_s_fft_exec_c2c_1d(fw_fft, tmp1, tmp1_f)
      call sll_s_fft_exec_c2c_1d(fw_fft, tmp2, tmp2_f)

      tmp1_f(1:ntau-1) = - tmp1_f(1:ntau-1)/ltau(1:ntau-1)
      tmp2_f(1:ntau-1) = - tmp2_f(1:ntau-1)/ltau(1:ntau-1)

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

      tmp1_f(1:ntau-1) = -tmp1_f(1:ntau-1)/ltau(1:ntau-1)
      tmp2_f(1:ntau-1) = -tmp2_f(1:ntau-1)/ltau(1:ntau-1)

      tmp1_f(0) = sll_p_i0
      tmp2_f(0) = sll_p_i0

      call sll_s_fft_exec_c2c_1d(bw_fft, tmp1_f, tmp1)
      call sll_s_fft_exec_c2c_1d(bw_fft, tmp2_f, tmp2)

      tmp1 = tmp1 - sum(tmp1_f)
      tmp2 = tmp2 - sum(tmp2_f)

      xi1(:,i,j) = x1(i) + eps*real(tmp1)
      xi2(:,i,j) = x2(j) + eps*real(tmp2)

    enddo
  enddo

  do l=0,Ntau-1
    ctau = cos(tau(l))
    stau = sin(tau(l))
    do j=1,n+1
      do i=1,n+1
        x(i,j)=ctau*xi1(l,i,j)+stau*xi2(l,i,j)
      end do
    end do
    call fsl_interp_1d(spl_1d, n, x, ens(:,l), gn(l,:,:))
  enddo

  do j=1,n+1
    do i=1,n+1

      do l=0,ntau-1
        F1(l) = rfct1(tau(l),1,xi1(l,i,j),xi2(l,i,j),gn(l,i,j))
        F2(l) = rfct2(tau(l),1,xi1(l,i,j),xi2(l,i,j),gn(l,i,j))
      enddo

      tmp1 = xi1(:,i,j) + 0.5_f64*k*F1
      tmp2 = xi2(:,i,j) + 0.5_f64*k*F2

      call sll_s_fft_exec_c2c_1d(fw_fft, tmp1, tmp1_f)
      call sll_s_fft_exec_c2c_1d(fw_fft, tmp2, tmp2_f)

      tmp1_f = tmp1_f/(1.0d0-0.5_f64*k*ltau/eps)/cmplx(ntau,0.0,f64)
      tmp2_f = tmp2_f/(1.0d0-0.5_f64*k*ltau/eps)/cmplx(ntau,0.0,f64)

      call sll_s_fft_exec_c2c_1d(bw_fft, tmp1_f, tmp1)
      call sll_s_fft_exec_c2c_1d(bw_fft, tmp2_f, tmp2)

      ftilde1(:,i,j) = tmp1
      ftilde2(:,i,j) = tmp2

      !----------Insert half step evaluation--------

      eta1 = real(sum(tmp1_f*exp(-0.5_f64*ltau*k/eps)))
      eta2 = real(sum(tmp2_f*exp(-0.5_f64*ltau*k/eps)))

      call apply_bc()

      eta1feet(i,j)=eta1
      eta2feet(i,j)=eta2

    enddo
  enddo

  call sll_s_deposit_value_2d(eta1feet,eta2feet,spl_2d,fh_fsl) !function value at the half time

  r        = solver%r
  r(n/2+1) = 1.0_f64

  do l=0,ntau-1
  
    call fsl_interp_2d(solver,fh_fsl,tau(l),n,fvr)
  
    do i=1,n
      tmp = cmplx(fvr(i,:),0.0,f64)
      call sll_s_fft_exec_c2c_1d(solver%fw, tmp, tmp)
      sum0(i)=tmp(1)*solver%r(i) !r*int_R fdv
    enddo
  
    sum0 = sum0 * 2.0_f64 * solver%L / cmplx(n,0.0, f64)
  
    call sll_s_fft_exec_c2c_1d(solver%fw, sum0, tmp)
  
    tmp(1)=cmplx(0.0d0,0.0d0,kind=f64)
  
    tmp(2:n) = tmp(2:n) / solver%lx(2:n)
  
    call sll_s_fft_exec_c2c_1d(solver%bw, tmp, tmp)
    
    Ens(1:n,l)=real(tmp-tmp(n/2+1),f64)/r(:)
   
  enddo

  !------End evaluation and continue 2nd solver-------------

  do l=0,Ntau-1
    ctau = cos(tau(l))
    stau = sin(tau(l))
    do j=1,n+1
      do i=1,n+1
        x(i,j)=ctau*real(ftilde1(l,i,j))+stau*real(ftilde2(l,i,j))
      end do
    end do
    call fsl_interp_1d(spl_1d, n, x, ens(:,l), gn(l,:,:))
  enddo

  do j=1,n+1
    do i=1,n+1

      do l=0,ntau-1
        F1(l) = cfct1(tau(l), 1, ftilde1(l,i,j), ftilde2(l,i,j), cmplx(gn(l,i,j),0.0,f64))
        F2(l) = cfct2(tau(l), 1, ftilde1(l,i,j), ftilde2(l,i,j), cmplx(gn(l,i,j),0.0,f64))
      enddo

      call sll_s_fft_exec_c2c_1d(fw_fft, F1, tmp1_f)
      call sll_s_fft_exec_c2c_1d(fw_fft, F2, tmp2_f)

      tmp1 = cmplx(xi1(:,i,j),0.0,f64)
      tmp2 = cmplx(xi2(:,i,j),0.0,f64)

      call sll_s_fft_exec_c2c_1d(fw_fft, tmp1, F1)
      call sll_s_fft_exec_c2c_1d(fw_fft, tmp2, F2)

      tmp1 = (F1*(1.0_f64+0.5_f64*k/eps*ltau)+k*tmp1_f) &
                /(1.0_f64-0.5_f64*k*ltau/eps)
      tmp2 = (F2*(1.0_f64+0.5_f64*k/eps*ltau)+k*tmp2_f) &
                /(1.0_f64-0.5_f64*k*ltau/eps)

      sumup1 = sum(tmp1*exp(-ltau*k/eps))/cmplx(ntau,0.0,f64)
      sumup2 = sum(tmp2*exp(-ltau*k/eps))/cmplx(ntau,0.0,f64)

      eta1 = real(sumup1)
      eta2 = real(sumup2)

      call apply_bc()

      eta1feet(i,j)=eta1
      eta2feet(i,j)=eta2

    enddo
  enddo

  call sll_s_deposit_value_2d(eta1feet,eta2feet,spl_2d,fh_fsl)

  print"('Step =', i6, ' Time = ', g15.3)", step, real(step*k/eps)
  !call sll_o_gnuplot_2d(n, x1, n, x2, fvr, 'fh', step, ierr)
  !call sll_s_xdmf_rect2d_nodes( 'fh', fvr, 'fh', x1(1:n), x2(1:n), &
  !                              'HDF5', step) 

  tau = tau + k/eps

enddo

!---------------end time solve-------------------------

call sll_o_delete(spl_1d)
call sll_o_delete(spl_2d)
call sll_s_fft_free(fw_fft)
call sll_s_fft_free(bw_fft)
call fsl_interp_2d(solver,fh_fsl,real(nb_step*k/eps,f64),n,fvr)
call free_fsl_solver(solver)

!$OMP END PARALLEL

call cpu_time(tend)
print"('CPU time = ', g15.3)", tend - tstart


error = 0.0
open(newunit = ref_id, file='fh.ref')
do i=1,n
  do j=1,n
    read(ref_id,*) err, l, fdum
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

function cfct1( tau, ntau, xi1, xi2, gn )

  sll_comp64 :: cfct1
  sll_real64 :: tau
  sll_int32  :: ntau
  sll_comp64 :: xi1
  sll_comp64 :: xi2
  sll_comp64 :: gn

  sll_comp64 :: ctau

  ctau = cmplx(tau,0.0,f64)

  cfct1 = ( - cos(2.0*ctau)**2 * ( cmplx(0.5,0.0,f64)*sin(2.0*ctau)*xi1 &
            + sin(ctau)**2*xi2) - sin(ctau)*gn)/cmplx(ntau,0.0,f64)

end function cfct1

function cfct2( tau, ntau, xi1, xi2, gn )

  sll_comp64 :: cfct2
  sll_real64 :: tau
  sll_int32  :: ntau
  sll_comp64 :: xi1
  sll_comp64 :: xi2
  sll_comp64 :: gn

  sll_comp64 :: ctau

  ctau = cmplx(tau,0.0,f64)

  cfct2 = (  cos(2.0*ctau)**2 * ( cos(ctau)**2*xi1 &
           + cmplx(0.5,0.0,f64)*sin(2.0*ctau)*xi2)  &
           + cos(ctau)*gn )/cmplx(ntau,0.0,f64)

end function cfct2

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

subroutine init_fsl_solver( self, xmin, xmax, num_cells)

type(fsl_solver)        :: self
sll_int32 , intent(in)  :: num_cells
sll_real64, intent(in)  :: xmin
sll_real64, intent(in)  :: xmax

sll_int32               :: i, j, m
sll_comp64, allocatable :: tmp(:)

self%num_cells = num_cells
m = num_cells/2
allocate(self%lx(num_cells))
self%lx = [ (cmplx(j,0.,f64), j=0,m-1), (cmplx(j,0.,f64), j=-m,-1 ) ]
self%lx = self%lx * 2.0d0*sll_p_pi/(xmax-xmin) * sll_p_i1 * num_cells

self%L = 4.0_f64

allocate(tmp(num_cells))
call sll_s_fft_init_c2c_1d(self%fw, num_cells, tmp, tmp, sll_p_fft_forward)
call sll_s_fft_init_c2c_1d(self%bw, num_cells, tmp, tmp, sll_p_fft_backward)
deallocate(tmp)

self%spl_2d => sll_f_new_cubic_spline_2D(num_cells+1, &
                                         num_cells+1, &
                                                xmin, &
                                                xmax, &
                                                xmin, &
                                                xmax, &
                                      SLL_P_PERIODIC, &
                                      SLL_P_PERIODIC)

allocate(self%r(num_cells))
do i = 1, num_cells
  self%r(i) = xmin + (i-1) * (xmax-xmin)/real(num_cells,f64)
end do

end subroutine init_fsl_solver

subroutine free_fsl_solver( self )

type(fsl_solver)   :: self

deallocate(self%lx)
call sll_s_fft_free(self%fw)
call sll_s_fft_free(self%bw)
call sll_o_delete(self%spl_2d)

end subroutine free_fsl_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsl_interp_2d(self, fh_fsl,t,num_cells,fvr)

type(fsl_solver)            :: self
sll_int32,  intent(in)      :: num_cells
sll_real64, intent(in)      :: fh_fsl(:,:)
sll_real64, intent(in)      :: t
sll_real64, intent(inout)   :: fvr(num_cells,num_cells)

sll_real64                  :: ct, st
sll_real64                  :: x, y
sll_real64                  :: xj, yj
sll_int32                   :: i, j

call sll_s_compute_cubic_spline_2d(fh_fsl, self%spl_2d)

ct = cos(t)
st = sin(t)
do j=1,num_cells
  xj = st*self%r(j)
  yj = ct*self%r(j)
  do i=1,num_cells
    x = ct*self%r(i)-xj
    y = st*self%r(i)+yj
    if (abs(x)<self%L .and. abs(y)<self%L) then
      fvr(i,j)=sll_f_interpolate_value_2d(x,y,self%spl_2d)
    else
      fvr(i,j)=0.0_f64
    endif
  enddo
enddo

end subroutine fsl_interp_2d

subroutine fsl_interp_1d(spl, n, x, e, g)

type(sll_t_cubic_spline_1d), pointer :: spl
sll_int32,  intent(in)               :: n
sll_real64, intent(in)               :: x(:,:)
sll_real64, intent(in)               :: e(:)
sll_real64, intent(out)              :: g(:,:)

sll_int32 :: i, j

call sll_s_compute_cubic_spline_1d(e, spl)
do j=1,n+1
  do i=1,n+1
    if (eta1_min < x(i,j) .and. x(i,j) < eta1_max ) then
      g(i,j) = sll_f_interpolate_from_interpolant_value(x(i,j), spl)
    else
      g(i,j) = 0.0_f64
    end if
  enddo
enddo

end subroutine fsl_interp_1d

end program test_deposit_cubic_splines

