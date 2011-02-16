!**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!     Chantal Passeron (chantal.passeron@cea.fr)
!     Guillaume Latu (guillaume.latu@cea.fr)
!     Xavier Garbet (xavier.garbet@cea.fr)
!     Philippe Ghendrih (philippe.ghendrih@cea.fr)
!     Yanick Sarazin (yanick.sarazin@cea.fr)
!  
!  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
!  is a 5D gyrokinetic global full-f code for simulating 
!  the plasma turbulence in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
      
!------------------------------------------------------------
! file : poisson.f90
! date : 08/04/2001
! . used to solve the Poisson equation
!    -[d^2(Phi)/dr^2+(1/r)*dPhi/dr+(1/r^2)*d^2(Phi)/dth^2] 
!        = rho(phi,r,th)
!   where Phi is considered with periodic conditions in
!   theta direction and Dirichlet conditions in r direction
!--------------------------------------------------------------
module poisson_class
  use prec_const
  use OMPutils_module, only : omp1dvector, omp1cvector, omp2cvector
  use geometry_class
  use init_profile_class
  use init_magnetic_class
  use fdistribu5d_class
  use fft_module
  use mem_alloc_module
  use fft_NRF90_module
  use integration_module
  use clock_module
  
  implicit none
      
  private
  public :: new_poisson, del_poisson, solve_poisson, &
    compute_nbparticles
      
  type, public :: poisson
     real(RKIND)   , dimension(:,:,:), pointer :: Phi       
     real(RKIND)   , dimension(:)    , pointer :: kth, kphi
     !-> vector used for the LU factorisation 
     !->  of the system for (m,n) = (0,0)
     ! ---> diagonal terms of matrix
     complex(CKIND), dimension(:,:)  , pointer :: diag_m0n0    
     ! ---> lower diagonal terms of matrix
     complex(CKIND), dimension(:,:)  , pointer :: ldiag_m0n0
     ! ---> upper diagonal terms of matrix
     complex(CKIND), dimension(:,:)  , pointer :: udiag_m0n0   
     ! -> vector used for the LU factorisation 
     ! ->  of the system for (m,n) != (0,0)
     ! ---> diagonal terms of matrix
     complex(CKIND), dimension(:,:)  , pointer :: diag
     ! ---> lower diagonal terms of matrix
     complex(CKIND), dimension(:,:)  , pointer :: ldiag
     ! ---> upper diagonal terms of matrix        
     complex(CKIND), dimension(:,:)  , pointer :: udiag
  end type poisson
      
  type(omp1cvector), dimension (:), pointer :: omp_rhs
  type(omp2cvector), dimension (:), pointer :: omp_data_i
  type(omp2cvector), dimension (:), pointer :: omp_tmp
      
  real(RKIND), dimension(:), pointer, public :: Phi00_diag
      
  include "Bstar_inline.h"
  !******************************
  contains
  !******************************
#include "Bstar_inline.f90"
      
  !----------------------------------------------------- 
  ! constructor for the poisson type
  !-----------------------------------------------------   
  subroutine new_poisson(pthis,geom,init_prof)
    use globals, only : memory_test
    type(poisson)     , intent(out) :: pthis
    type(geometry)    , intent(in)  :: geom 
    type(init_profile), intent(in)  :: init_prof
      
    !*** local variables ***
    integer     :: i, j, k           ! loop index
    integer     :: N2m1, N3m1        ! dimension in theta and phi
    integer     :: min_m, max_m, min_n, max_n
    integer     :: tid
      
    !*** set local variables ***
    N2m1  = geom%Ntheta           !size(geom%thetag)-1
    N3m1  = geom%Nphi             !size(geom%phig)-1
      
    !*** memory allocation ***
    call array_allocation(pthis,geom)
      
    if (.not.memory_test) then
      !*** compute kth for the FFT in theta direction ***
      call NR_fft_vector(N2m1,geom%dtheta,pthis%kth)
      
      !*** compute kphi for the FFT in phi direction ***
      call NR_fft_vector(N3m1,geom%dphi,pthis%kphi)
      
      !*** assemble and factorize poisson matrix for all k ***
      call LU_factorisation(pthis,init_prof,geom)
    end if
  end subroutine new_poisson
      
  !----------------------------------------------------- 
  ! Memory allocation (used in the constructor)
  !-----------------------------------------------------   
  subroutine array_allocation(pthis,geom)
    use globals, only : Nbthread, plocal_id, cible
    type(poisson) , intent(inout) :: pthis
    type(geometry), intent(in)    :: geom  
      
    integer :: N1, N2, N3, N2m1, N3m1
    integer :: tid
      
    N1   = geom%Nr+1       !size(geom%rg)
    N2   = geom%Ntheta+1   !size(geom%thetag)
    N2m1 = N2-1
    N3   = geom%Nphi+1     !size(geom%phig)
    N3m1 = N3-1
      
    call glob_allocate(pthis%kth,1,N2m1,'pthis%kth')
    call glob_allocate(pthis%kphi,1,N3m1,'pthis%kphi')
      
    ! -> array for the solution : Phi 
    call glob_allocate(pthis%Phi,0,geom%Nr,0,geom%Ntheta, &
      0,geom%Nphi,'pthis%Phi')
      
    ! -> arrays for LU factorisation and 
    ! ->   solving with (m,n) = (0,0) 
    call glob_allocate(pthis%diag_m0n0,1,N2m1,1,N1-2, &
      'pthis%diag_m0n0')
    call glob_allocate(pthis%ldiag_m0n0,1,N2m1,1,N1-3, &
      'pthis%ldiag_m0n0')
    call glob_allocate(pthis%udiag_m0n0,1,N2m1,1,N1-3, &
      'pthis%udiag_m0n0')
      
    ! -> arrays for LU factorisation and solving
    !    with (m,n) != (0,0)
    call glob_allocate(pthis%diag,1,N2m1,1,N1-2,'pthis%diag')
    call glob_allocate(pthis%ldiag,1,N2m1,1,N1-3,'pthis%ldiag')
    call glob_allocate(pthis%udiag,1,N2m1,1,N1-3,'pthis%udiag')
      
    !*** allocation of arrays for OMP parallelisation ***
    allocate(omp_rhs(1:Nbthread))
    allocate(omp_tmp(1:Nbthread))
    allocate(omp_data_i(1:Nbthread))
    do tid = 1,Nbthread
       call glob_allocate(omp_rhs(tid)%val,1,N1,'pthis%omp_rhs')
       call glob_allocate(omp_tmp(tid)%val,1,N3m1,1,N2m1, &
         'pthis%omp_tmp')
       call glob_allocate(omp_data_i(tid)%val,1,N2m1,1,N3m1, &
         'pthis%omp_data_i')
    enddo
      
    !*** allocation of arrays for output saving ***
    !-> for radial profile saving of the Phi (0,0) Fourier mode
    call glob_allocate(Phi00_diag,0,geom%Nr,'Phi00_diag')
  end subroutine array_allocation
      
  !----------------------------------------------------- 
  ! LU factorisation (used in the constructor)
  !-----------------------------------------------------   
  subroutine LU_factorisation(pthis,init_prof,geom)
    use globals, only : Zi, zonal_flow, coeff_pol, &
      Phi00BCr_Neumann
    type(poisson)     , intent(inout) :: pthis
    type(init_profile), intent(in)    :: init_prof
    type(geometry)    , intent(in)    :: geom
      
    integer     :: i, j, N1, N2, N2m1
    real(RKIND) :: lambda, inv_ZiTe, inv_dr2
    real(RKIND) :: ri, rip1, coef_ri, coef_rip1
    real(RKIND) :: alpha_ri, alpha_rip1
    real(RKIND) :: alpha_r1, coef_r1, lr1, coef1_lr1, coef1_lr1_00
    real(RKIND) :: rN1m2, alpha_rN1m2
    real(RKIND) :: coef_rN1m2, crN1m2, coef1_crN1m2_00
    real(RKIND) :: ldiag_tmp, ddiag_tmp, udiag_tmp
        
    !*** Choice between Dirichlet or Neumann conditions ***
    !***  for the mode (0,0)                            ***
    if (Phi00BCr_Neumann) then
      coef1_lr1       = 0._RKIND      
      coef1_lr1_00    = 1._RKIND   ! Neumann   at r=rmin
      coef1_crN1m2_00 = 0._RKIND   ! Dirichlet at r=rmax
    else  
      !-> Neumann conditions are imposed on the axis 
      !    in the case rmin=0, otherwise Dirichlet    
      !    conditions are imposed                     
      if (geom%rg(0).ne.0_RKIND) then
        coef1_lr1 = 0._RKIND
      else
        coef1_lr1 = 1._RKIND
      end if
      coef1_lr1_00    = 0._RKIND
      coef1_crN1m2_00 = 0._RKIND
    end if
      
    !*** Initialisation of lambda ***
    !***  lambda = 0 if zonal_flow = false ***
    !***  lambda = 1 if zonal_flow = true  ***
    if (zonal_flow) then
      lambda = 1._RKIND
    else
      lambda = 0._RKIND
    end if
      
    !*** set local variables ***
    N1   = size(geom%rg)
    N2   = size(geom%thetag)
    N2m1 = N2-1
      
    !*** assemble the system matrices ***
    inv_dr2 = 1._RKIND/(geom%dr*geom%dr)
    !-> initialisation of the lower and upper diagonal 
    !->  of the matrix
    do i = 1,N1-3
      ri         = geom%rg(i)
      alpha_ri   = 1._RKIND/ri+init_prof%dn0dr(i)/init_prof%n0(i)
      coef_ri    = alpha_ri/(2._RKIND*geom%dr)
      rip1       = geom%rg(i+1)
      alpha_rip1 = 1._RKIND/rip1 + &
        init_prof%dn0dr(i+1)/init_prof%n0(i+1)
      coef_rip1  = alpha_rip1/(2._RKIND*geom%dr)
      do j = 1,N2m1
        udiag_tmp = -1._RKIND*coeff_pol*(inv_dr2 + coef_ri)
        ldiag_tmp = -1._RKIND*coeff_pol*(inv_dr2 - coef_rip1)
        ! -> forall (m,n) != (0,0)
        pthis%udiag(j,i) = udiag_tmp
        pthis%ldiag(j,i) = ldiag_tmp
        ! -> for (m,n) = (0,0)
        pthis%udiag_m0n0(j,i) = udiag_tmp
        pthis%ldiag_m0n0(j,i) = ldiag_tmp
      end do
    end do
    !-> initialisation of the diagonal of the matrix
    do i = 1,N1-3
      ri         = geom%rg(i)
      inv_ZiTe   = 1._RKIND/(real(Zi)*init_prof%Te(i)) 
      do j = 1,N2m1
        ddiag_tmp = coeff_pol * (2._RKIND * inv_dr2 &
          + (pthis%kth(j)*pthis%kth(j))/(ri*ri))
        ! -> for (m,n) != (0,0)
        pthis%diag(j,i) = ddiag_tmp + inv_ZiTe
        ! -> for (m,n) = (0,0)
        pthis%diag_m0n0(j,i) = ddiag_tmp + &
          (1._RKIND-lambda)*inv_ZiTe
      end do
    end do
      
    !*** boundary conditions at r = rmin ***
    !***  (first matrix row -> i = 1)    ***
    alpha_r1 = 1._RKIND/geom%rg(1) + &
      init_prof%dn0dr(1)/init_prof%n0(1)
    coef_r1  = alpha_r1/(2._RKIND*geom%dr)
    lr1      = -1._RKIND*(inv_dr2 - coef_r1)
    ! -> for (m,n) != (0,0)
    do j = 1,N2m1
      pthis%diag(j,1) = coeff_pol*(coef1_lr1*lr1) + pthis%diag(j,1)
    end do
    ! -> for (m,n) = (0,0)
    do j = 1,N2m1
      pthis%diag_m0n0(j,1) = coeff_pol*(coef1_lr1_00*lr1) + &
        pthis%diag_m0n0(j,1)
    end do
      
    !*** boundary conditions at r = rmax       ***
    !***  (last matrix row -> i = N1-2 = Nr-1) ***
    !***  case i = N1-2                        ***
    rN1m2       = geom%rg(N1-2)
    inv_ZiTe    = 1._RKIND/(real(Zi)*init_prof%Te(N1-2))
    alpha_rN1m2 = 1._RKIND/rN1m2+init_prof%dn0dr(1)/init_prof%n0(1)
    coef_rN1m2  = alpha_rN1m2/(2._RKIND*geom%dr)
    crN1m2      = -1._RKIND*(inv_dr2 + coef_rN1m2)
    do j = 1,N2m1
      ddiag_tmp = coeff_pol * (2._RKIND * inv_dr2 &
        + (pthis%kth(j)*pthis%kth(j))/(rN1m2*rN1m2))
      ! -> for (m,n) != (0,0)
      pthis%diag(j,N1-2) = ddiag_tmp + inv_ZiTe
      ! -> for (m,n) = (0,0)
      pthis%diag_m0n0(j,N1-2) = coeff_pol*coef1_crN1m2_00*crN1m2 + &
        ddiag_tmp + (1._RKIND-lambda)*inv_ZiTe
    end do
  end subroutine LU_factorisation
      
  !----------------------------------------------------- 
  ! destructor for the poisson type
  !-----------------------------------------------------   
  subroutine del_poisson(pthis)
    use globals, only : Nbthread, plocal_id, cible
    type(poisson), intent(inout) :: pthis
      
    integer :: tid
      
    call glob_deallocate(pthis%kth)
    call glob_deallocate(pthis%kphi)
    call glob_deallocate(pthis%Phi)
      
    !*** Matrix for (m,n) = (0,0) ***
    call glob_deallocate(pthis%diag_m0n0)
    call glob_deallocate(pthis%ldiag_m0n0)
    call glob_deallocate(pthis%udiag_m0n0)
      
    !*** Matrix for (m,n) != (0,0) ***
    call glob_deallocate(pthis%diag)
    call glob_deallocate(pthis%ldiag)
    call glob_deallocate(pthis%udiag)
      
    do tid = 1,Nbthread
       call glob_deallocate(omp_rhs(tid)%val)
       call glob_deallocate(omp_data_i(tid)%val)
       call glob_deallocate(omp_tmp(tid)%val)
    enddo
      
    !*** deallocation of arrays for OMP parallelisation ***
    deallocate(omp_rhs)
    deallocate(omp_data_i)
    deallocate(omp_tmp)
      
    !*** deallocation of arrays for output saving ***
    call glob_deallocate(Phi00_diag)
  end subroutine del_poisson
      
  !*****************************************************
  ! PARALLEL VERSION
  !*****************************************************
  !-----------------------------------
  ! Thomas algorithm for LU inversion
  !-----------------------------------
  subroutine thomas_c(N,ud,ld,diag,rhs,err) 
    integer                     , intent(in)  :: N
    real(RKIND)   , dimension(:), pointer     :: ud
    real(RKIND)   , dimension(:), pointer     :: ld
    real(RKIND)   , dimension(:), pointer     :: diag
    complex(CKIND), dimension(:), pointer     :: rhs
    integer                     , intent(out) :: err
      
    integer :: nm1, i, i1
      
    nm1 = N-1
    err = 0
    do i = 1,nm1
       i1 = i+1
       if (diag(i) .eq. 0._RKIND) then
          err = 1
          return
       endif
       ld(i)    = ld(i)/diag(i)
       diag(i1) = diag(i1) - ld(i)*ud(i)
       rhs(i1)  = rhs(i1) - ld(i)*rhs(i)
    end do
      
    rhs(N) = rhs(N)/diag(N)
      
    do i = nm1,1,-1
       rhs(i) = (rhs(i)-ud(i)*rhs(i+1))/diag(i)
    end do
  end subroutine thomas_c
      
  !----------------------------------------------
  ! Thomas algorithm for LU inversion
  !  (the same than before but using a temporary
  !  array 'tmpv' to keep the initial values of
  !  'ld' and 'diag')
  !  => case where 'rhs' is real 
  !----------------------------------------------
  subroutine lusolver(N,ud,ld,diag,rhs,tmpv,err) 
    integer                  , intent(in)  :: N
    real(RKIND), dimension(:), pointer     :: ud
    real(RKIND), dimension(:), pointer     :: ld
    real(RKIND), dimension(:), pointer     :: diag
    real(RKIND), dimension(:), pointer     :: rhs
    real(RKIND), dimension(:), pointer     :: tmpv
    integer                  , intent(out) :: err
      
    integer     :: nm1, i
    real(RKIND) :: denom_tmp
    
    nm1     = N-1
    err     = 0
    tmpv(1) = 1._RKIND/diag(1)
    rhs(1)  = rhs(1)*tmpv(1)
    do i = 1,nm1
      denom_tmp = (diag(i+1) - ld(i)*ud(i)*tmpv(i))
      if (abs(denom_tmp).le.1.e-12) then
        print*,'Problem denom_tmp = ',denom_tmp
        err = 1
        return
      end if
      tmpv(i+1) = 1._RKIND/denom_tmp
      rhs(i+1) = (rhs(i+1) - ld(i)*rhs(i))*tmpv(i+1)
    end do
      
    do i = nm1,1,-1
       rhs(i)= (rhs(i)-ud(i)*tmpv(i)*rhs(i+1))
    end do
  end subroutine lusolver
      
  !----------------------------------------------
  ! Thomas algorithm for LU inversion
  !  (the same than before but using a temporary
  !  array 'tmpv' to keep the initial values of
  !  'ld' and 'diag')
  !  => case where 'rhs' is complex 
  !----------------------------------------------
  subroutine lusolvec(N,ud,ld,diag,rhs,tmpv,err) 
    integer                     , intent(in)  :: N
    real(RKIND),    dimension(:), pointer     :: ud
    real(RKIND),    dimension(:), pointer     :: ld
    real(RKIND),    dimension(:), pointer     :: diag
    complex(RKIND), dimension(:), pointer     :: rhs
    real(RKIND),    dimension(:), pointer     :: tmpv
    integer                     , intent(out) :: err
      
    integer :: nm1, i
    
    nm1     = N-1
    err     = 0
    tmpv(1) = 1._RKIND/diag(1)
    rhs(1)  = rhs(1)*tmpv(1)
    do i = 1,nm1
       tmpv(i+1) = 1._RKIND/(diag(i+1) - ld(i)*ud(i)*tmpv(i))
       if (tmpv(i+1) .eq. 0._RKIND) then
          err = 1
          return
       endif
       rhs(i+1) = (rhs(i+1) - ld(i)*rhs(i))*tmpv(i+1)
    end do
    do i = nm1,1,-1
       rhs(i)= (rhs(i)-ud(i)*tmpv(i)*rhs(i+1))
    end do
  end subroutine lusolvec
      
  !----------------------------------------------------------------
  ! computation of rho which represents the density for
  !  the RHS of the quasi-neutrality equation
  !   rho(r,theta,phi) = (1/n0(r)) * 
  !                      (nGi(r,theta,phi)-nGi_eq(r,theta))
  !  where nGi = \int J0(sqrt{2*mu}).f Jv(r,theta,vpar) dvpar dmu
  !            = \int J0(sqrt{2*mu}).[\int f Jv dvpar] dmu 
  !  with  nGi_eq(r,theta) = \int J0.feq(r,theta,vpar,mu) 
  !                          Jv(r,theta,vpar)dmu dvpar
  !  Jv being the Jacobian in velocity space :
  !   . Jv(r,theta,vpar) = 2*pi*Bstar in the 5D case
  !                      = 1          in the 4D case
  !
  ! Rk : nGi is computed as the sum of 2 integrals:
  !   nGi =  \int dmu (I0*B_norm + 
  !            I1*mi*mu0*vec_J.vec_b/(e*B_norm))
  !   where I0 = J0.(\int f dvpar) and          
  !         I1 = J0.(\int f*vpar dvpar)
  !   this by using the expresion of Bstar
  !    Bstar = B_norm + mi*vpar*mu0*vec_J.vec_b/(e*B_norm)
  !---------------------------------------------------------------
  subroutine compute_rho_FFT2Dhybrid(J0,f,Sfmu_eq, geom, &
    init_prof,init_magnet,init_curr,rho)
    use globals, only : Nbproc_r, Nbproc_loc, Nbproc_mu, &
      dom_r, dom_theta, Nr, Ntheta, Nphi, Nmu, mumin, Zi, &
      istart, iend, jstart, jend, mu_id, mpi_comm_intermu, &
      Rarray1_NrNtheta, Rarray1_NrNthetaDomphi, &
      Rarray_PNrPNthetaNphi_nbM2, Rarray_NrNthetamuphi_nbM2, &
      integration_CS
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp_scoefvpar, &
      Romp1_0Nvpar, Romp2_0Nvpar, &
      Romp1_1Nrp1_1Nthetap1, Romp2_1Nrp1_1Nthetap1, &
      Comp1_1Nrp1_1Ntheta
    use gyroaverage_class
    include "mpiperso.h"
    type(J0operator)             , intent(in) :: J0
    type(fdistribu5d)            , intent(in) :: f
    real(RKIND), dimension(:,:,:), pointer    :: Sfmu_eq
    type(geometry)               , intent(in) :: geom
    type(init_profile)           , intent(in) :: init_prof
    type(init_magnetic)          , intent(in) :: init_magnet
    type(init_current)           , intent(in) :: init_curr
    real(RKIND), dimension(:,:,:), pointer    :: rho
      
    integer, parameter :: nbmoments = 2
    real(RKIND)   , dimension(:)    , pointer :: fmu_1D
    real(RKIND)   , dimension(:)    , pointer :: fmu_v_1D
    real(RKIND)   , dimension(:)    , pointer :: scoefvpar
    real(RKIND)   , dimension(:,:)  , pointer :: J0_intdv_f
    real(RKIND)   , dimension(:,:)  , pointer :: J0_intvdv_f
    real(RKIND)   , dimension(:,:,:), pointer :: M0_loc
    complex(CKIND), dimension(:)    , pointer :: crhs
    complex(CKIND), dimension(:,:)  , pointer :: Acomp
      
    integer     :: ierr
    integer     :: ir, itheta, iphi, ivpar, imu, iphiloc
    logical     :: case4D
    real(RKIND) :: fmu_tmp
    real(RKIND) :: Bnorm_ij, scalprod_mu0Jb_ij
    real(RKIND) :: coeff_intdmu
      
    ! -> variable for parallelization
    integer :: iproc, nbreqr, nbbcast_r, nbbcast_s
    integer :: local_id, dep, base, global_id, startbuf
    integer :: l_istart, l_iend, l_jstart, l_jend
    integer :: itag_phi, nbelements, tid
      
#ifdef TIMER
    integer(TIMEPREC) :: tdeb, tfin, tth1, tth2
    real(RKIND)       :: tbar(12), tcmp(12), tcom(12), tthd
    call clck_time(tdeb)
    tcmp(1:4) = 0.
    tcom(1:4) = 0.
#endif
    
    case4D = .false.
    if ((geom%Nmu.eq.0).and.(mumin.eq.0._RKIND)) &
      case4D = .true.
      
!baoter (ATTENTION redimensionner localement en tenant 
!baoter compte de la taille locale en phi)
    M0_loc => Rarray1_NrNthetaDomphi
!eaoter
    tid = 1
      
    !*** initialization ***
    do iphi = iphistart,iphistart+dom_mapphi-1
      do itheta = 0,geom%Ntheta
        do ir = 0,geom%Nr
          M0_loc(ir,itheta,iphi) = 0._RKIND
        end do
      end do
    end do
      
    !*** computation ***
!$OMP PARALLEL private(tid,iphi,itheta,ir,ivpar,&
!$OMP scoefvpar,fmu_tmp,fmu_1D,fmu_v_1D)  default(shared)
#ifdef _OPENMP
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    scoefvpar => Romp_scoefvpar(tid)%val
    fmu_1D    => Romp1_0Nvpar(tid)%val
    fmu_v_1D  => Romp2_0Nvpar(tid)%val
!$OMP DO SCHEDULE(STATIC)
    do iphi = 0,f%n3-1
      do itheta = f%jstart,f%jend
        do ir = f%istart,f%iend
          do ivpar = 0,f%n4
            fmu_tmp         = f%values(ir,itheta,iphi,ivpar) - &
              Sfmu_eq(ir,itheta,ivpar)
            fmu_1D(ivpar)   = fmu_tmp
            fmu_v_1D(ivpar) = fmu_tmp*geom%vparg(ivpar)
          end do
          !*** computation of the integrals in vpar direction ***
          if (integration_CS) then
            !-> \int f dvpar
            call compute_omp_vpar_integral_CS(f%n4,f%h4,fmu_1D, &
              f%BCvpar_left,f%BCvpar_right, &
              f%nspline1d_vpar(tid),scoefvpar, &
              Rarray_PNrPNthetaNphi_nbM2(0,ir,itheta,iphi))
            !-> \int f vpar dvpar
            call compute_omp_vpar_integral_CS(f%n4,f%h4,fmu_v_1D, &
              f%BCvpar_left,f%BCvpar_right, &
              f%nspline1d_vpar(tid),scoefvpar, &
              Rarray_PNrPNthetaNphi_nbM2(1,ir,itheta,iphi))
          else
            !-> \int f dvpar
            call compute_vpar_integral_colloc(fmu_1D(0:f%n4), &
              geom,Rarray_PNrPNthetaNphi_nbM2(0,ir,itheta,iphi))
            !-> \int f vpar dvpar
            call compute_vpar_integral_colloc(fmu_v_1D(0:f%n4), &
              geom,Rarray_PNrPNthetaNphi_nbM2(1,ir,itheta,iphi))
          end if
        enddo
      enddo
    end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
      
#ifdef TIMER
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcmp(1))
    tdeb = tfin
      
    call ppbarrier()
      
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tbar(1))
    tdeb = tfin
#endif
      
    !*** distribution of \int f dvpar and \int f vpar dvpar ***
    !-> message receive
    if (iphistart .ne. -1) then
      nbreqr  = 0
      iphiloc = 0
      do iphi = 0,f%n3-1
        if (moments_mapphi(iphi).eq.pglobal_id) then
          do imu = 0, Nmu
            do local_id = 0, Nbproc_loc-1
              global_id  = imu*Nbproc_loc + local_id
              startbuf   = dom_theta * (local_id + Nbproc_loc * &
                (imu + (Nmu+1) * iphiloc)) 
              itag_phi   = 4000 + iphi
              nbelements = nbmoments * dom_r * dom_theta
              call MPI_IRECV(&
                Rarray_NrNthetamuphi_nbM2(0,startbuf), &
                nbelements,MPI_REAL8,global_id,itag_phi, &
                MPI_COMM_WORLD,moments_reqr(nbreqr),ierr)
              nbreqr = nbreqr + 1
            enddo
          enddo
          iphiloc = iphiloc + 1
        endif
      enddo
    endif
      
    !-> message sending 
    do iphi = 0,f%n3-1
      itag_phi   = 4000 + iphi
      nbelements = nbmoments*dom_r*dom_theta
      call MPI_SEND( &
        Rarray_PNrPNthetaNphi_nbM2(0,istart,jstart,iphi),&
        nbelements,MPI_REAL8,moments_mapphi(iphi),itag_phi,&
        MPI_COMM_WORLD,ierr)
    end do
    if (iphistart .ne. -1) &
      call MPI_WAITALL(nbreqr,moments_reqr(0), &
      moments_status(1,0),ierr)
      
#ifdef TIMER
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcom(1))
    tdeb = tfin
#endif
   
    if (iphistart .ne. -1) then
!$OMP PARALLEL private(tid,imu,iphi,coeff_intdmu, &
!$OMP J0_intdv_f,J0_intvdv_f, &
!$OMP Bnorm_ij,scalprod_mu0Jb_ij, &
!$OMP nbelements,local_id,base,dep, &
!$OMP l_istart,l_iend,l_jstart,l_jend,startbuf,&
!$OMP itheta,ir,iphiloc,Acomp)  default(shared)
!$OMP BARRIER
#ifdef _OPENMP
      tid = 1+omp_get_thread_num()
#else
      tid = 1
#endif
      J0_intdv_f  => Romp1_1Nrp1_1Nthetap1(tid)%val
      J0_intvdv_f => Romp2_1Nrp1_1Nthetap1(tid)%val
      Acomp       => Comp1_1Nrp1_1Ntheta(tid)%val
!$OMP DO SCHEDULE(STATIC,1)
      do iphi = 0,geom%Nphi-1
        if (moments_mapphi(iphi) .eq. pglobal_id) then
          iphiloc = iphi - iphistart
          do imu = 0,Nmu
            coeff_intdmu = geom%coeff_intdmu(imu)
            do itheta = 1,geom%Ntheta+1
              do ir = 1,geom%Nr+1
                J0_intdv_f(ir,itheta)  = 0._RKIND
                J0_intvdv_f(ir,itheta) = 0._RKIND
              end do
            end do
            do local_id = 0,Nbproc_loc-1
              base     = (local_id/Nbproc_r)
              dep      = mod(local_id,Nbproc_r)              
              l_istart = dep *  dom_r
              l_iend   = l_istart + dom_r - 1
              l_jstart = base * dom_theta
              l_jend   = l_jstart + dom_theta - 1
              startbuf = dom_theta * (local_id + Nbproc_loc * &
                (imu + (Nmu+1) * iphiloc)) 
              !*** J0_intdv_f  = (\int f dvpar) and  ***
              !*** J0_intvdv_f = (\int f vpar dvpar) *** 
              do itheta = l_jstart, l_jend
                do ir = l_istart, l_iend
                  J0_intdv_f(1+ir,1+itheta) = &
                    Rarray_NrNthetamuphi_nbM2(nbmoments * &
                    (ir-l_istart)+0,startbuf+itheta-l_jstart) 
                  J0_intvdv_f(1+ir,1+itheta) = &
                    Rarray_NrNthetamuphi_nbM2(nbmoments * &
                    (ir-l_istart)+1,startbuf+itheta-l_jstart) 
                enddo
              enddo
            enddo
            !*** compute of the gyroaverages ****             
            !-> J0.(\int f dvpar) 
            call omp_compute_gyrofunction_2D(J0,geom,imu, &
              Acomp,J0_intdv_f)
            if (.not.case4D) then
              !-> J0.(\int f vpar dvpar)
              call omp_compute_gyrofunction_2D(J0,geom,imu, &
                Acomp,J0_intvdv_f)
            end if
            !******************************************************
            !*** computation of the integrals in mu             ***
            !***  nGi is computed as the sum of 2 integrals:    ***
            !***   nGi =  \int dmu (I0*B_norm +                 ***
            !***            I1*mi*mu0*vec_J.vec_b/(e*B_norm))   ***
            !***   where I0 = J0.(\int f dvpar) and             ***
            !***         I1 = J0.(\int f*vpar dvpar)            ***
            !***   this by using the expresion of Bstar         ***
            !***    Bstar = B_norm +                            ***
            !***            mi*vpar*mu0*vec_J.vec_b/(e*B_norm)  ***
            !******************************************************
            if (.not.case4D) then
              do itheta = 0,f%n2-1
                do ir = 0,f%n1
                  Bnorm_ij = init_magnet%B_norm(ir,itheta) 
                  call compute_scalprod_mu0Jb(geom,init_magnet, &
                    init_curr,ir,itheta,scalprod_mu0Jb_ij)
                  M0_loc(ir,itheta,iphi) = &
                    M0_loc(ir,itheta,iphi) + &
                    TWOPI*coeff_intdmu * &
                    (J0_intdv_f(1+ir,1+itheta)*Bnorm_ij + &
                    J0_intvdv_f(1+ir,1+itheta)*Zi * &
                    scalprod_mu0Jb_ij/Bnorm_ij)
                end do
              end do
            else
              do itheta = 0,f%n2-1
                do ir = 0,f%n1
                  M0_loc(ir,itheta,iphi) = &
                    M0_loc(ir,itheta,iphi) + &
                    coeff_intdmu*J0_intdv_f(1+ir,1+itheta)
                end do
              end do
            end if
          end do
          !*** periodic conditions in theta ***
          M0_loc(0:f%n1,f%n2,iphi) = M0_loc(0:f%n1,0,iphi)
        endif
      end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
    end if
      
#ifdef TIMER
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcmp(2))
    call ppbarrier()
    tdeb = tfin
#endif
      
    !**** distribution of rho ***
    !-> message receive
    nbbcast_r = 0
    do iphi = 0,f%n3-1
      itag_phi   = 7000 + iphi
      nbelements = (Nr+1)*(Ntheta+1)
      call MPI_IRECV(rho(0,0,iphi),nbelements,MPI_REAL8,&
        moments_mapphi(iphi),itag_phi, &
        MPI_COMM_WORLD,bcast_reqr(nbbcast_r),ierr)
      nbbcast_r = nbbcast_r + 1
    enddo
      
    !-> message sending
    if (iphistart .ne. -1) then
      nbbcast_s = 0
      do iphi = 0,geom%Nphi-1
        if (moments_mapphi(iphi) .eq. pglobal_id) then
          nbelements = (Nr+1)*(Ntheta+1)
          itag_phi   = 7000 + iphi
          do iproc = 0,Nbproc_tot-1
            call MPI_SEND(M0_loc(0,0,iphi),nbelements,MPI_REAL8,&
              iproc,itag_phi,MPI_COMM_WORLD,ierr)
            nbbcast_s = nbbcast_s+1
          enddo
        endif
      end do
    endif
    call MPI_WAITALL(nbbcast_r,bcast_reqr(0), &
      moments_status(1,0),ierr)
      
    
    !*** periodic conditions in phi ***
    do itheta = 0,f%n2
      do ir = 0,f%n1
        rho(ir,itheta,f%n3) = rho(ir,itheta,0)
      end do
    end do
    
    !*** compute                                             *** 
    !***  rho(r,theta,phi) = (1/n0(r))*                      ***
    !***    (\int J0(ni) J(r,theta) dmu - nGi_eq(r,theta))   ***
    !*** where                                               ***
    !***  nGi_eq(r,theta) = \int J0.feq J(r,theta) dmu dvpar ***
    do iphi = 0,geom%Nphi
      do itheta = 0,geom%Ntheta
        do ir = 0,geom%Nr
          rho(ir,itheta,iphi) = rho(ir,itheta,iphi)/init_prof%n0(ir)
        end do
      end do
    end do
      
#ifdef TIMER
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcom(2))
    tdeb = tfin
      
    write(6,'(I4,A,3F13.6,A,4F13.6,A,F13.6)') pglobal_id, &
      " Temps poisson1A, tcomm ",tcom(1:3)," tcomp ",tcmp(1:4), &
      " tbarrier ",tbar(1)
#endif
  end subroutine compute_rho_FFT2Dhybrid
      
  !-------------------------------------------------------
  !  4th order Poisson solver with periodic 
  !  boundary conditions in r and Dirichlet 0 boundary 
  !  conditions in th. phi is considered as a parameter.
  !-------------------------------------------------------
  subroutine solve_poisson_FFT2Dhybrid(pthis,geom, &
    init_prof,init_magnet,init_curr,J0,f,Sfmu_eq)
    use globals, only : mpi_comm_mu, plocal_id, &
      Rarray1_NrNthetaNphi, Carray_1Nrp11Ntheta1Nphi, &
      m, n, m0n0_eq0, single_m, single_n, &
      filter_choice, filter_deltam, &
      zonal_flow, iter_glob, Phi00BCr_Neumann, &
      nbelectrons_diag, nbions_diag
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp1_1Nrp1, Romp2_1Nrp1, &
      Romp3_1Nrp1
    use gyroaverage_class
    include "mpiperso.h"
    type(poisson)      , intent(inout) :: pthis
    type(geometry)     , intent(in)    :: geom
    type(init_profile) , intent(in)    :: init_prof
    type(init_magnetic), intent(in)    :: init_magnet
    type(init_current) , intent(in)    :: init_curr
    type(J0operator)   , intent(in)    :: J0
    type(fdistribu5d)  , intent(in)    :: f
    real(RKIND), &
      dimension(:,:,:) , pointer       :: Sfmu_eq
      
    !*** local variables ***
    integer     :: ierr
    integer     :: i, j, k                   ! loop indices
    integer     :: ir, itheta, iphi          ! loop indices
    integer     :: N1, N2, N3, N2m1, N3m1    ! domain dimensions
    integer     :: stopit, err               ! error flags
    real(RKIND) :: mask, coef_fftinv
    integer     :: tid, isign
      
    real(RKIND)   , dimension(:,:,:), pointer :: rho
    complex(CKIND), dimension(:,:,:), pointer :: FFT_Phi
    real(RKIND)   , dimension(:)    , pointer :: ldiag, udiag, diag
      
    !*** variables for the diagonal filter ***
    real(RKIND) :: ktheta, kphi, ktheta0, ktheta_min, ktheta_max
    real(RKIND) :: qmax, qmin, qloc
      
#ifdef TIMER
    integer(TIMEPREC) :: tdeb, tfin, tth1, tth2
    real(RKIND)       :: tbar(12), tcmp(12), tcom(12), tthd
    call clck_time(tdeb)
    tcmp(1:4) = 0.
    tcom(1:4) = 0.
#endif
      
    !*** initialize local variables ***
    N1   = size(geom%rg)
    N2   = size(geom%thetag)
    N3   = size(geom%phig)
    N2m1 = N2-1
    N3m1 = N3-1
   
    rho     => Rarray1_NrNthetaNphi
    FFT_Phi => Carray_1Nrp11Ntheta1Nphi 
    
    !*** compute the RHS of Poisson equation   ***
    !***  rho = 1/n0(\int(J0(f)-J0(feq)dv)     ***
    call compute_rho_FFT2Dhybrid(J0,f,Sfmu_eq,geom, &
      init_prof,init_magnet,init_curr,rho)
      
#ifdef TIMER
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcmp(1))
    tdeb = tfin
#endif
      
    !*** perform FFT 2D in theta and phi  ***
    isign  = 1
    stopit = 0    
!$OMP PARALLEL private(tid,ldiag,udiag,diag,k,j,i,err,mask, &
!$OMP ktheta,kphi,ktheta0,ktheta_min,ktheta_max,qmin,qmax,qloc) &
!$OMP firstprivate(isign) default(shared)
!$OMP BARRIER
#ifdef _OPENMP
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    ldiag => Romp1_1Nrp1(tid)%val
    udiag => Romp2_1Nrp1(tid)%val
    diag  => Romp3_1Nrp1(tid)%val
      
!$OMP DO SCHEDULE(STATIC)
    do i = 1,N1 
      do k = 1,N3m1
        do j = 1,N2m1
          omp_data_i(tid)%val(j,k) = rho(i-1,j-1,k-1)
        end do
      end do
      call fourrow(omp_data_i(tid)%val,isign)
      do k = 1,N3m1
        do j = 1,N2m1
          omp_tmp(tid)%val(k,j) = omp_data_i(tid)%val(j,k)
        end do
      end do
      call fourrow(omp_tmp(tid)%val,isign)
      do k = 1,N3m1
        do j = 1,N2m1
          FFT_Phi(i,j,k) = omp_tmp(tid)%val(k,j)
        end do
      end do
    enddo
!$OMP END DO
!$OMP BARRIER
    
!$OMP DO SCHEDULE(STATIC)
    do k = 1,N3m1
      !*** solving of tridiagonal systems L Phi^{m,n}=rho^{m,n} ***
      do j = 1,N2m1
        omp_rhs(tid)%val(1:N1-2) = FFT_Phi(2:N1-1,j,k)
        !*** solve the LU system with LAPACK routine ***
        if ((j.eq.1).and.(k.eq.1)) then
          do i = 1,N1-3
            ldiag(i) = pthis%ldiag_m0n0(j,i)
            udiag(i) = pthis%udiag_m0n0(j,i)
          end do
          do i = 1,N1-2
            diag(i) = pthis%diag_m0n0(j,i)
          end do
          call thomas_c(N1-2,udiag,ldiag,diag,omp_rhs(tid)%val,err)
        else
          do i = 1,N1-3
            ldiag(i) = pthis%ldiag(j,i)
            udiag(i) = pthis%udiag(j,i)
          end do
          do i = 1,N1-2
            diag(i)  = pthis%diag(j,i)
          end do
          call thomas_c(N1-2,udiag,ldiag,diag,omp_rhs(tid)%val,err)
        end if
        if (err.ne.0) then
          print *, 'solve_poisson_par: ', &
            'problem in Poisson solving for k = ', k
!$OMP CRITICAL
          stopit = err
!$OMP END CRITICAL
        endif
        !*** copy solution in FFT_Phi ***
        do i = 2,N1-1
          FFT_Phi(i,j,k) = omp_rhs(tid)%val(i-1)
        end do
      
        !*** boundary conditions ***
        if (Phi00BCr_Neumann) then
          if (j.eq.1) then 
            !-> treatment of the mode (0,0)
            !---> Neumann at r=rmin
            FFT_Phi(1,j,k)  = FFT_Phi(2,j,k)
            !---> Dirichlet at r=rmax
            FFT_Phi(N1,j,k) = 0._CKIND
          else
            !-> treatment of the modes (m,n).ne.0
            !---> Dirichlet at r=rmin
            FFT_Phi(1,j,k)  = 0._CKIND
            !---> Dirichlet at r=rmax 
            FFT_Phi(N1,j,k) = 0._CKIND
          end if
        else
          if (j.eq.1) then
            !-> treatment of the mode (0,0)
            if (geom%rg(0).ne.0) then 
              ! Dirichlet at r=rmin if rmin.ne.0
              FFT_Phi(1,j,k) = 0._CKIND
            else
              ! Neumann   at r=rmin if rmin=0
              FFT_Phi(1,j,k) = FFT_Phi(2,j,k)
            end if
            ! Dirichlet at r=rmax
            FFT_Phi(N1,j,k) = 0._CKIND
          else
            !-> treatment of the modes (m,n).ne.0
            ! Dirichlet at r=rmin
            FFT_Phi(1,j,k)  = 0._CKIND
            ! Dirichlet at r=rmin
            FFT_Phi(N1,j,k) = 0._CKIND
          end if
        end if
      end do
    end do
!$OMP END DO
      
!$OMP BARRIER
!$OMP MASTER
    !*** Diagonal filter [ m = n*q(r) +- delta_m ]          ***
    !*** with m0 = n*q(r), then FFT_PHI is put equal to 0   ***
    !***  if abs(m-m0)>delta_m                              ***
    !*** Rk : this test is equivalent to                    ***
    !***  ktheta = kphi*q(r)*(Lphi/Ltheta)                  ***
    !***           +- 2*pi*delta_m/Ltheta                   ***
    !***   with ktheta = 2*pi*m/Ltheta and                  ***
    !***        kphi = 2*pi*n/Lphi                          ***
    if (filter_choice.eq.1) then
      ktheta_min = minval(pthis%kth)
      ktheta_max = maxval(pthis%kth)
      do k = 1,N3m1
        kphi = pthis%kphi(k)
        do j = 1,N2m1
          ktheta = pthis%kth(j)
          do i = 1,N1
            ktheta0 = -kphi*geom%Lphi / &
              (init_prof%iota(i-1)*geom%Ltheta)
            ktheta0 = max(ktheta0,ktheta_min)
            ktheta0 = min(ktheta0,ktheta_max)
            if ( abs(ktheta-ktheta0).gt.filter_deltam ) then
              FFT_Phi(i,j,k) = 0._CKIND
            end if
          end do
        end do
      end do
      
    !*** Non res. modes filter           ***
    else if (filter_choice .eq. 2) then
!AS!      ktheta_min = minval(pthis%kth)
!AS!      ktheta_max = maxval(pthis%kth)
      qmax = maxval(1._RKIND/init_prof%iota)
      qmin = minval(1._RKIND/init_prof%iota)
      do k = 1,N3m1
        kphi = pthis%kphi(k)
        do j = 1,N2m1
          ktheta = pthis%kth(j)
          if (kphi .eq. 0) then
            if (ktheta .ne. 0) FFT_Phi(i,j,k) = 0._CKIND
          else
            qloc = - ktheta/kphi
            if ((qloc .gt. qmax).or.(qloc .lt. qmin)) then
              do i = 1,N1
                FFT_Phi(i,j,k) = 0._CKIND
              end do
            end if
          end if
        end do
      end do
    end if
      
    !*** the single n mode is retained ***
    if (single_n) then
      do k = 1,N3/2
        mask =  ON - min(ON,abs(ON*(k-1-n)))
        do j = 1,N2m1
          do i = 1,N1
            FFT_Phi(i,j,k) = mask*FFT_Phi(i,j,k)
          end do
        end do
      end do
      do k = N3/2+1,N3m1
        mask =  ON - min(ON,abs(ON*(N3m1-k+1-n)))
        do j = 1,N2m1
          do i = 1,N1
            FFT_Phi(i,j,k) = mask*FFT_Phi(i,j,k)
          end do
        end do
      end do
    end if
      
    !*** the single m mode is retained ***
    if (single_m) then
      do j = 1,N2/2
        mask =  ON - min(ON,abs(ON*(j-1-m)))
        do k = 1,N3m1
          do i = 1,N1
            FFT_Phi(i,j,k) = mask*FFT_Phi(i,j,k)
          end do
        end do
      end do
      do j = N2/2+1,N2m1
        mask =  ON - min(ON,abs(ON*(N2m1-j+1-m)))
        do k = 1,N3m1
          do i = 1,N1
            FFT_Phi(i,j,k) = mask*FFT_Phi(i,j,k)
          end do
        end do
      end do
    end if
      
    !*** the mode (0,0) is forced to 0 ***
    if (m0n0_eq0) then
      do i = 1,N1
        FFT_Phi(i,1,1) = 0._CKIND
      end do
    end if
      
    !*** Saving of the mode (0,0) (zonal flows) ***
    do i = 0,geom%Nr
      Phi00_diag(i) = real(FFT_Phi(i+1,1,1))/(N2m1*N3m1)
    end do
!$OMP END MASTER
!$OMP BARRIER
      
    !*** Perform inverse FFT 2D of the system solution ***
    isign       = -1
    coef_fftinv = 1._RKIND/(geom%Ntheta*geom%Nphi)
!$OMP DO SCHEDULE(STATIC)
    do i = 1,N1 
      do k = 1,N3m1
        do j = 1,N2m1
          omp_data_i(tid)%val(j,k) = FFT_Phi(i,j,k)
        end do
      end do
      call fourrow(omp_data_i(tid)%val,isign)
      do k = 1,N3m1
        do j = 1,N2m1
          omp_tmp(tid)%val(k,j) = omp_data_i(tid)%val(j,k)
        end do
      end do
      call fourrow(omp_tmp(tid)%val,isign)
      do k = 0,geom%Nphi-1
        do j = 0,geom%Ntheta-1
          pthis%Phi(i-1,j,k) = &
            dreal(omp_tmp(tid)%val(k+1,j+1))*coef_fftinv
        enddo
      enddo
      
      !*** duplicate periodic value ***
      do itheta = 0,geom%Ntheta-1
        pthis%Phi(i-1,itheta,geom%Nphi) = pthis%Phi(i-1,itheta,0)
      end do
      do iphi = 0,geom%Nphi
        pthis%Phi(i-1,geom%Ntheta,iphi) = pthis%Phi(i-1,0,iphi)
      end do
    end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
    if (stopit.ne.0) then
      print*,'solve_poisson : problem in Poisson solving for k = '
      stop
    end if
      
    !**********************************************
    !*** computation of the number of electrons ***
    !***  and the number of ions                ***
    !**********************************************
    call compute_nbparticles(geom,init_prof,rho, &
      pthis%Phi,Phi00_diag,nbelectrons_diag,nbions_diag)
      
#ifdef TIMER
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcmp(2))
    tdeb = tfin
    write(6,'(I4,A,2F13.6)') pglobal_id, &
         " Temps poisson1B,  tcomp ",tcmp(1:2)
#endif
  end subroutine solve_poisson_FFT2Dhybrid
      
  !----------------------------------------------------------------
  ! computation of rho which represents the density for
  !  the RHS of the quasi-neutrality equation
  !   rho(r,theta,phi) = (1/n0(r)) * 
  !                      (nGi(r,theta,phi)-nGi_eq(r,theta))
  !  where nGi = \int J0(sqrt{2*mu}).f Jv(r,theta,vpar) dvpar dmu
  !            = \int J0(sqrt{2*mu}).[\int f Jv dvpar] dmu 
  !  with  nGi_eq(r,theta) = \int J0.feq(r,theta,vpar,mu) 
  !                          Jv(r,theta,vpar)dmu dvpar
  !  Jv being the Jacobian in velocity space :
  !   . Jv(r,theta,vpar) = 2*pi*Bstar in the 5D case
  !                      = 1          in the 4D case
  !
  ! Rk : nGi is computed as the sum of 2 integrals:
  !   nGi =  \int dmu (I0*B_norm + 
  !            I1*mi*mu0*vec_J.vec_b/(e*B_norm))
  !   where I0 = J0.(\int f dvpar) and          
  !         I1 = J0.(\int f*vpar dvpar)
  !   this by using the expresion of Bstar
  !    Bstar = B_norm + mi*vpar*mu0*vec_J.vec_b/(e*B_norm)
  !----------------------------------------------------------------
  subroutine compute_rho_FFT1Dhybrid(pthis,J0,f,Sfmu_eq,geom, &
    init_prof,init_magnet,init_curr,rho,orho,ophi)
    use globals, only : Nbproc_r, Nbproc_loc, Nbproc_mu, &
      dom_r, dom_theta, Nr, Ntheta, Nphi, Nmu, mumin, Zi, &
      istart, iend, jstart, jend, mu_id, mpi_comm_intermu, &
      Rarray1_NrNtheta, Rarray1_NrNthetaDomphi, &
      Rarray_PNrPNthetaNphi_nbM2, Rarray_NrNthetamuphi_nbM2, &
      Phi00BCr_Neumann, integration_CS
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp_scoefvpar, &
      Romp1_0Nr, Comp1_1Nrp1_1Ntheta, &
      Romp1_1Nrp1_1Nthetap1, Romp2_1Nrp1_1Nthetap1, &
      Romp1_1Nrp1, Romp2_1Nrp1, Romp3_1Nrp1, Romp4_1Nrp1, &
      Romp1_0Nvpar, Romp2_0Nvpar
    use gyroaverage_class
    include "mpiperso.h"
    type(poisson)            , intent(in) :: pthis
    type(J0operator)         , intent(in) :: J0
    type(fdistribu5d)        , intent(in) :: f
    real(RKIND), &
      dimension(:,:,:)       , pointer    :: Sfmu_eq
    type(geometry)           , intent(in) :: geom
    type(init_profile)       , intent(in) :: init_prof
    type(init_magnetic)      , intent(in) :: init_magnet
    type(init_current)       , intent(in) :: init_curr
    real(RKIND), dimension(:), pointer    :: orho
    real(RKIND), dimension(:), pointer    :: ophi
      
    real(RKIND), dimension(:)    , pointer :: ldiag
    real(RKIND), dimension(:)    , pointer :: udiag
    real(RKIND), dimension(:)    , pointer :: diag
    real(RKIND), dimension(:)    , pointer :: tmpv
    real(RKIND), dimension(:,:,:), pointer :: rho
    real(RKIND), dimension(:)    , pointer :: mrho
      
    integer, parameter :: nbmoments = 2
    real(RKIND)   , dimension(:)    , pointer :: fmu_1D
    real(RKIND)   , dimension(:)    , pointer :: fmu_v_1D
    real(RKIND)   , dimension(:)    , pointer :: scoefvpar
    real(RKIND)   , dimension(:,:)  , pointer :: J0_intdv_f
    real(RKIND)   , dimension(:,:)  , pointer :: J0_intvdv_f
    real(RKIND)   , dimension(:,:,:), pointer :: M0_loc
    complex(CKIND), dimension(:)    , pointer :: crhs
    complex(CKIND), dimension(:,:)  , pointer :: Acomp
      
    integer     :: ierr
    integer     :: ir, itheta, iphi, ivpar, imu, iphiloc
    integer     :: i, N1
    logical     :: case4D
    real(RKIND) :: fmu_tmp
    real(RKIND) :: Bnorm_ij, scalprod_mu0Jb_ij
    real(RKIND) :: coeff_intdmu
    real(RKIND) :: vpar, coeff_dthdphi
      
    ! -> variable for parallelization
    integer     :: iproc, nbreqr, nbbcast_r, nbbcast_s
    integer     :: local_id, dep, base, global_id, startbuf
    integer     :: l_istart, l_iend, l_jstart, l_jend
    integer     :: nbelements, itag_phi, tid
      
#ifdef TIMER
    integer(TIMEPREC) :: tdeb, tfin, tth1, tth2
    real(RKIND)       :: tbar(12), tcmp(12), tcom(12), tthd
    call clck_time(tdeb)
    tcmp(1:4) = 0.
    tcom(1:4) = 0.
#endif
    
    case4D = .false.
    if ((geom%Nmu.eq.0).and.(mumin.eq.0._RKIND)) &
      case4D = .true.
      
!baoter (ATTENTION redimensionner localement en tenant 
!baoter compte de la taille locale en phi)
    M0_loc => Rarray1_NrNthetaDomphi
!eaoter
      
!$OMP PARALLEL private(tid,imu,iphi,itheta,ir,ivpar,vpar,&
!$OMP scoefvpar,fmu_tmp,fmu_1D,fmu_v_1D)  default(shared)
#ifdef _OPENMP
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    scoefvpar => Romp_scoefvpar(tid)%val
    fmu_1D    => Romp1_0Nvpar(tid)%val
    fmu_v_1D  => Romp2_0Nvpar(tid)%val
      
    !********************************************
    !***  Computation of the two integrals:   ***
    !***    . \int f dvpar                    ***
    !***    . \int f*vpar dvpar               ***
    !********************************************
!$OMP DO SCHEDULE(STATIC)
    do iphi = 0,f%n3-1
      do itheta = f%jstart,f%jend
        do ir = f%istart,f%iend
          do ivpar = 0,f%n4
            fmu_tmp         = f%values(ir,itheta,iphi,ivpar) - &
              Sfmu_eq(ir,itheta,ivpar)
            fmu_1D(ivpar)   = fmu_tmp
            fmu_v_1D(ivpar) = fmu_tmp*geom%vparg(ivpar)
          end do
          !*** computation of the integrals in vpar direction ***
          if (integration_CS) then
            !-> \int f dvpar
            call compute_omp_vpar_integral_CS(f%n4,f%h4,fmu_1D, &
              f%BCvpar_left,f%BCvpar_right, &
              f%nspline1d_vpar(tid),scoefvpar, &
              Rarray_PNrPNthetaNphi_nbM2(0,ir,itheta,iphi))
            !-> \int f vpar dvpar
            call compute_omp_vpar_integral_CS(f%n4,f%h4,fmu_v_1D, &
              f%BCvpar_left,f%BCvpar_right, &
              f%nspline1d_vpar(tid),scoefvpar, &
              Rarray_PNrPNthetaNphi_nbM2(1,ir,itheta,iphi))
          else
            !-> \int f dvpar
            call compute_vpar_integral_colloc(fmu_1D(0:f%n4), &
              geom,Rarray_PNrPNthetaNphi_nbM2(0,ir,itheta,iphi))
            !-> \int f vpar dvpar
            call compute_vpar_integral_colloc(fmu_v_1D(0:f%n4), &
              geom,Rarray_PNrPNthetaNphi_nbM2(1,ir,itheta,iphi))
          end if
        enddo
      enddo
    end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
      
#ifdef TIMER
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcmp(1))
    tdeb = tfin
#endif
      
    call ppbarrier()
      
#ifdef TIMER
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tbar(1))
    tdeb = tfin
#endif
      
    !*** distribution of \int f dvpar and \int f vpar dvpar ***
    !-> message receive
    if (iphistart .ne. -1) then
      nbreqr  = 0
      iphiloc = 0
      do iphi = 0,f%n3-1
        if (moments_mapphi(iphi).eq.pglobal_id) then
          do imu = 0,Nmu
            do local_id = 0,Nbproc_loc-1
              global_id  = imu*Nbproc_loc + local_id
              startbuf   = dom_theta * &
                (local_id + Nbproc_loc * (imu + (Nmu+1) * iphiloc)) 
              itag_phi   = 4000 + iphi
              nbelements = nbmoments * dom_r * dom_theta
              call MPI_IRECV( &
                Rarray_NrNthetamuphi_nbM2(0,startbuf),nbelements,&
                MPI_REAL8,global_id,itag_phi,MPI_COMM_WORLD, &
                moments_reqr(nbreqr),ierr)
              nbreqr = nbreqr + 1
            enddo
          enddo
          iphiloc = iphiloc + 1
        endif
      enddo
    endif
      
    !-> message sending 
    do iphi = 0,f%n3-1
      itag_phi   = 4000 + iphi
      nbelements = nbmoments*dom_r*dom_theta
      call MPI_SEND( &
        Rarray_PNrPNthetaNphi_nbM2(0,istart,jstart,iphi),&
        nbelements,MPI_REAL8,moments_mapphi(iphi),itag_phi,&
        MPI_COMM_WORLD,ierr)
    end do
    if (iphistart .ne. -1) &
      call MPI_WAITALL(nbreqr,moments_reqr(0), &
      moments_status(1,0),ierr)
      
#ifdef TIMER
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcom(1))
    tdeb = tfin
#endif
      
    if (iphistart .ne. -1) then
!$OMP PARALLEL private(mrho,tid,imu,iphi, &
!$OMP coeff_intdmu,J0_intdv_f,J0_intvdv_f,nbelements, &
!$OMP Bnorm_ij,scalprod_mu0Jb_ij, &
!$OMP local_id,base,dep,l_istart,l_iend,l_jstart,l_jend,startbuf, &
!$OMP itheta,ir,iphiloc,Acomp) default(shared)
!$OMP BARRIER
#ifdef _OPENMP
      tid = 1+omp_get_thread_num()
#else
      tid = 1
#endif
      mrho        => Romp1_0Nr(tid)%val
      J0_intdv_f  => Romp1_1Nrp1_1Nthetap1(tid)%val
      J0_intvdv_f => Romp2_1Nrp1_1Nthetap1(tid)%val
      Acomp       => Comp1_1Nrp1_1Ntheta(tid)%val
      
      !***********************************************
      !***  Computation of the gyroaverage applied ***
      !***   to the two integrals \int f dvpar and ***
      !***   \int f*vpar dvpar, i.e:               ***
      !***      . J0.(\int f dvpar)                ***
      !***      . J0.(\int f*vpar dvpar)           ***
      !***********************************************
!$OMP DO SCHEDULE(STATIC,1)
      do iphi = 0,geom%Nphi-1
        if (moments_mapphi(iphi) .eq. pglobal_id) then
          do itheta = 0,geom%Ntheta
            do ir = 0,geom%Nr
              M0_loc(ir,itheta,iphi) = 0._RKIND
            end do
          end do
          iphiloc = iphi - iphistart
          do imu = 0,Nmu
            coeff_intdmu = geom%coeff_intdmu(imu)
            do itheta = 1,geom%Ntheta+1
              do ir = 1,geom%Nr+1
                J0_intdv_f(ir,itheta)  = 0._RKIND
                J0_intvdv_f(ir,itheta) = 0._RKIND
              end do
            end do
            do local_id = 0,Nbproc_loc-1
              base     = (local_id/Nbproc_r)
              dep      = mod(local_id,Nbproc_r)              
              l_istart = dep *  dom_r
              l_iend   = l_istart + dom_r - 1
              l_jstart = base * dom_theta
              l_jend   = l_jstart + dom_theta - 1
              startbuf = dom_theta * &
                (local_id + Nbproc_loc * (imu + (Nmu+1) * iphiloc)) 
              !*** J0_intdv_f  = (\int f dvpar) and  ***
              !*** J0_intvdv_f = (\int f vpar dvpar) *** 
              do itheta = l_jstart, l_jend
                do ir = l_istart, l_iend
                  J0_intdv_f(1+ir,1+itheta)  = &
                    Rarray_NrNthetamuphi_nbM2( &
                    nbmoments*(ir-l_istart)+0, &
                    startbuf+itheta-l_jstart) 
                  J0_intvdv_f(1+ir,1+itheta) = &
                    Rarray_NrNthetamuphi_nbM2( &
                    nbmoments*(ir-l_istart)+1, &
                    startbuf+itheta-l_jstart) 
                enddo
              enddo
            enddo
            !*** compute of the gyroaverages ****             
            !-> J0.(\int f dvpar)
            call omp_compute_gyrofunction_2D(J0,geom,imu, &
              Acomp,J0_intdv_f)
            if (.not.case4D) then
              !-> J0.(\int f vpar dvpar)
              call omp_compute_gyrofunction_2D(J0,geom,imu, &
                Acomp,J0_intvdv_f)
            end if
            !******************************************************
            !*** computation of the integrals in mu             ***
            !***  nGi is computed as the sum of 2 integrals:    ***
            !***   nGi =  \int dmu (I0*B_norm +                 ***
            !***            I1*mi*mu0*vec_J.vec_b/(e*B_norm))   ***
            !***   where I0 = J0.(\int f dvpar) and             ***
            !***         I1 = J0.(\int f*vpar dvpar)            ***
            !***   this by using the expresion of Bstar         ***
            !***    Bstar = B_norm +                            ***
            !***            mi*vpar*mu0*vec_J.vec_b/(e*B_norm)  ***
            !******************************************************
            if (.not.case4D) then
              do itheta = 0,f%n2-1
                do ir = 0,f%n1
                  Bnorm_ij = init_magnet%B_norm(ir,itheta) 
                  call compute_scalprod_mu0Jb(geom,init_magnet, &
                    init_curr,ir,itheta,scalprod_mu0Jb_ij)
                  M0_loc(ir,itheta,iphi) = &
                    M0_loc(ir,itheta,iphi)+TWOPI*coeff_intdmu * &
                    (J0_intdv_f(1+ir,1+itheta)*Bnorm_ij + &
                    J0_intvdv_f(1+ir,1+itheta)*Zi * &
                    scalprod_mu0Jb_ij/Bnorm_ij)
                end do
              end do
            else
              do itheta = 0,f%n2-1
                do ir = 0,f%n1
                  M0_loc(ir,itheta,iphi) = &
                    M0_loc(ir,itheta,iphi) + coeff_intdmu * &
                    J0_intdv_f(1+ir,1+itheta)
                end do
              end do
            end if
          end do
          
          !*** compute rho(r,theta,phi) = (1/n0(r))*            ***
          !***  (\int J0(ni) J(r,theta) dmu-nGi_eq(r,theta))    ***
          !*** where                                            ***
          !*** nGi_eq(r,theta)=\int J0.feq J(r,theta) dmu dvpar ***
          !***  AND compute the averaged value of rho           ***
          !***  over theta in mrho                              *** 
          mrho(0:f%n1) = 0._RKIND
          do itheta = 0,f%n2-1
            do ir = 0,f%n1
              M0_loc(ir,itheta,iphi) = M0_loc(ir,itheta,iphi) / &
                init_prof%n0(ir)          
              mrho(ir) = mrho(ir) + M0_loc(ir,itheta,iphi)
            end do
          end do
          !*** Store the average of rho over theta into   ***
          !***  the boundaries of M0_loc in order to send ***
          !***  it to other processors                    ***
          M0_loc(0:f%n1,f%n2,iphi) = mrho(0:f%n1)
        endif
      end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
    endif
      
#ifdef TIMER
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcmp(2))
    call ppbarrier()
    tdeb = tfin
#endif
      
    !*** distibution of rho ***
    !-> message receive
    nbbcast_r = 0
    do iphi = 0,f%n3-1
      itag_phi   = 7000 + iphi
      nbelements = (Nr+1)*(Ntheta+1)
      call MPI_IRECV(rho(0,0,iphi),nbelements,MPI_REAL8, &
        moments_mapphi(iphi),itag_phi,MPI_COMM_WORLD, &
        bcast_reqr(nbbcast_r),ierr)
      nbbcast_r = nbbcast_r + 1
    enddo
      
    !-> message sending
    if (iphistart .ne. -1) then
      nbbcast_s = 0
      do iphi = 0,geom%Nphi-1
        if (moments_mapphi(iphi) .eq. pglobal_id) then
          nbelements = (Nr+1)*(Ntheta+1)
          itag_phi   = 7000 + iphi
          do iproc = 0,Nbproc_tot-1
            call MPI_SEND(M0_loc(0,0,iphi),nbelements,MPI_REAL8, &
              iproc,itag_phi,MPI_COMM_WORLD,ierr)
            nbbcast_s = nbbcast_s+1
          enddo
        endif
      end do
    endif
    call MPI_WAITALL(nbbcast_r,bcast_reqr(0), &
      moments_status(1,0),ierr)
      
#ifdef TIMER
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcom(2))
    tdeb = tfin
#endif
      
    orho(0:f%n1) = 0._RKIND
    do iphi = 0,geom%Nphi-1
      !*** accumulate the averaged value of rho ***
      !***  over theta (packed into rho)        ***
      orho(0:f%n1) = orho(0:f%n1) + rho(0:f%n1,f%n2,iphi)
      !*** restore periodic conditions in theta ***
      rho(0:f%n1,f%n2,iphi) = rho(0:f%n1,0,iphi)
    enddo
    coeff_dthdphi   = 1._RKIND/(geom%Nphi*geom%Ntheta)
    orho(0:geom%Nr) = orho(0:geom%Nr)*coeff_dthdphi
      
    !*** periodic conditions in phi ***
    do itheta = 0,f%n2
      do ir = 0,f%n1
        rho(ir,itheta,f%n3) = rho(ir,itheta,0)
      end do
    end do
    
    ! === phi solve
    tid = 1
    ldiag => Romp1_1Nrp1(tid)%val
    udiag => Romp2_1Nrp1(tid)%val
    diag  => Romp3_1Nrp1(tid)%val
    tmpv  => Romp4_1Nrp1(tid)%val
      
    N1             = f%n1+1
    itheta         = 1
    do i = 1,N1-3
      ldiag(i)  = pthis%ldiag_m0n0(itheta,i) 
      udiag(i)  = pthis%udiag_m0n0(itheta,i)
    end do
    do i = 1,N1-2
      diag(i)   = pthis%diag_m0n0(itheta,i) 
      ophi(i)   = orho(i)
    end do
    call lusolver(N1-2,udiag,ldiag,diag,ophi,tmpv,ierr) 
    if (Phi00BCr_Neumann) then
      ophi(0)    = ophi(1)
      ophi(N1-1) = 0._RKIND
    else
      if (geom%rg(0).eq.0_RKIND) then 
        ophi(0) = ophi(1)
      else
        ophi(0) = 0._RKIND
      end if
      ophi(N1-1) = 0._RKIND
    end if
      
#ifdef TIMER
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcmp(3))
    tdeb = tfin
      
    write(6,'(I4,A,3F13.6,A,4F13.6,A,F13.6)') pglobal_id, &
      " Temps poisson2A, tcomm ", tcom(1:3)," tcomp ",tcmp(1:4), &
      " tbarrier ",tbar(1)
#endif
  end subroutine compute_rho_FFT1Dhybrid
      
  !-------------------------------------------------------
  ! Choice of the QN solver
  !-------------------------------------------------------
  subroutine solve_poisson(pthis,geom,init_prof, &
    init_magnet,init_curr,J0,f,Sfmu_eq)
    use globals, only : &
      solve_poiss
    use gyroaverage_class
    use clock_module
      
    type(poisson)      , intent(inout) :: pthis
    type(geometry)     , intent(in)    :: geom
    type(init_profile) , intent(in)    :: init_prof
    type(init_magnetic), intent(in)    :: init_magnet
    type(init_current) , intent(in)    :: init_curr
    type(J0operator)   , intent(inout) :: J0
    type(fdistribu5d)  , intent(inout) :: f
    real(RKIND), &
      dimension(:,:,:) , pointer       :: Sfmu_eq
      
!R3 #include "r3_info.h" !R3
      
!R3 call r3_info_begin (r3_info_index_0, 'CPU_poisson')      !R3
      
#ifdef TIMER
    call ppbarrier()
#endif 
    call clck_time(bclock_poisson)
    if (solve_poiss) then
          call solve_poisson_FFT2Dhybrid(pthis,geom, &
            init_prof,init_magnet,init_curr,J0,f,Sfmu_eq)
    else
      pthis%Phi     = 0._RKIND
    end if
      
    call clck_time(eclock_poisson)  
    call clck_diff(bclock_poisson,eclock_poisson, &
      global_time_poisson)
!R3 call r3_info_end (r3_info_index_0) !R3
  end subroutine solve_poisson
      
  !---------------------------------------------------------- 
  ! Computes the number of electrons and ions as:
  !  . nbelectrons = nbions_eq + 
  !                  \int ne(r,theta,phi) Js dr dtheta dphi
  !     with 
  !       --> nbions_eq = \int J0.feq dtau
  !        where dtau = Js Jv dr dtheta d3v
  !
  !       --> ne(r,theta,phi) =  n0(r)/(Te*Zi) * 
  !                             (Phi-lambda*<Phi>_{FS})
  !     where: . lambda = 1 if zonal flows 
  !                    = 0 otherwise
  !            . <Phi>_{FS} = \int Phi Js dtheta dphi /
  !                           \int Js(r,theta) dtheta dphi
  !
  !  . nions = \int J0.feq dtau computed as
  !          = \int J0.(f-feq) dtau + \int J0.feq dtau
  !          = \int QN_RHS*n0 Js dr dtheta dphi + nbions_eq
  !---------------------------------------------------------- 
  subroutine compute_nbparticles(geom,init_prof,QN_RHS, &
    Phi,Phi_avg,Snbelectrons,Snbions)
    use OMPutils_module
    use geometry_class
    use init_profile_class
    use globals           , only : Zi, nbions_eq, &
      zonal_flow, Rarray1_Nr, Nbthread, &
      pglobal_id, outputproc, Rarray2_NrNthetaNphi, diag_targ
    use coord_system_class, only : jacobian_space
    include "mpiperso.h"    
    type(geometry)                  , intent(in)  :: geom
    type(init_profile)              , intent(in)  :: init_prof
    real(RKIND), dimension(:,:,:)   , pointer     :: QN_RHS
    real(RKIND), dimension(0:,0:,0:), intent(in)  :: Phi
    real(RKIND), &
             dimension(0:), optional, intent(in)  :: Phi_avg
    real(RKIND)                     , intent(out) :: Snbelectrons
    real(RKIND)                     , intent(out) :: Snbions
      
    !--> local variables 
    integer     :: status, ierr, request_tmp
    integer, dimension(MPI_STATUS_SIZE) :: status_tmp
    integer     :: ir, itheta, iphi, tid
    real(RKIND) :: jacob_space_tmp
    real(RKIND) :: coeff_intdr_tmp, coeff_intdtheta_tmp
    real(RKIND) :: coeff_intdphi_tmp, intdspace
    real(RKIND), dimension(:), pointer :: Phi_avg_tmp
    real(RKIND) :: nbspecies_tmp(2), nbspecies_loc(2)
    real(RKIND) :: nbions_loc(1:Nbthread), nbelec_loc(1:Nbthread)
      
    !*** initialisation of <Phi>_avg ***
    Phi_avg_tmp => Rarray1_Nr
    if ( present(Phi_avg) .and. zonal_flow ) then
      do ir = 0,geom%Nr
        Phi_avg_tmp(ir) = Phi_avg(ir)
      end do
    else
      do ir = 0,geom%Nr
        Phi_avg_tmp(ir) = 0._RKIND
      end do
    end if
      
    !*** computation of the number of electrons ***
    !***  and the number of ions                ***
#ifdef _OPENMP
!$OMP PARALLEL private(tid,ir,itheta,iphi, &
!$OMP coeff_intdr_tmp,coeff_intdtheta_tmp, &
!$OMP coeff_intdphi_tmp,jacob_space_tmp,intdspace) &
!$OMP default(shared) 
!$OMP BARRIER
      tid = 1+omp_get_thread_num()
#else
      tid = 1
#endif
    nbelec_loc(tid) = 0._RKIND
    nbions_loc(tid) = 0._RKIND
      
    if (iphistart .ne. -1) then
!$OMP DO SCHEDULE(STATIC)
       do iphi = 0,geom%Nphi-1
          if (moments_mapphi(iphi) .eq. pglobal_id) then
             coeff_intdphi_tmp = geom%coeff_intdphi(iphi)
             do itheta = 0,geom%Ntheta-1
                coeff_intdtheta_tmp = geom%coeff_intdtheta(itheta)
                do ir = 0,geom%Nr
                   coeff_intdr_tmp = geom%coeff_intdr(ir)
                   jacob_space_tmp = jacobian_space(ir,itheta)
                   intdspace       = jacob_space_tmp * &
                        coeff_intdr_tmp*coeff_intdtheta_tmp * &
                        coeff_intdphi_tmp
                   !--> number of electrons
                   nbelec_loc(tid) = nbelec_loc(tid) + &
                        init_prof%n0(ir)/(init_prof%Te(ir)*real(Zi)) * &
                        (Phi(ir,itheta,iphi)-Phi_avg_tmp(ir)) * &
                        intdspace
                   !--> number of ions ( using the fact that
                   !-->  rhoi = 1/n0*\int (J0f-J0.feq) Jv d3v )
                   nbions_loc(tid)  = nbions_loc(tid) + &
                        QN_RHS(ir,itheta,iphi)*init_prof%n0(ir) * &
                        intdspace
                end do
             end do
          end if
       end do
!$OMP END DO
    end if
!$OMP BARRIER
!$OMP END PARALLEL
    nbspecies_loc(1) = SUM(nbelec_loc(1:Nbthread))
    nbspecies_loc(2) = SUM(nbions_loc(1:Nbthread))
    nbspecies_tmp(1) = 0._RKIND
    nbspecies_tmp(2) = 0._RKIND
    call MPI_REDUCE(nbspecies_loc,nbspecies_tmp,2, &
         MPI_REAL8,MPI_SUM,diag_targ(1),MPI_COMM_WORLD,ierr)
    nbspecies_loc(1:2) = nbspecies_tmp(1:2)
      
    ! processor with id 'diag_targ(1)' send the computed
    ! densities to processor with id 'outputproc'.
    if (pglobal_id .eq. outputproc) &
         call MPI_IRECV(nbspecies_tmp,2,MPI_REAL8,diag_targ(1), &
         1000,MPI_COMM_WORLD,request_tmp,ierr)
    if (pglobal_id .eq. diag_targ(1)) &
         call  MPI_SEND(nbspecies_loc,2,MPI_REAL8,outputproc, &
         1000,MPI_COMM_WORLD,ierr)
    if (pglobal_id .eq. outputproc) &
         call MPI_wait(request_tmp,status_tmp, ierr) 
   
    Snbelectrons = nbions_eq + nbspecies_tmp(1)
    Snbions      = nbions_eq + nbspecies_tmp(2)
    if (Snbelectrons.ne.Snbelectrons .or. Snbions.ne.Snbions) then
       if (pglobal_id.eq.outputproc) then
          print *,' NAN detected in compute_nbparticles', &
            'so STOP the run !'
          status = 1
          call ppexit()
          call exit(status)
       end if
    end if
  end subroutine compute_nbparticles
end module poisson_class
