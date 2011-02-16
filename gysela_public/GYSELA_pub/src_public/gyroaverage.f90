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
      
!-------------------------------------------------------
! file : gyroaverage.f90
! date : 07/30/2003
! - compute the gyroaverage of the distribution functio
!  and the electric potential
!-------------------------------------------------------
module gyroaverage_class
  use prec_const
  use globals
  use mem_alloc_module
  use OMPutils_module, only : omp2dvector
  use clock_module
  use geometry_class
  use fft_module
      
  implicit none
 
  private
  public :: new_J0operator, del_J0operator
  public :: omp_compute_gyrofunction_2D, compute_gyrofunction_2D, &
    compute_gyrofunction_3D
      
  type, public :: J0operator
     !-> matrix system (1:N1,1:N2)
     integer                                      :: N1, N2      
     real(RKIND)      , dimension(:)    , pointer :: kth
     !-> vector used for the LU factorisation 
     !     of the tridiagonal matricial system
     !--> diagonal terms of matrix
     real(RKIND)      , dimension(:,:,:), pointer :: diag
     !--> lower diagonal terms of matrix
     real(RKIND)      , dimension(:,:,:), pointer :: ldiag
     !--> upper diagonal terms of matrix
     real(RKIND)      , dimension(:,:,:), pointer :: udiag       
     !--> used for OpenMP parallelization
     type(omp2dvector), dimension (:)   , pointer :: omp_diag
     type(omp2dvector), dimension (:)   , pointer :: omp_ldiag
     type(omp2dvector), dimension (:)   , pointer :: omp_udiag
  end type J0operator
      
  !******************************
  contains
  !******************************
      
  !----------------------------------------------------- 
  ! constructor for the J0operator type
  !-----------------------------------------------------   
  subroutine new_J0operator(jthis,geom)
    use globals, only : memory_test
    type(J0operator), intent(out) :: jthis
    type(geometry)  , intent(in)  :: geom 
      
    !*** local variables ***
    integer :: i, j              ! loop index
    integer :: N2m1
   
    !*** system dimension initialization ***
    jthis%N1 = geom%Nr-1             ! size(geom%rg)-2
    jthis%N2 = geom%Ntheta
    N2m1     = jthis%N2-1
      
    !*** memory allocation ***
    call J0_array_allocation(jthis,geom)
      
    if (.not.memory_test) then
      !*** compute kth for the FFT in theta direction ***
      call NR_fft_vector(N2m1,geom%dtheta,jthis%kth)
      
      !*** initialization of the matrix system ***
      call J0_matrix_construction(jthis,geom)
    end if
  end subroutine new_J0operator
      
  !----------------------------------------------------- 
  ! Memory allocation (used in the constructor)
  !-----------------------------------------------------   
  subroutine J0_array_allocation(jthis,geom)
    type(J0operator), intent(inout) :: jthis
    type(geometry)  , intent(in)    :: geom  
      
    integer :: tid
      
    call glob_allocate(jthis%kth,1,jthis%N2,'jthis%kth')
      
    ! -> arrays for LU factorisation
    call glob_allocate(jthis%diag,1,jthis%N1,1,jthis%N2, &
      0,geom%Nmu,'jthis%diag')
    call glob_allocate(jthis%ldiag,1,jthis%N1-1,1,jthis%N2, &
      0,geom%Nmu,'jthis%ldiag')
    call glob_allocate(jthis%udiag,1,jthis%N1-1,1,jthis%N2, &
      0,geom%Nmu,'jthis%udiag')
    allocate(jthis%omp_diag(1:Nbthread))
    allocate(jthis%omp_ldiag(1:Nbthread))
    allocate(jthis%omp_udiag(1:Nbthread))
    do tid = 1,Nbthread
      call glob_allocate(jthis%omp_diag(tid)%val, &
        1,jthis%N1,1,jthis%N2,'jthis%omp_diag(tid)%val')
      call glob_allocate(jthis%omp_ldiag(tid)%val, &
        1,jthis%N1-1,1,jthis%N2,'jthis%omp_ldiag(tid)%val')
      call glob_allocate(jthis%omp_udiag(tid)%val, &
        1,jthis%N1-1,1,jthis%N2,'jthis%omp_udiag(tid)%val')
    enddo
  end subroutine J0_array_allocation
      
  !----------------------------------------------
  ! Factorisation of the matrix associated to the 
  !  system J0(sqrt(2*mu)).Abar = A
  !----------------------------------------------   
  subroutine J0_matrix_construction(jthis,geom)
    type(J0operator), intent(inout) :: jthis
    type(geometry)  , intent(in)    :: geom
      
    integer     :: i, j, imu
    integer     :: err                  ! error flag
    real(RKIND) :: dr2inv, drinv_half, mu_half
    real(RKIND) :: rip1, rip1inv, rip2inv
    
    !*** assemble the matricial system ***
    dr2inv      = 1._RKIND/(geom%dr*geom%dr)
    drinv_half  = 0.5_RKIND/geom%dr
    do imu = 0,geom%Nmu
      mu_half = 0.5_RKIND*geom%mug(imu)
      do j = 1,jthis%N2
        do i = 1,jthis%N1-1
          rip1                 = geom%rg(i) ! rg begins at i = 0
          rip1inv              = 1._RKIND/rip1
          rip2inv              = 1._RKIND/geom%rg(i+1)
          jthis%udiag(i,j,imu) = -mu_half * &
            (dr2inv+drinv_half*rip1inv)
          jthis%ldiag(i,j,imu) = -mu_half * &
            (dr2inv-drinv_half*rip2inv)
          jthis%diag(i,j,imu)  = 1._RKIND + &
            mu_half * (2._RKIND*dr2inv &
            + (jthis%kth(j)*jthis%kth(j))/(rip1*rip1))
        end do
        !*** -> Newmann coditions for boundary conditions ***
        jthis%diag(1,j,imu)        = jthis%diag(1,j,imu) - &
          mu_half*(dr2inv-drinv_half/geom%rg(1))
        rip1                       = geom%rg(jthis%N1)
        jthis%diag(jthis%N1,j,imu) = 1._RKIND + mu_half * &
          (2._RKIND*dr2inv + &
          (jthis%kth(j)*jthis%kth(j))/(rip1*rip1)) - &
          mu_half*(dr2inv+drinv_half/rip1)
      end do
    end do
  end subroutine J0_matrix_construction
      
  !----------------------------------------------------- 
  ! destructor for the poisson type
  !-----------------------------------------------------   
  subroutine del_J0operator(jthis)
    type(J0operator), intent(inout) :: jthis
      
    integer :: tid
      
    !*** case with gyroaverage ***
    call glob_deallocate(jthis%kth)
    call glob_deallocate(jthis%diag)
    call glob_deallocate(jthis%ldiag)
    call glob_deallocate(jthis%udiag)
      
    do tid = 1,Nbthread
      call glob_deallocate(jthis%omp_diag(tid)%val)
      call glob_deallocate(jthis%omp_ldiag(tid)%val)
      call glob_deallocate(jthis%omp_udiag(tid)%val)
    enddo
    deallocate(jthis%omp_diag)
    deallocate(jthis%omp_ldiag)
    deallocate(jthis%omp_udiag)
  end subroutine del_J0operator
      
  !----------------------------------------------------------
  ! Thomas algorithm used for the LU inversion of the
  !  matricial system for the Padde approximation computation 
  !----------------------------------------------------------
  subroutine thomas_algorithm(N1,N2,ldiag,diag,udiag,rhs) 
    integer                       , intent(in) :: N1, N2
    real(RKIND)   , dimension(:,:), pointer    :: ldiag
    real(RKIND)   , dimension(:,:), pointer    :: diag
    real(RKIND)   , dimension(:,:), pointer    :: udiag
    complex(CKIND), dimension(:,:), pointer    :: rhs
      
    integer                      :: i1, i2
    real(RKIND), dimension(1:N1) :: alpha
    real(RKIND)                  :: betai
      
    !*** LU solving ***
    do i2 = 1,N2
      alpha    = 0._RKIND
      alpha(1) = diag(1,i2)
      do i1 = 2,N1
        betai      = ldiag(i1-1,i2)/alpha(i1-1)
        alpha(i1)  = diag(i1,i2) - betai*udiag(i1-1,i2)
        rhs(i1,i2) = rhs(i1,i2) - betai*rhs(i1-1,i2)
      end do
      rhs(N1,i2) = rhs(N1,i2)/alpha(N1)
      do i1 = N1-1,1,-1
        rhs(i1,i2) = (rhs(i1,i2)-udiag(i1,i2) * &
          rhs(i1+1,i2))/alpha(i1)
      end do
    end do
  end subroutine thomas_algorithm
      
  !----------------------------------------------------------
  !  Computation of the gyroaverage for all function
  !   g(r,\theta):
  !   gyro_g(r,\theta)=J0(sqrt(2*mu)).g(r,\theta)
  !----------------------------------------------------------
  subroutine compute_gyrofunction_2D(jthis,geom,imu,A,Abar)
    use OMPutils_module, only : Comp1_1Nrp1_1Ntheta
    type(J0operator)             , intent(in) :: jthis
    type(geometry)               , intent(in) :: geom
    integer                      , intent(in) :: imu
    real(RKIND), dimension(0:,0:), intent(in)  :: A
    real(RKIND), dimension(0:,0:), intent(out) :: Abar
      
    !*** variables locales ***
    integer :: i, j            ! loop indices
    integer :: nr, nth, nthm1  ! domain dimensions
    integer :: midth           ! center of theta dimension
    integer :: rem
      
    complex(CKIND), dimension(:,:), pointer :: A_c
    real(RKIND)   , dimension(:,:), pointer :: ldiag
    real(RKIND)   , dimension(:,:), pointer :: diag
    real(RKIND)   , dimension(:,:), pointer :: udiag
      
    !*** initialize local variables ***
    nr     = size(geom%rg)-1
    nth    = size(geom%thetag)-1
    nthm1  = nth-1
    midth  = nth / 2
      
    !*** in the case of mu=0 J0(mu)=1 ***
    if ((.not.gyroaverage).or.(geom%mug(imu).eq.0._RKIND)) then
      do j = 0,nth
        do i = 0,nr
          Abar(i,j) = A(i,j)
        end do
      end do
    else
       A_c => Comp1_1Nrp1_1Ntheta(1)%val
      
       !*** Perform FFT 1D in theta direction of ***
       !***  the system solution                 ***
       rem = modulo(log(real(nth)),real(log(2._RKIND)))
       if (rem.ne.0._RKIND) then
          write(6,*) 'Warning in fft1_2D : n2 = ', &
            nth, ' is not a power of 2'
          stop
       end if
      
       do j = 1,nth
         do i = 1,nr-1
           A_c(i,j) = A(i,j-1)
         end do
       end do
      
       call fourrow(A_c,1)
      
       !*** algorithm de Thomas for LU matrix system ***
       ldiag => jthis%ldiag(:,:,imu)
       diag  => jthis%diag(:,:,imu)
       udiag => jthis%udiag(:,:,imu)
       ! -> solving of the tridiagonal system only 
       !     on half the Fourier modes
       call thomas_algorithm(nr-1,midth+1,ldiag,diag,udiag,A_c)
       ! -> copy of the results fot the Fourier modes '-j'
       do j = 1,midth+1
          if (j.ne.1 .and. j.ne.midth+1) then 
             !*** copy the result for the fourier mode '-j' ***
             A_c(1:nr-1,nth+2-j) = conjg(A_c(1:nr-1,j))
          endif
       end do
      
       !*** Perform FFT 1D inverse ***
       call fourrow(A_c,-1)
      
       do j = 1,nth
         do i = 1,nr-1
           Abar(i,j-1) = dreal(A_c(i,j))/geom%Ntheta
         end do
       end do
      
       !***  Neumann boundary conditions ***
       Abar(0,:)  = Abar(1,:)
       Abar(nr,:) = Abar(nr-1,:)
      
       !*** duplicate periodic value ***
       Abar(:,nth) = Abar(:,0)
    end if
  end subroutine compute_gyrofunction_2D
      
  !----------------------------------------------------------
  !  Computation of the gyroaverage for all function
  !   g(r,\theta):
  !   gyro_g(r,\theta)=J0(sqrt(2*mu)).g(r,\theta)
  ! (used in a OpenMP parallelized loop)
  ! RK : - Acomp is the temporary array
  !      - Amat received g in input and send J0.g in output
  !----------------------------------------------------------
  subroutine omp_compute_gyrofunction_2D(jthis,geom,imu,Acomp,Amat)
    type(J0operator)              , intent(in) :: jthis
    type(geometry)                , intent(in) :: geom
    integer                       , intent(in) :: imu
    complex(CKIND), dimension(:,:), pointer    :: Acomp
    real(RKIND)   , dimension(:,:), pointer    :: Amat
      
    !*** variables locales ***
    integer :: i, j              ! loop indices
    integer :: nrp1, nth, nthp1  ! domain dimensions
    integer :: midth             ! center of theta dimension
    integer :: rem
      
    real(RKIND), dimension(:,:), pointer :: ldiag
    real(RKIND), dimension(:,:), pointer :: diag
    real(RKIND), dimension(:,:), pointer :: udiag
      
    !*** in the case of mu=0 J0(mu)=1 ***
    if ((gyroaverage).and.(geom%mug(imu).ne.0._RKIND)) then
       !*** initialize local variables ***
       nrp1   = size(geom%rg)
       nthp1  = size(geom%thetag)
       nth    = nthp1-1
       midth  = nth / 2
      
       !*** Perform FFT 1D in theta direction of ***
       !***   the system solution                ***
       rem = modulo(log(real(nth)),real(log(2._RKIND)))
       if (rem.ne.0._RKIND) then
          write(6,*) 'Warning in fft1_2D : n2 = ', &
            nth, ' is not a power of 2'
          stop
       end if
      
       do j = 1,nth
         do i = 1,nrp1-2
           Acomp(i,j) = Amat(i+1,j)
         end do
       end do
      
       call fourrow(Acomp,1)
      
       !*** algorithm de Thomas for LU matrix system ***
       ldiag => jthis%ldiag(:,:,imu)
       diag  => jthis%diag(:,:,imu)
       udiag => jthis%udiag(:,:,imu)
       ! -> solving of the tridiagonal system only on 
       !     half the Fourier modes
       call thomas_algorithm(nrp1-2,midth+1,ldiag,diag,udiag,Acomp)
       ! -> copy of the results fot the Fourier modes '-j'
       do j = 1,midth+1
          if (j.ne.1 .and. j.ne.midth+1) then 
             !*** copy the result for the fourier mode '-j' ***
             Acomp(:,nth+2-j) = conjg(Acomp(:,j))
          endif
       end do
      
       !*** Perform FFT 1D inverse ***
       call fourrow(Acomp,-1)
      
       do j = 1,nth
         do i = 1,nrp1-2
           Amat(i+1,j) = dreal(Acomp(i,j))/geom%Ntheta
         end do
       end do
      
       !***  Neumann boundary conditions ***
       Amat(1,:)    = Amat(2,:)
       Amat(nrp1,:) = Amat(nrp1-1,:)
         
       !*** duplicate periodic value ***
       Amat(:,nthp1) = Amat(:,1)
    end if
  end subroutine omp_compute_gyrofunction_2D
      
  !----------------------------------------------------------
  !  Computation of the gyroaverage for all function
  !   g(r,\theta,phi):
  !   gyro_g(r,\theta,phi)=J0(sqrt(2*mu)).g(r,\theta,phi)
  !----------------------------------------------------------
  subroutine compute_gyrofunction_3D(jthis,geom,imu,A,Abar)
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp1_1Nrp1_1Nthetap1, &
      Comp1_1Nrp1_1Ntheta
    type(J0operator)                , intent(in)  :: jthis
    type(geometry)                  , intent(in)  :: geom
    integer                         , intent(in)  :: imu
    real(RKIND), dimension(0:,0:,0:), intent(in)  :: A     
    real(RKIND), dimension(0:,0:,0:), intent(out) :: Abar 
      
    !*** variables locales ***
    integer :: i, j, iphi      ! loop indices
    integer :: nr, nth, nthm1  ! domain dimensions
    integer :: midth           ! center of theta dimension
    integer :: rem
      
    real(RKIND)   , dimension(:,:), pointer :: ldiag
    real(RKIND)   , dimension(:,:), pointer :: diag
    real(RKIND)   , dimension(:,:), pointer :: udiag
    real(RKIND)   , dimension(:,:), pointer :: Amat
    complex(CKIND), dimension(:,:), pointer :: Acomp
    integer                                 :: tid
!R3 #include "r3_info.h" !R3
      
!R3 call r3_info_begin (r3_info_index_0, 'CPU_gyroaverage')      !R3
    call clck_time(bclock_compJ0)
      
    !*** initialize local variables ***
    nr     = size(geom%rg)-1
    nth    = size(geom%thetag)-1
    nthm1  = nth-1
    midth  = nth / 2
      
    !*** in the case of mu=0 J0(mu)=1 ***
    if ((.not.gyroaverage).or.(geom%mug(imu).eq.0._RKIND)) then
       do iphi = 0,geom%Nphi
        do j = 0,nth
          do i = 0,nr
            Abar(i,j,iphi) = A(i,j,iphi)
          end do
        end do
      end do
    else
      !*** Perform FFT 1D in theta direction of the ***
      !***   system solution                        ***
      rem = modulo(log(real(nth)),real(log(2._RKIND)))
      if (rem.ne.0._RKIND) then
        write(6,*) 'Warning in fft1_2D : n2 = ', &
          nth, ' is not a power of 2'
        stop
      end if
!$OMP PARALLEL private(tid,Amat,Acomp,iphi,ldiag,diag,udiag,i,j), &
!$OMP default(shared)
!$OMP BARRIER
#ifdef _OPENMP
      tid = 1+omp_get_thread_num()
#else
      tid = 1
#endif
      Amat  => Romp1_1Nrp1_1Nthetap1(tid)%val
      Acomp => Comp1_1Nrp1_1Ntheta(tid)%val
      ldiag => jthis%omp_ldiag(tid)%val
      diag  => jthis%omp_diag(tid)%val
      udiag => jthis%omp_udiag(tid)%val
      do j = 1,jthis%N2 
         do i = 1,jthis%N1 
            diag(i,j)  = jthis%diag(i,j,imu)
         end do
         do i = 1,jthis%N1-1 
            ldiag(i,j) = jthis%ldiag(i,j,imu)
            udiag(i,j) = jthis%udiag(i,j,imu)
         end do
      end do
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC)
      do iphi = 0,geom%Nphi-1
         do j = 1,nth
            do i = 1,nr-1
               Acomp(i,j) = A(i,j-1,iphi)
            end do
         end do
        call fourrow(Acomp,1)
      
        !*** algorithm de Thomas for LU matrix system ***
        ! -> solving of the tridiagonal system only on 
        !     half the Fourier modes
        call thomas_algorithm(nr-1,midth+1,ldiag,diag,udiag,Acomp)
        ! -> copy of the results fot the Fourier modes '-j'
        do j = 1,midth+1
           if (j.ne.1 .and. j.ne.midth+1) then 
              !*** copy the result for the fourier mode '-j' ***
              do i = 1,nr-1
                 Acomp(i,nth+2-j) = conjg(Acomp(i,j))
              end do
           endif
        end do
      
        !*** Perform FFT 1D inverse ***
        call fourrow(Acomp,-1)
        do j = 0,nthm1
          do i = 1,nr-1
            Abar(i,j,iphi) = dreal(Acomp(i,j+1))/geom%Ntheta
          end do
        end do
      
        !***  Neumann boundary conditions ***
        do j = 0,nthm1
           Abar(0,j,iphi)  = Abar(1,j,iphi)
           Abar(nr,j,iphi) = Abar(nr-1,j,iphi)
        end do
      
        !*** periodic conditions in theta ***
        do i = 0,nr
           Abar(i,nth,iphi) = Abar(i,0,iphi)
        end do
      end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
      !*** periodic conditions in phi ***
      do j = 0,nth
        do i = 0,nr
          Abar(i,j,geom%Nphi) = Abar(i,j,geom%Nphi-1)
        end do
      end do
    end if
    call clck_time(eclock_compJ0)
    call clck_diff(bclock_compJ0,eclock_compJ0,global_time_compJ0)
!R3 call r3_info_end (r3_info_index_0)      !R3
  end subroutine compute_gyrofunction_3D
end module gyroaverage_class
