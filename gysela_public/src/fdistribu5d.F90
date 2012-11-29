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
! file : fdistribu5d.f90
! date : 16/09/2005
! - 5D distribution function which is assumed
!  - no-periodic in the first direction   (r),
!  - periodic in the second direction     (theta),
!  - periodic in the third direction      (phi),
!  - no-periodic in the fourth direction  (vparallel),
!  - no-periodic in the fifth direction   (mu).
!-------------------------------------------------------
module fdistribu5d_class
  use prec_const
  use MPIutils_module
  use mem_alloc_module
  use spline1d_class
  use geometry_class
  use init_profile_class
  use local_spline_module
  use fequil4d_module
      
  implicit none
  public :: new_f5D, del_f5D
      
  !*** f(r,theta,phi,vpar) parametrized by mu_val ***
  type :: fdistribu5d
    integer     :: muid
    real(RKIND) :: mu_value
    integer     :: n1, n2, n3, n4
    integer     :: istart, istart_buf, iend, iend_buf
    integer     :: jstart, jstart_buf, jend, jend_buf
    integer     :: istart_modif, iend_modif
    integer     :: rstart_modif, rend_modif
    real(RKIND) :: h1, h2, h3, h4
      
    !*** spline initialization ***
    integer                                 :: BCr_left
    integer                                 :: BCr_right
    integer                                 :: BCvpar_left
    integer                                 :: BCvpar_right
    type(nspline1d)                         :: nspline1d_r
    type(nspline1d), dimension (:), pointer :: nspline1d_vpar
    type(pspline1d)                         :: pspline1d_theta
    type(pspline1d)                         :: pspline1d_phi
    type(splinehh) , dimension(:), pointer  :: hhspline2d_rtheta
    type(splinehh), dimension(:), pointer   :: hhspline2d_vparphi
    type(splinehh) , dimension(:), pointer  :: uxspline2d_rtheta
    type(splinehh) , dimension(:), pointer  :: uyspline2d_rtheta
      
    ! values (4D array parametrized by a value of mu)
    real(RKIND), dimension(:,:,:,:), pointer  :: values
  end type fdistribu5d
      
  ! for the initial perturbation
  real(RKIND), private, &
    dimension(:,:), pointer :: fperturb_thetaphi
  real(RKIND), private, &
    dimension(:,:), pointer :: random_phase    
      
  !******************************
  contains
  !******************************       
      
  !-------------------------------------------------------
  ! check if the point must be ignore for the computation
  !  in the radial direction
  !------------------------------------------------------- 
  subroutine ignore_r_boundary(ir,bound_r)
    use globals, only : Nr
    integer, intent(in)  :: ir
    logical, intent(out) :: bound_r 
    
    bound_r = .false.
    if (ir.eq.0) then
        bound_r = .true.
    end if
    if (ir.eq.Nr) bound_r = .true.
  end subroutine ignore_r_boundary
      
  !----------------------------------------------------------
  ! Constructor of the 5D distribution function 
  !   f(r,theta,phi,vpar,mu)=feq(1+delta f)
  !   where feq is the equilibrium part and delta f the 
  !   perturbated part.
  !   This function is parallelized on all the processors.
  !----------------------------------------------------------
  subroutine new_f5d(fthis,geom)
    use globals, only : mu_id, istart, istart_buf, &
      iend, iend_buf, jstart, jstart_buf, jend, jend_buf, &
      dom_r, dom_theta, Nbthread, m, n, &
      BC_Hermite, transpose4D 
    type(fdistribu5d), intent(out) :: fthis
    type(geometry)   , intent(in)  :: geom
      
    real(RKIND) :: rmin_loc, thetamin_loc
    integer     :: ithread
      
    !*** dimension initialization ***
    fthis%n1         = geom%Nr
    fthis%n2         = geom%Ntheta
    fthis%n3         = geom%Nphi
    fthis%n4         = geom%Nvpar   
    fthis%h1         = geom%dr
    fthis%h2         = geom%dtheta
    fthis%h3         = geom%dphi
    fthis%h4         = geom%dvpar
    fthis%istart     = istart
    fthis%istart_buf = istart_buf
    fthis%jstart     = jstart
    fthis%jstart_buf = jstart_buf
    fthis%iend       = iend
    fthis%iend_buf   = iend_buf
    fthis%jend       = jend
    fthis%jend_buf   = jend_buf
    fthis%muid       = mu_id
    fthis%mu_value   = 0
    if (.not.memory_test) &
      fthis%mu_value   = geom%mug(mu_id)
       
    !*** boundary conditions initialization ***
      fthis%istart_modif = max(1,fthis%istart)
      fthis%iend_modif   = min(geom%Nr-1,fthis%iend)
      fthis%rstart_modif = 1
      fthis%rend_modif   = geom%Nr-1
      fthis%BCr_left     = BC_Hermite
      fthis%BCr_right    = BC_Hermite
    ! Neumann boundary conditions in vpar direction
    fthis%BCvpar_left    = 1
    fthis%BCvpar_right   = 1
      
    !*** array allocation for distribution function values ***
    call glob_allocate(fthis%values, &
      fthis%istart_buf,fthis%iend_buf,fthis%jstart_buf, &
      fthis%jend_buf,0,fthis%n3,0,fthis%n4, &
      'fthis%values')
      
    !*** spline initialization ***
    call new_spline1d_natural(fthis%nspline1d_r,fthis%n1,fthis%h1)
    call new_spline1d_period(fthis%pspline1d_theta, &
      fthis%n2,fthis%h2)
    call new_spline1d_period(fthis%pspline1d_phi,fthis%n3,fthis%h3)
    allocate(fthis%nspline1d_vpar(1:Nbthread))
    allocate(fthis%hhspline2d_rtheta(1:Nbthread))
      allocate(fthis%uxspline2d_rtheta(1:Nbthread))
      allocate(fthis%uyspline2d_rtheta(1:Nbthread))
    do ithread = 1,Nbthread
      call new_spline1d_natural(fthis%nspline1d_vpar(ithread), &
        fthis%n4,fthis%h4)
      if (transpose4D) then
        rmin_loc     = geom%rmin + (-1)*fthis%h1
        thetamin_loc = geom%thetamin + (-1)*fthis%h2
          call new_splinehh(fthis%hhspline2d_rtheta(ithread), &
            fthis%n1+2,fthis%n2+1,rmin_loc,thetamin_loc, &
            fthis%h1,fthis%h2)
          call new_splinehh(fthis%uxspline2d_rtheta(ithread), &
            fthis%n1+2,fthis%n2+1,rmin_loc,thetamin_loc, &
            fthis%h1,fthis%h2)
          call new_splinehh(fthis%uyspline2d_rtheta(ithread), &
            fthis%n1+2,fthis%n2+1,rmin_loc,thetamin_loc, &
            fthis%h1,fthis%h2)
      else
        rmin_loc     = geom%rmin + (istart-1)*fthis%h1
        thetamin_loc = geom%thetamin + (jstart-1)*fthis%h2
        call new_splinehh(fthis%hhspline2d_rtheta(ithread), &
          dom_r+1,dom_theta+1,rmin_loc,thetamin_loc, &
          fthis%h1,fthis%h2)
        call new_splinehh(fthis%uxspline2d_rtheta(ithread), &
          dom_r+1,dom_theta+1,rmin_loc,thetamin_loc, &
          fthis%h1,fthis%h2)
        call new_splinehh(fthis%uyspline2d_rtheta(ithread), &
          dom_r+1,dom_theta+1,rmin_loc,thetamin_loc, &
          fthis%h1,fthis%h2)
      end if
    enddo
      
    !*** allocation for the initial perturbation ***
    call glob_allocate(random_phase,1,m, &
      1,n,'random_phase')
    call glob_allocate(fperturb_thetaphi,0,fthis%n2, &
      0,fthis%n3,'fperturb_thetaphi')
  end subroutine new_f5d
      
  
  !--------------------------------- 
  ! 5D spline destruction 
  !---------------------------------
  subroutine del_f5d(fthis)
    use globals, only : Nbthread
    type(fdistribu5d), intent(out) :: fthis
      
    integer :: ithread
      
    call glob_deallocate(fthis%values)
      
    call del_spline1d_natural(fthis%nspline1d_r)
    call del_spline1d_period(fthis%pspline1d_theta)
    call del_spline1d_period(fthis%pspline1d_phi)
    do ithread = 1,Nbthread
      call del_spline1d_natural(fthis%nspline1d_vpar(ithread))
      call del_splinehh(fthis%hhspline2d_rtheta(ithread))
        call del_splinehh(fthis%uxspline2d_rtheta(ithread))
        call del_splinehh(fthis%uyspline2d_rtheta(ithread))
    enddo
    deallocate(fthis%nspline1d_vpar)
    deallocate(fthis%hhspline2d_rtheta)
      deallocate(fthis%uxspline2d_rtheta)
      deallocate(fthis%uyspline2d_rtheta)
    call glob_deallocate(random_phase)
    call glob_deallocate(fperturb_thetaphi)
  end subroutine del_f5d
      
  !************************************************
  ! AT THE FIRST STEP THE DISTRIBUTION FUNCTION 
  !  IS GIVEN BY THE INITIAL CONDITIONS
  !************************************************
  !------------------------------------------------------------
  ! Initialisation of the perturbed part (delta f) of the 
  !   distribution function
  !   delta_f(theta,phi) = epsilon*cos(k*phi+m*theta)) or
  !         = sum_{m,n} epsilon*cos(k*phi+m*theta+pert_{m,n})
  !  where pert_{m,n} has random values
  !------------------------------------------------------------
  subroutine init_fperturb(geom,init_prof,Sfperturb_thetaphi)
    use globals, only : init_perturb, m, n, &
      epsilon, Rosenbluth_Hinton
    type(geometry)               , intent(in)  :: geom
    type(init_profile)           , intent(in)  :: init_prof
    real(RKIND), dimension(0:,0:), intent(out) :: Sfperturb_thetaphi
    
    real(RKIND) :: toroidal_ratio, fperturb_tmp
    integer     :: itheta, iphi
    integer     :: m_max, n_max, in, im, in_min, in_max
      
    !*** f initialisation *** 
    toroidal_ratio = TWOPI/geom%Lphi
    if (Rosenbluth_Hinton) then
      do iphi = 0,geom%Nphi
        do itheta = 0,geom%Ntheta
          Sfperturb_thetaphi(itheta,iphi) = epsilon
        end do
      end do
    else
      select case (init_perturb)
      case (1)
        ! -> initialization with the mode (m,n)
        do iphi = 0,geom%Nphi
          do itheta = 0,geom%Ntheta      
            Sfperturb_thetaphi(itheta,iphi) = &
              epsilon*cos(m*geom%thetag(itheta) - &
              n*geom%phig(iphi)*toroidal_ratio)
          end do
        end do
      case (2)
        ! -> initialization with several modes 
        !    (1 < im < m) and (1 < in < n) 
        m_max = m
        n_max = n
        do iphi = 0,geom%Nphi
          do itheta = 0,geom%Ntheta
            fperturb_tmp = 0._RKIND
            do in = 1,n_max
              do im = 1,m_max
                fperturb_tmp = fperturb_tmp + &
                  cos(im*geom%thetag(itheta) - &
                  in*geom%phig(iphi)*toroidal_ratio + &
                  random_phase(im,in))
              end do
            end do
            fperturb_tmp = epsilon*fperturb_tmp / &
              real(m_max*n_max)
            if ( abs(fperturb_tmp).le.1.e-30_RKIND ) &
              fperturb_tmp = 1.e-30_RKIND
            Sfperturb_thetaphi(itheta,iphi) = fperturb_tmp
          end do
        end do
      case (3)
        ! -> initialization with several m modes 
        !     for one unique n mode
        in    = 1
        m_max = m
        do iphi = 0,geom%Nphi
          do itheta = 0,geom%Ntheta
            fperturb_tmp = 0._RKIND        
            do im = 1,m_max
              fperturb_tmp = fperturb_tmp + &
                cos(im*geom%thetag(itheta) - &
                n*geom%phig(iphi)*toroidal_ratio + &
                random_phase(im,in))
            end do
            fperturb_tmp = epsilon*fperturb_tmp
            if ( abs(fperturb_tmp).le.1.e-30_RKIND ) &
              fperturb_tmp = 1.e-30_RKIND
            Sfperturb_thetaphi(itheta,iphi) = fperturb_tmp
          end do
        end do
      case (4)
        ! -> initialization with several n modes 
        !     for one unique m mode
        im    = 1
        n_max = n
        do iphi = 0,geom%Nphi
          do itheta = 0,geom%Ntheta
            fperturb_tmp = 0._RKIND        
            do in = 1,n_max
              fperturb_tmp = fperturb_tmp + &
                cos(m*geom%thetag(itheta) - &
                in*geom%phig(iphi)*toroidal_ratio + &
                random_phase(im,in))
            end do
            fperturb_tmp = epsilon*fperturb_tmp
            if ( abs(fperturb_tmp).le.1.e-30_RKIND ) &
              fperturb_tmp = 1.e-30_RKIND
            Sfperturb_thetaphi(itheta,iphi) = fperturb_tmp
          end do
        end do
      case default
        if (pglobal_id.eq.0) then
          print*, 'init_perturb = ', &
            init_perturb, ' must be 1, 2, 3 or 4'
          stop
        end if        
      end select
    end if
  end subroutine init_fperturb
      
  !-------------------------------------------------------------
  ! Calculate another form of the initial distribution function
  !  which is an equilibrium solution of the gyrokinetic
  !  equation with collisions when phi is small and gradTi = 0
  !-------------------------------------------------------------
  subroutine init_fPhi_analyt(geom,init_prof,init_magnet, &
    ir,itheta,iphi,ivpar,imu,fanalyt)
    use globals, only : init_perturb, ipow
    use bessel_module, only : compute_modified_bessel
    use init_magnetic_class
    type(geometry)     , intent(in)  :: geom
    type(init_profile) , intent(in)  :: init_prof
    type(init_magnetic), intent(in)  :: init_magnet
    integer            , intent(in)  :: ir, itheta, iphi
    integer            , intent(in)  :: ivpar, imu
    real(RKIND)        , intent(out) :: fanalyt
      
    real(RKIND)    :: rmin, rmax
    real(RKIND)    :: denom, Ti_r
    real(RKIND)    :: energy, phi_init
    real(RKIND)    :: hamiltonian, neq0, neq1, Dn, npic
    real(RKIND)    :: CC, CC2, rmaxmod
    complex(CKIND) :: I00, I01, I02, nouseI0
    complex(CKIND) :: K00, K01, K02, nouseK0
    complex(CKIND) :: cst1, cst2, crmin, crmaxmod, cCC2
    integer        :: ii
      
    !*** f initialisation *** 
    rmin       = geom%rg(0)
    rmax       = geom%rg(geom%Nr)
    Ti_r       = init_prof%Ti(ir)
    denom      = sqrt(2._RKIND*PI*Ti_r)*(2._RKIND*PI*Ti_r)**ipow
    ! finding the extrema of the initial profile
    do ii = 0,geom%Nr
      neq0 = -1.e10_RKIND
      neq1 = 1.e10_RKIND
      if (neq0 < init_prof%n0(ii)) then
        neq0  = init_prof%n0(0)
      end if
      if (neq1 > init_prof%n0(ii)) then 
        neq1  = init_prof%n0(geom%Nr)
      end if
    end do
    Dn        = HF*abs(neq1 - neq0)
    npic      = HF*(neq1 + neq0)
    !*** LIMIT CONDITIONS ***
    CC        = sqrt( (npic+Dn)/(npic-Dn) )
    rmaxmod   = rmax*CC
    ! for input in Bessel function
    crmin     = dcmplx(rmin)
    crmaxmod  = dcmplx(rmaxmod)
    call compute_modified_bessel(-crmin,nouseI0,K00)
    call compute_modified_bessel(crmaxmod,I00,nouseK0)
    call compute_modified_bessel(-crmaxmod,nouseI0,K01)
    call compute_modified_bessel(crmin,I01,nouseK0)
    cst1 = -TW*Dn*K00 / ( (npic+Dn)*( I00*K00 - I01*K01) )
    cst2 =  TW*Dn*I01 / ( (npic+Dn)*( I00*K00 - I01*K01) )
    !*** GENERAL EXPRESSION ***
    CC2  = geom%rg(ir)*sqrt( (npic+Dn)/init_prof%n0(ir) )
    ! for input in Bessel function
    cCC2 = dcmplx(CC2)
    call compute_modified_bessel(-cCC2,nouseI0,K02)
    call compute_modified_bessel(cCC2,I02,nouseK0)
    phi_init    = dble(- log( init_prof%n0(ir)/(npic+Dn) ) + &
      cst1*I02 + cst2*K02 )
    energy      = 0.5_RKIND*geom%vparg(ivpar)**2 + &
      geom%mug(imu)*init_magnet%B_norm(ir,itheta)
    hamiltonian = energy + phi_init
    ! what is here called "fanalyt" is in fact what 
    ! we want as initial analytical solution for f
    fanalyt     = (npic+Dn)/denom * exp(-hamiltonian/Ti_r)
      
    if ( abs(fanalyt).le.1.e-30_RKIND ) fanalyt = 1.e-30_RKIND
  end subroutine init_fPhi_analyt
      
  !------------------------------------------------------------- 
  ! Initialisation of the distribution function
  !  f(r,theta,phi,vpar,mu) = 
  !    feq(r,theta,vpar,mu)*(1+delta f(r,theta,phi,vpar))
  !  where delta f = g(r)*h(vpar)*delta p(theta,phi) with
  !   - g(r) is defined with polynomial functions and
  !   - h(vpar) is defined with a Gaussian shape
  !  such that g and h vanishes at both radial boundaries.
  !   - delta p = epsilon*cos(k*phi+m*theta)) or
  !             = sum_{m,n} epsilon*cos(k*phi+m*theta+pert_{m,n})
  !  where pert_{m,n} has random values
  !-------------------------------------------------------------
  subroutine f_initialisation(fthis,geom,init_prof,Sfmu_eq)
    use globals, only : mu_id, m, n, epsilon, &
      Rosenbluth_Hinton, hvpar_in_fperturb, &
      deltarTi, profile_choice, uout_res, &
      Rarray1_Nr, Rarray1_Nvpar, Rarray1_NrNthetaNphi, &
      outputproc
    use utils_module, only : init_random_seed
    include "mpiperso.h"
!R3 #include "r3_info.h" !R3
    type(fdistribu5d)            , intent(inout) :: fthis
    type(geometry)               , intent(in)    :: geom
    type(init_profile)           , intent(in)    :: init_prof
    real(RKIND), dimension(:,:,:), pointer       :: Sfmu_eq
      
    integer     :: i1, i2, i3, i4, ir, ivpar
    integer     :: i, j, m_max, n_max, ierr
    real(RKIND) :: rmin, rmax, vparmin, vparmax
    real(RKIND) :: max_gr_pol, rbuff, lambda, gr_pol_tmp, hvpar_tmp
    real(RKIND) :: Lr_loc, Lvpar_loc
    real(RKIND) :: feq_value, fperturb_tmp
    real(RKIND), dimension(:)  , pointer :: gr_pol, hvpar
      
    if ((.not.Rosenbluth_Hinton).and.((n .eq. 0) .or. (m .eq. 0))) then
      !*** there is normally no perturbation, i.e. epsilon = 0
      if ((epsilon .ne. 0._RKIND) .and. &
        (pglobal_id .eq. outputproc)) then
        write (uout_res,*) '---> epsilon is set to 0 ', &
          'as n or m is equal to 0' 
        write (uout_res,*) ' '
      end if
      !*** computation of the initial values for f ***
!$OMP PARALLEL private(i1,i2,i3,i4,feq_value) default(shared)
!$OMP DO SCHEDULE(static)
    do i4 = 0,fthis%n4
      do i3 = 0,fthis%n3
        do i2 = fthis%jstart,fthis%jend
          do i1 = fthis%istart,fthis%iend
            feq_value  = Sfmu_eq(i1,i2,i4)
            fthis%values(i1,i2,i3,i4) = feq_value
          end do
        end do
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL
      
  else
      gr_pol => Rarray1_Nr
      hvpar  => Rarray1_Nvpar
      
      !*** computation of the initial perturbation delta f *** 
!R3 call r3_info_begin (r3_info_index_0, 'CPU_f_initialisation') !R3
      !-> computation of g(r)
      Lr_loc = geom%Lr
      rmin   = geom%rg(0)
      rmax   = geom%rg(geom%Nr)
      if (Rosenbluth_Hinton) then
        do ir = 0,geom%Nr
          gr_pol(ir) = &
            sin(TWOPI*(geom%rg(ir)-rmin)/(2._RKIND*Lr_loc))
        end do
      else
        if (profile_choice.ne.6) then
          max_gr_pol = (0.5_RKIND*Lr_loc)**6 * &
            (rmin-rmax+0.5_RKIND*Lr_loc)**6
          do ir = 0,geom%Nr
            gr_pol(ir) = (geom%rg(ir)-rmin)**6 * &
              (geom%rg(ir)-rmax)**6/max_gr_pol
          end do
        else
          lambda = deltarTi*Lr_loc
          rbuff  = 0.15_RKIND*Lr_loc
          do ir = 0,geom%Nr
            gr_pol(ir) = 0.5_RKIND * &
              (tanh((geom%rg(ir)-rmin-rbuff)/lambda) - &
              tanh((geom%rg(ir)-rmax+rbuff)/lambda))
          end do
        end if
      end if
      
      !-> computation of h(vpar)
      if (hvpar_in_fperturb) then
        vparmin    = geom%vparg(0)
        vparmax    = geom%vparg(geom%Nvpar)
        Lvpar_loc  = geom%Lvpar
        do ivpar = 0,geom%Nvpar
          hvpar(ivpar) = exp(-((geom%vparg(ivpar)-&
            (vparmin+Lvpar_loc/2._RKIND))/(0.1*Lvpar_loc))**2)
        end do
      else
        do ivpar = 0,geom%Nvpar
          hvpar(ivpar) = 1._RKIND
        end do
      end if
      
      !-> computation of the random phase
      m_max = m 
      n_max = n
      call init_random_seed()
      if (pglobal_id .eq. 0) then
        do j = 1,n_max
          do i = 1,m_max
            call random_number(random_phase(i,j))
          end do
        end do
      end if
      call MPI_BCAST(random_phase(1,1),n_max*m_max, &
        MPI_REAL8,0,MPI_COMM_WORLD,ierr)
      random_phase = TWOPI * random_phase
      
      !-> computation of the perturbation delta_p(theta,phi)
      call init_fperturb(geom,init_prof,fperturb_thetaphi)
      
      !*** computation of the initial values for f ***
!$OMP PARALLEL private(i1,i2,i3,i4,feq_value,gr_pol_tmp,hvpar_tmp, &
!$OMP  fperturb_tmp) default(shared)
!$OMP DO SCHEDULE(static)
      do i4 = 0,fthis%n4
        hvpar_tmp = hvpar(i4)
        do i3 = 0,fthis%n3
          do i2 = fthis%jstart,fthis%jend
            fperturb_tmp = fperturb_thetaphi(i2,i3)
            do i1 = fthis%istart,fthis%iend
              gr_pol_tmp = gr_pol(i1)
              feq_value  = Sfmu_eq(i1,i2,i4)
              fthis%values(i1,i2,i3,i4) = feq_value*(1._RKIND + &
                fperturb_tmp*gr_pol_tmp*hvpar_tmp)
            end do
          end do
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
!R3 call r3_info_end (r3_info_index_0) !R3
      
    end if
      
  end subroutine f_initialisation
  
      
  !************************************************
  ! CUBIC SPLINE INTERPOLATION
  !************************************************
  !------------------------------------------------------- 
  ! Compute the 2D cubic spline in the directions r and
  !  theta. These spline coefficients are used for the
  !  shift in (r,theta)
  !-------------------------------------------------------     
  subroutine compute_spline_rtheta(fthis,rhs,scoef2d)
    use globals, only : Rarray1_m1Nrp1m1Nthetap1
    type(fdistribu5d)              , intent(inout) :: fthis
    real(RKIND), dimension(0:,0:)  , intent(in)    :: rhs
    real(RKIND), dimension(-1:,-1:), intent(inout) :: scoef2d
      
    real(RKIND), dimension(:,:), pointer :: gamma	
    integer                              :: i1, i2
      
    scoef2d = 0._RKIND
      
    !*** solving in r direction with Dirichlet conditions ***
    gamma => Rarray1_m1Nrp1m1Nthetap1
      
    do i2 = 0,fthis%n2
      call natural_spline_coef(fthis%nspline1d_r,rhs(0:,i2),&
        fthis%BCr_left,fthis%BCr_right)
      do i1 = -1,fthis%n1+1
        gamma(i1,i2) = fthis%nspline1d_r%scoef(i1)
      end do
    end do
    
    !*** solving in theta direction with ***
    !***  periodic boundary conditions      ***
    do i1=-1,fthis%n1+1
      call period_spline_coef(fthis%pspline1d_theta,gamma(i1,0:))
      do i2 = -1,fthis%n2+1
        scoef2d(i1,i2) = fthis%pspline1d_theta%scoef(i2)
      end do
    end do
  end subroutine compute_spline_rtheta
      
end module fdistribu5d_class
