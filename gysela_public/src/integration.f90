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
! file : integration.f90
! date : 07/06/2006
! - computation of the integrals with cubic splines or
!  with a collocation method
!-------------------------------------------------------
module integration_module
  use prec_const
  use geometry_class
  use spline1D_class
  use OMPutils_module, only : omp1dvector
  
  implicit none
      
  !-------------------------------------
  !  PRIVATE VARIABLES
  !-------------------------------------
  !-> used for OpenMP parallelisation
  type(omp1dvector), dimension(:), pointer, &
    private :: Romp_funcr_0Nr
  type(omp1dvector), dimension(:), pointer, &
    private :: Romp_intrdr_0Ntheta
      
  !******************************
  contains
  !******************************
      
  !-------------------------------------------------
  ! Constructor 
  !-------------------------------------------------
  subroutine new_integration(geom)
    use globals, only : Nbthread
    type(geometry), intent(in) :: geom
      
    integer :: tid
      
    !*** array allocation for OpenMP parallelization ***
    allocate(Romp_funcr_0Nr(1:Nbthread))
    allocate(Romp_intrdr_0Ntheta(1:Nbthread))
    do tid = 1,Nbthread
      call glob_allocate(Romp_funcr_0Nr(tid)%val,0,geom%Nr)
      call glob_allocate(Romp_intrdr_0Ntheta(tid)%val,0,geom%Ntheta)
    end do
  end subroutine new_integration
      
  !-------------------------------------------------
  ! Destructor 
  !-------------------------------------------------
  subroutine delete_integration
    use globals, only : Nbthread
      
    integer :: tid
      
    !*** array deallocation for OpenMP parallelization ***
    do tid = 1,Nbthread
      call glob_deallocate(Romp_funcr_0Nr(tid)%val)
      call glob_deallocate(Romp_intrdr_0Ntheta(tid)%val)
    end do
    deallocate(Romp_funcr_0Nr)
    deallocate(Romp_intrdr_0Ntheta)
  end subroutine delete_integration
      
  !***********************************************
  !  INTEGRATION WITH CUBIC SPLINES
  !***********************************************
  !---------------------------------------------------------
  ! Computes for all function H(x) the integral \int(H(x))dx
  ! where H(x) is a periodic function
  !---------------------------------------------------------
  subroutine compute_phase_integral_pCS(Nx,dx, &
    Hx,Hx_pspline1d,integral_value)
    use spline1d_class, only : integration_coef_boundaries
    integer                   , intent(in)    :: Nx
    real(RKIND)               , intent(in)    :: dx
    real(RKIND), dimension(0:), intent(in)    :: Hx
    type(pspline1d)           , intent(inout) :: Hx_pspline1d
    real(RKIND)               , intent(out)   :: integral_value
    
    integer     :: idx, ix
    real(RKIND) :: Hx_int_tmp
      
    !*** case for periodic boundary condition ***
    call period_spline_coef(Hx_pspline1d,Hx)
      
    !*** \int(H(x))dx computation ***
    Hx_int_tmp = sum(Hx_pspline1d%scoef(2:Nx-2))
    Hx_int_tmp = 6._RKIND*dx*Hx_int_tmp
    do idx = 1,6
      Hx_int_tmp = Hx_int_tmp + Hx_pspline1d%factor_int(idx)*&
        Hx_pspline1d%scoef(Hx_pspline1d%indx(idx))
    end do
    integral_value = Hx_int_tmp
  end subroutine compute_phase_integral_pCS
      
  !********************************************************
  !  COMPUTATION OF INTEGRALS WITH CUBIC SPLINES
  !********************************************************
  !---------------------------------------------------------
  ! Computes the integral \int(H(r))dr where 
  !   H(r) is a non-periodic function
  !---------------------------------------------------------
  subroutine compute_r_integral_CS(Nr,dr, &
    Hr,BCr_left,BCr_right,Hr_nspline1d, &
    integral_value,Sderiv_Hr)
    use spline1d_class, only : integration_coef_boundaries
    integer                    , intent(in)    :: Nr
    real(RKIND)                , intent(in)    :: dr
    real(RKIND), dimension(0:) , intent(in)    :: Hr
    integer                    , intent(in)    :: BCr_left
    integer                    , intent(in)    :: BCr_right
    type(nspline1d)            , intent(inout) :: Hr_nspline1d
    real(RKIND)                , intent(out)   :: integral_value
    real(RKIND), dimension(0:1), intent(in), & 
                                      optional :: Sderiv_Hr
    
    real(RKIND) :: Hr_int_tmp
    integer     :: idr, ir
      
    !*** case for natural boundary condition ***
    call natural_spline_coef(Hr_nspline1d,Hr, &
      BCr_left,BCr_right,Sderiv_Hr)
      
    !*** \int(H(r))dr computation ***
    Hr_int_tmp = sum(Hr_nspline1d%scoef(2:Nr-2))
    Hr_int_tmp = 6._RKIND*dr*Hr_int_tmp
    do idr = 1,6
      Hr_int_tmp = Hr_int_tmp + &
        Hr_nspline1d%factor_int(idr)*&
        Hr_nspline1d%scoef(Hr_nspline1d%indx(idr))
    end do
    integral_value = Hr_int_tmp
  end subroutine compute_r_integral_CS
      
  !---------------------------------------------------------
  ! Computes for all function H(theta) the integral 
  !  \int(H(theta))dtheta where H(theta) is a 
  !  periodic function
  !---------------------------------------------------------
  subroutine compute_theta_integral_CS(Ntheta,dtheta, &
    Htheta,Htheta_pspline1d,integral_value)
    use spline1d_class, only : integration_coef_boundaries
    integer                   , intent(in)    :: Ntheta
    real(RKIND)               , intent(in)    :: dtheta
    real(RKIND), dimension(0:), intent(in)    :: Htheta
    type(pspline1d)           , intent(inout) :: Htheta_pspline1d
    real(RKIND)               , intent(out)   :: integral_value
    
    real(RKIND) :: Htheta_int_tmp
    integer     :: idtheta, itheta
      
    !*** case for periodic boundary condition ***
    call period_spline_coef(Htheta_pspline1d,Htheta)
      
    !*** \int(H(theta))dtheta computation ***
    Htheta_int_tmp = sum(Htheta_pspline1d%scoef(2:Ntheta-2))
    Htheta_int_tmp = 6._RKIND*dtheta*Htheta_int_tmp
    do idtheta = 1,6
      Htheta_int_tmp = Htheta_int_tmp + &
        Htheta_pspline1d%factor_int(idtheta)*&
        Htheta_pspline1d%scoef(Htheta_pspline1d%indx(idtheta))
    end do
    integral_value = Htheta_int_tmp
  end subroutine compute_theta_integral_CS
      
  !---------------------------------------------------------
  ! Computes for all function H(phi) the integral 
  !  \int(H(phi))dphi where H(phi) is a 
  !  periodic function
  !---------------------------------------------------------
  subroutine compute_phi_integral_CS(Nphi,dphi, &
    Hphi,Hphi_pspline1d,integral_value)
    use spline1d_class, only : integration_coef_boundaries
    integer                   , intent(in)    :: Nphi
    real(RKIND)               , intent(in)    :: dphi
    real(RKIND), dimension(0:), intent(in)    :: Hphi
    type(pspline1d)           , intent(inout) :: Hphi_pspline1d
    real(RKIND)               , intent(out)   :: integral_value
    
    real(RKIND) :: Hphi_int_tmp
    integer     :: idphi, iphi
      
    !*** case for periodic boundary condition ***
    call period_spline_coef(Hphi_pspline1d,Hphi)
      
    !*** \int(H(phi))dphi computation ***
    Hphi_int_tmp = sum(Hphi_pspline1d%scoef(2:Nphi-2))
    Hphi_int_tmp = 6._RKIND*dphi*Hphi_int_tmp
    do idphi = 1,6
      Hphi_int_tmp = Hphi_int_tmp + &
        Hphi_pspline1d%factor_int(idphi)*&
        Hphi_pspline1d%scoef(Hphi_pspline1d%indx(idphi))
    end do
    integral_value = Hphi_int_tmp
  end subroutine compute_phi_integral_CS
      
  !---------------------------------------------------------
  ! Computes the integral \int(H(vpar))dvpar where 
  !   H(vpar) is a non-periodic function
  !---------------------------------------------------------
  subroutine compute_vpar_integral_CS(Nvpar,dvpar, &
    Hvpar,BCvpar_left,BCvpar_right,Hvpar_nspline1d, &
    integral_value,Sderiv_Hvpar)
    use spline1d_class, only : integration_coef_boundaries
    integer                    , intent(in)    :: Nvpar
    real(RKIND)                , intent(in)    :: dvpar
    real(RKIND), dimension(0:) , intent(in)    :: Hvpar
    integer                    , intent(in)    :: BCvpar_left
    integer                    , intent(in)    :: BCvpar_right
    type(nspline1d)            , intent(inout) :: Hvpar_nspline1d
    real(RKIND)                , intent(out)   :: integral_value
    real(RKIND), dimension(0:1), intent(in), & 
                                      optional :: Sderiv_Hvpar
    
    real(RKIND) :: Hvpar_int_tmp
    integer     :: idvpar, ivpar
      
    !*** case for natural boundary condition ***
    call natural_spline_coef(Hvpar_nspline1d,Hvpar, &
      BCvpar_left,BCvpar_right,Sderiv_Hvpar)
      
    !*** \int(H(vpar))dvpar computation ***
    Hvpar_int_tmp = sum(Hvpar_nspline1d%scoef(2:Nvpar-2))
    Hvpar_int_tmp = 6._RKIND*dvpar*Hvpar_int_tmp
    do idvpar = 1,6
      Hvpar_int_tmp = Hvpar_int_tmp + &
        Hvpar_nspline1d%factor_int(idvpar)*&
        Hvpar_nspline1d%scoef(Hvpar_nspline1d%indx(idvpar))
    end do
    integral_value = Hvpar_int_tmp
  end subroutine compute_vpar_integral_CS
      
  !--------------------------------------------------
  ! Computes for all function H(r,theta,phi) the integral
  !  I(r) = \int(H(r,theta,phi)dthetadphi
  !--------------------------------------------------
  subroutine compute_intdthetadphi_CS(H,geom, &
    pspline1d_theta,pspline1d_phi,H_intdthetadphi)
    use globals, only : Rarray_NrNphi
    real(RKIND), &
      dimension(0:,0:,0:), intent(in)    :: H
    type(geometry)       , intent(in)    :: geom
    type(pspline1d)      , intent(inout) :: pspline1d_theta
    type(pspline1d)      , intent(inout) :: pspline1d_phi
    real(RKIND), &
      dimension(0:)      , intent(out)   :: H_intdthetadphi
      
    integer                              :: ir, itheta, iphi
    real(RKIND)                          :: H_inttheta_tmp
    real(RKIND), dimension(:,:), pointer :: H_inttheta
      
    H_inttheta => Rarray_NrNphi
      
    do ir = 0,geom%Nr
      do iphi = 0,geom%Nphi
        !*** integral according to theta direction ***
        call compute_theta_integral_CS(geom%Ntheta,geom%dtheta, &
          H(ir,0:,iphi),pspline1d_theta,H_inttheta_tmp)
        H_inttheta(ir,iphi) = H_inttheta_tmp
      end do
      !*** integral according to phi direction ***
      call compute_phi_integral_CS(geom%Nphi,geom%dphi, &
        H_inttheta(ir,0:),pspline1d_phi,H_intdthetadphi(ir))
    end do
  end subroutine compute_intdthetadphi_CS
      
  !--------------------------------------------------------
  ! Computes for all function H(r,theta,phi) the integral
  !  \int(H(r,theta,phi)rdrdthetadphi
  !--------------------------------------------------------
  subroutine compute_intrdrdthetadphi_CS(H,geom, &
    BCr_left,BCr_right,nspline1d_r,pspline1d_theta, &
    pspline1d_phi,H_intrdrdthdphi)
    use globals, only : Rarray1_NrNtheta, &
      Rarray1_Nr, Rarray1_Ntheta
    real(RKIND), &
      dimension(0:,0:,0:), intent(in)    :: H
    type(geometry)       , intent(in)    :: geom
    integer              , intent(in)    :: BCr_left, BCr_right
    type(nspline1d)      , intent(inout) :: nspline1d_r
    type(pspline1d)      , intent(inout) :: pspline1d_theta
    type(pspline1d)      , intent(inout) :: pspline1d_phi
    real(RKIND)          , intent(out)   :: H_intrdrdthdphi
      
    integer                              :: ir, itheta, idr
    real(RKIND), dimension(:,:), pointer :: H_intdphi
    real(RKIND), dimension(:), pointer   :: rhs_r
    real(RKIND), dimension(:), pointer   :: intrdr
    real(RKIND)                          :: intrdr_tmp
      
    H_intdphi => Rarray1_NrNtheta
    rhs_r     => Rarray1_Nr
    intrdr    => Rarray1_Ntheta
      
    !*** integral according to phi direction ***
    do itheta = 0,geom%Ntheta
      do ir = 0,geom%Nr
        call compute_phi_integral_CS(geom%Nphi,geom%dphi, &
          H(ir,itheta,0:geom%Nphi), &
          pspline1d_phi,H_intdphi(ir,itheta))
      end do
    end do
      
    !*** integral according to r direction ***
    do itheta = 0,geom%Ntheta-1
      do ir = 0,geom%Nr
        rhs_r(ir) = geom%rg(ir) * H_intdphi(ir,itheta)
      end do
      call natural_spline_coef(nspline1d_r,rhs_r,BCr_left,BCr_right)
      !*** \int(H(r))dr computation ***
      intrdr_tmp = sum(nspline1d_r%scoef(2:geom%Nr-2))
      intrdr_tmp = 6._RKIND*geom%dr*intrdr_tmp
      do idr = 1,6
        intrdr_tmp = intrdr_tmp + nspline1d_r%factor_int(idr)*&
          nspline1d_r%scoef(nspline1d_r%indx(idr))
      end do
      intrdr(itheta) = intrdr_tmp
    end do
      
    !*** periodic conditions in theta ****
    intrdr(geom%Ntheta) = intrdr(0)
      
    !*** integral according to theta direction ***
    call compute_theta_integral_CS(geom%Ntheta,geom%dtheta, &
      intrdr(0:geom%Ntheta),pspline1d_theta,H_intrdrdthdphi)
  end subroutine compute_intrdrdthetadphi_CS
      
  !********************************************************
  !  COMPUTATION OF INTEGRALS WITH CUBIC SPLINES + OPENMP
  !********************************************************
  !---------------------------------------------------------
  ! Computes the integral \int(H(r))dr where 
  !   H(r) is a non-periodic function
  !---------------------------------------------------------
  subroutine compute_omp_r_integral_CS(Nr,dr, &
    Hr,BCr_left,BCr_right,Hr_nspline1d, &
    scoef,integral_value,Sderiv_Hr)
    use spline1d_class, only : integration_coef_boundaries
    integer                     , intent(in)    :: Nr
    real(RKIND)                 , intent(in)    :: dr
    real(RKIND), dimension(:)   , pointer       :: Hr
    integer                     , intent(in)    :: BCr_left
    integer                     , intent(in)    :: BCr_right
    type(nspline1d)             , intent(inout) :: Hr_nspline1d
    real(RKIND) , dimension(:)  , pointer       :: scoef
    real(RKIND)                 , intent(out)   :: integral_value
    real(RKIND), &
      dimension(0:1), intent(in), optional      :: Sderiv_Hr
    
    real(RKIND) :: Hr_int_tmp
    integer     :: idr, ir
      
    !*** case for natural boundary condition ***
    call natural_omp_spline_coef(Hr_nspline1d,Hr, &
      BCr_left,BCr_right,scoef,Sderiv_Hr)
      
    !*** \int(H(r))dr computation ***
    Hr_int_tmp = sum(scoef(2:Nr-2))
    Hr_int_tmp = 6._RKIND*dr*Hr_int_tmp
    do idr = 1,6
      Hr_int_tmp = Hr_int_tmp + Hr_nspline1d%factor_int(idr)* &
        scoef(Hr_nspline1d%indx(idr))
    end do
    integral_value = Hr_int_tmp
  end subroutine compute_omp_r_integral_CS
      
  !---------------------------------------------------------
  ! Computes for all function H(theta) the integral 
  !  \int(H(theta))dtheta where H(theta) is a 
  !  periodic function
  !---------------------------------------------------------
  subroutine compute_omp_theta_integral_CS(Ntheta,dtheta, &
    Htheta,Htheta_pspline1d,scoef,integral_value)
    use spline1d_class, only : integration_coef_boundaries
    integer                  , intent(in)    :: Ntheta
    real(RKIND)              , intent(in)    :: dtheta
    type(pspline1d)          , intent(inout) :: Htheta_pspline1d
    real(RKIND), dimension(:), pointer       :: Htheta
    real(RKIND), dimension(:), pointer       :: scoef
    real(RKIND)              , intent(out)   :: integral_value
    
    real(RKIND) :: Htheta_int_tmp
    integer     :: idtheta, itheta
      
    !*** case for periodic boundary condition ***
    call period_omp_spline_coef(Htheta_pspline1d,Htheta,scoef)
      
    !*** \int(H(theta))dtheta computation ***
    Htheta_int_tmp = sum(scoef(2:Ntheta-2))
    Htheta_int_tmp = 6._RKIND*dtheta*Htheta_int_tmp
    do idtheta = 1,6
      Htheta_int_tmp = Htheta_int_tmp + &
        Htheta_pspline1d%factor_int(idtheta) * &
        scoef(Htheta_pspline1d%indx(idtheta))
    end do
    integral_value = Htheta_int_tmp
  end subroutine compute_omp_theta_integral_CS
      
  !---------------------------------------------------------
  ! Computes for all function H(phi) the integral 
  !  \int(H(phi))dphi where H(phi) is a 
  !  periodic function
  !---------------------------------------------------------
  subroutine compute_omp_phi_integral_CS(Nphi,dphi, &
    Hphi,Hphi_pspline1d,scoef,integral_value)
    use spline1d_class, only : integration_coef_boundaries
    integer                  , intent(in)    :: Nphi
    real(RKIND)              , intent(in)    :: dphi
    real(RKIND), dimension(:), pointer       :: Hphi
    type(pspline1d)          , intent(inout) :: Hphi_pspline1d
    real(RKIND), dimension(:), pointer       :: scoef
    real(RKIND)              , intent(out)   :: integral_value
    
    real(RKIND) :: Hphi_int_tmp
    integer     :: idphi, iphi
      
    !*** case for periodic boundary condition ***
    call period_omp_spline_coef(Hphi_pspline1d,Hphi,scoef)
      
    !*** \int(H(phi))dphi computation ***
    Hphi_int_tmp = sum(scoef(2:Nphi-2))
    Hphi_int_tmp = 6._RKIND*dphi*Hphi_int_tmp
    do idphi = 1,6
      Hphi_int_tmp = Hphi_int_tmp + &
        Hphi_pspline1d%factor_int(idphi) * &
        scoef(Hphi_pspline1d%indx(idphi))
    end do
    integral_value = Hphi_int_tmp
  end subroutine compute_omp_phi_integral_CS
      
  !---------------------------------------------------------
  ! Computes the integral \int(H(vpar))dvpar where 
  !   H(vpar) is a non-periodic function
  !---------------------------------------------------------
  subroutine compute_omp_vpar_integral_CS(Nvpar,dvpar, &
    Hvpar,BCvpar_left,BCvpar_right,Hvpar_nspline1d,&
    scoef,integral_value,Sderiv_Hvpar)
    use spline1d_class, only : integration_coef_boundaries
    integer                      , intent(in)    :: Nvpar
    real(RKIND)                  , intent(in)    :: dvpar
    real(RKIND)    , dimension(:), pointer       :: Hvpar
    integer                      , intent(in)    :: BCvpar_left
    integer                      , intent(in)    :: BCvpar_right
    type(nspline1d)              , intent(in)    :: Hvpar_nspline1d
    real(RKIND)  , dimension(:)  , pointer       :: scoef
    real(RKIND)                  , intent(out)   :: integral_value
    real(RKIND), &
      dimension(0:1) , intent(in), optional :: Sderiv_Hvpar
    
    real(RKIND) :: Hvpar_int_tmp
    integer     :: idvpar, ivpar
      
    !*** case for natural boundary condition ***
    call natural_omp_spline_coef(Hvpar_nspline1d,Hvpar,&
      BCvpar_left,BCvpar_right,scoef,Sderiv_Hvpar)
      
    !*** \int(H(vpar))dvpar computation ***
    Hvpar_int_tmp = sum(scoef(2:Nvpar-2))
    Hvpar_int_tmp = 6._RKIND*dvpar*Hvpar_int_tmp
    do idvpar = 1,6
      Hvpar_int_tmp = Hvpar_int_tmp + &
        Hvpar_nspline1d%factor_int(idvpar) * &
        scoef(Hvpar_nspline1d%indx(idvpar))
    end do
    integral_value = Hvpar_int_tmp
  end subroutine compute_omp_vpar_integral_CS
      
  !--------------------------------------------------
  ! Computes for all function H(r,theta,phi) the integral
  !  I(r) = \int(H(r,theta,phi)dthetadphi
  !--------------------------------------------------
  subroutine compute_omp_intdthetadphi_CS(H,geom, &
    pspline1d_theta,pspline1d_phi,H_intdthetadphi)
    use globals        , only : Rarray_NrNphi
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp_scoeftheta, &
      Romp_scoefphi, Romp1_0Ntheta, Romp1_0Nphi, Romp1_0Nr_0Ntheta
    real(RKIND), &
      dimension(0:,0:,0:), intent(in)    :: H
    type(geometry)       , intent(in)    :: geom
    type(pspline1d)      , intent(inout) :: pspline1d_theta
    type(pspline1d)      , intent(inout) :: pspline1d_phi
    real(RKIND), &
      dimension(0:)      , intent(out)   :: H_intdthetadphi
      
    integer                              :: tid
    integer                              :: ir, itheta, iphi
    real(RKIND), dimension(:,:), pointer :: H_rtheta
    real(RKIND), dimension(:)  , pointer :: rhstheta
    real(RKIND), dimension(:)  , pointer :: scoeftheta
    real(RKIND), dimension(:,:), pointer :: H_rphi
    real(RKIND), dimension(:)  , pointer :: rhsphi
    real(RKIND), dimension(:)  , pointer :: scoefphi
      
    H_rphi => Rarray_NrNphi
      
#ifdef _OPENMP
!$OMP PARALLEL private(ir,itheta,iphi,tid, &
!$OMP rhstheta,scoeftheta,rhsphi,scoefphi,H_rtheta) default(shared)
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
      
    H_rtheta   => Romp1_0Nr_0Ntheta(tid)%val
    rhstheta   => Romp1_0Ntheta(tid)%val
    rhsphi     => Romp1_0Nphi(tid)%val
    scoeftheta => Romp_scoeftheta(tid)%val
    scoefphi   => Romp_scoefphi(tid)%val
      
!$OMP DO SCHEDULE(STATIC) 
    do iphi = 0,geom%Nphi
      do itheta = 0,geom%Ntheta
        do ir = 0,geom%Nr
          H_rtheta(ir,itheta) = H(ir,itheta,iphi)
        end do
      end do
      !*** integral according to theta direction ***
      do ir = 0,geom%Nr
        do itheta = 0,geom%Ntheta
          rhstheta(itheta) = H_rtheta(ir,itheta)
        end do
        call compute_omp_theta_integral_CS(geom%Ntheta, &
          geom%dtheta,rhstheta,pspline1d_theta, &
          scoeftheta,H_rphi(ir,iphi))
      end do
    end do
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC) 
    do ir = 0,geom%Nr
      do iphi = 0,geom%Nphi
        rhsphi(iphi) = H_rphi(ir,iphi)
      end do
      !*** integral according to phi direction ***
      call compute_omp_phi_integral_CS(geom%Nphi,geom%dphi, &
        rhsphi,pspline1d_phi,scoefphi,H_intdthetadphi(ir))
    end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
  end subroutine compute_omp_intdthetadphi_CS
      
  !--------------------------------------------------------
  ! Computes for all function H(r,theta,phi) the integral
  !  \int(H(r,theta,phi) Js dr dtheta dphi
  ! where Js is the jacobian in space, 
  ! by using cubic splines 
  !--------------------------------------------------------
  subroutine compute_omp_intdvolume_CS(H,geom, &
    BCr_left,BCr_right,nspline1d_r,pspline1d_theta, &
    pspline1d_phi,H_intdvolume)
    use globals, only : Rarray1_Nr
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp1_0Ntheta, Romp1_0Nphi, &
      Romp_scoefr, Romp_scoeftheta, Romp_scoefphi 
    use coord_system_class, only : jacobian_space
    real(RKIND), &
      dimension(0:,0:,0:), intent(in)    :: H
    type(geometry)       , intent(in)    :: geom
    integer              , intent(in)    :: BCr_left, BCr_right
    type(nspline1d)      , intent(inout) :: nspline1d_r
    type(pspline1d)      , intent(inout) :: pspline1d_theta
    type(pspline1d)      , intent(inout) :: pspline1d_phi
    real(RKIND)          , intent(out)   :: H_intdvolume
      
    integer :: tid
    integer :: ir, itheta, iphi
    real(RKIND), dimension(:), pointer :: Htmp_theta
    real(RKIND), dimension(:), pointer :: H_intdtheta
    real(RKIND), dimension(:), pointer :: H_intdthetadphi
    real(RKIND), dimension(:), pointer :: scoefr
    real(RKIND), dimension(:), pointer :: scoeftheta
    real(RKIND), dimension(:), pointer :: scoefphi
      
#ifdef _OPENMP
!$OMP PARALLEL private(ir,itheta,iphi,tid, &
!$OMP Htmp_theta,H_intdtheta,H_intdthetadphi, &
!$OMP scoefr,scoeftheta,scoefphi) default(shared)
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    scoefr          => Romp_scoefr(tid)%val
    scoeftheta      => Romp_scoeftheta(tid)%val
    scoefphi        => Romp_scoefphi(tid)%val
    Htmp_theta      => Romp1_0Ntheta(tid)%val
    H_intdtheta     => Romp1_0Nphi(tid)%val
    H_intdthetadphi => Rarray1_Nr
      
!$OMP DO SCHEDULE(STATIC)
    do ir = 0,geom%Nr
      !*** integral according to theta and phi directions ***
      do iphi = 0,geom%Nphi
        do itheta = 0,geom%Ntheta
          Htmp_theta(itheta) = H(ir,itheta,iphi) * &
            jacobian_space(ir,itheta)
        end do
        !-> integral according to theta direction 
        call compute_omp_theta_integral_CS(geom%Ntheta, &
          geom%dtheta,Htmp_theta,pspline1d_theta, &
          scoeftheta,H_intdtheta(iphi))
      end do
      !-> integral according to phi direction
      call compute_omp_phi_integral_CS(geom%Nphi,geom%dphi, &
        H_intdtheta,pspline1d_phi,scoefphi,H_intdthetadphi(ir))
    end do
!$OMP END DO
!$OMP BARRIER
      
!$OMP MASTER
!$OMP END MASTER
    !*** integral according to r direction ***
    call compute_omp_r_integral_CS(geom%Nr, &
      geom%dr,H_intdthetadphi,BCr_left,BCr_right, &
      nspline1d_r,scoefr,H_intdvolume)
!$OMP BARRIER
!$OMP END PARALLEL
  end subroutine compute_omp_intdvolume_CS
      
  !------------------------------------------------------------
  ! Computes the (theta-phi) average of
  !  all function H(r,theta,phi) the integral, i.e
  !  H_avg(r) = 1/(Ltheta*Lphi)*\int(H(r,theta,phi)dthetadphi
  ! by using cubic splines 
  !------------------------------------------------------------
  subroutine compute_omp_thetaphiavg_CS(H,geom, &
    pspline1d_theta,pspline1d_phi,H_avg)
    use globals        , only : Rarray_NrNphi
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp_scoeftheta, &
      Romp_scoefphi, Romp1_0Ntheta, Romp1_0Nphi, Romp1_0Nr_0Ntheta
    real(RKIND), &
      dimension(0:,0:,0:), intent(in)    :: H
    type(geometry)       , intent(in)    :: geom
    type(pspline1d)      , intent(inout) :: pspline1d_theta
    type(pspline1d)      , intent(inout) :: pspline1d_phi
    real(RKIND), &
      dimension(0:)      , intent(out)   :: H_avg
      
    integer                              :: tid
    integer                              :: ir, itheta, iphi
    real(RKIND)                          :: inv_LthetaLphi
    real(RKIND), dimension(:,:), pointer :: H_rtheta
    real(RKIND), dimension(:)  , pointer :: rhstheta
    real(RKIND), dimension(:)  , pointer :: scoeftheta
    real(RKIND), dimension(:,:), pointer :: H_rphi
    real(RKIND), dimension(:)  , pointer :: rhsphi
    real(RKIND), dimension(:)  , pointer :: scoefphi
      
    inv_LthetaLphi = 1._RKIND/(geom%Ltheta*geom%Lphi)
    H_rphi         => Rarray_NrNphi
      
#ifdef _OPENMP
!$OMP PARALLEL private(ir,itheta,iphi,tid, &
!$OMP rhstheta,scoeftheta,rhsphi,scoefphi,H_rtheta) default(shared)
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
      
    H_rtheta   => Romp1_0Nr_0Ntheta(tid)%val
    rhstheta   => Romp1_0Ntheta(tid)%val
    rhsphi     => Romp1_0Nphi(tid)%val
    scoeftheta => Romp_scoeftheta(tid)%val
    scoefphi   => Romp_scoefphi(tid)%val
      
!$OMP DO SCHEDULE(STATIC) 
    do iphi = 0,geom%Nphi
      do itheta = 0,geom%Ntheta
        do ir = 0,geom%Nr
          H_rtheta(ir,itheta) = H(ir,itheta,iphi)
        end do
      end do
      !*** integral according to theta direction ***
      do ir = 0,geom%Nr
        do itheta = 0,geom%Ntheta
          rhstheta(itheta) = H_rtheta(ir,itheta)
        end do
        call compute_omp_theta_integral_CS(geom%Ntheta, &
          geom%dtheta,rhstheta,pspline1d_theta, &
          scoeftheta,H_rphi(ir,iphi))
      end do
    end do
!$OMP END DO
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC) 
    do ir = 0,geom%Nr
      do iphi = 0,geom%Nphi
        rhsphi(iphi) = H_rphi(ir,iphi)
      end do
      !*** integral according to phi direction ***
      call compute_omp_phi_integral_CS(geom%Nphi,geom%dphi, &
        rhsphi,pspline1d_phi,scoefphi,H_avg(ir))
      H_avg(ir) = H_avg(ir)*inv_LthetaLphi
    end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
  end subroutine compute_omp_thetaphiavg_CS
      
  !--------------------------------------------------------------
  ! Compute the flux surface average of a 3D function 
  !  H(r,theta,phi) as :
  !   H_FSavg(r) =  \int H(r,theta,phi) Js(r,theta) dtheta dphi /
  !                 \int Js(r,theta) dtheta dphi
  ! where Js corresponds to the Jacobian in space
  ! by using cubic splines 
  !--------------------------------------------------------------
  subroutine compute_omp_fsurfavg_CS(H,geom, &
    pspline1d_theta,pspline1d_phi,H_Sflux_avg)
      
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp_scoeftheta, &
      Romp_scoefphi, Romp1_0Nr, Romp2_0Nr, &
      Romp1_0Ntheta, Romp2_0Ntheta, &
      Romp1_0Nphi, Romp2_0Nphi
    use coord_system_class, only : jacobian_space
    real(RKIND), &
      dimension(0:,0:,0:), intent(in)    :: H
    type(geometry)       , intent(in)    :: geom
    type(pspline1d)      , intent(inout) :: pspline1d_theta
    type(pspline1d)      , intent(inout) :: pspline1d_phi
    real(RKIND), &
      dimension(0:)      , intent(out)   :: H_Sflux_avg
      
    integer :: tid
    integer :: ir, itheta, iphi
    real(RKIND)                        :: jacob_space_tmp
    real(RKIND), dimension(:), pointer :: H_Js_theta
    real(RKIND), dimension(:), pointer :: Js_theta
    real(RKIND), dimension(:), pointer :: H_Js_intdtheta
    real(RKIND), dimension(:), pointer :: Js_intdtheta
    real(RKIND), dimension(:), pointer :: H_Js_intdthetadphi
    real(RKIND), dimension(:), pointer :: Js_intdthetadphi
    real(RKIND), dimension(:), pointer :: scoeftheta
    real(RKIND), dimension(:), pointer :: scoefphi
      
    H_Js_intdthetadphi => Romp1_0Nr(1)%val
    Js_intdthetadphi   => Romp2_0Nr(1)%val
      
#ifdef _OPENMP
!$OMP PARALLEL private(ir,itheta,iphi,tid, &
!$OMP H_Js_theta,H_Js_intdtheta,Js_theta,Js_intdtheta, &
!$OMP jacob_space_tmp,scoeftheta,scoefphi) default(shared)
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    scoeftheta         => Romp_scoeftheta(tid)%val
    scoefphi           => Romp_scoefphi(tid)%val
    H_Js_theta         => Romp1_0Ntheta(tid)%val
    Js_theta           => Romp2_0Ntheta(tid)%val
    H_Js_intdtheta     => Romp1_0Nphi(tid)%val
    Js_intdtheta       => Romp2_0Nphi(tid)%val
      
!$OMP DO SCHEDULE(STATIC)
    do ir = 0,geom%Nr
      !*** integral according to theta and phi directions ***
      do iphi = 0,geom%Nphi
        do itheta = 0,geom%Ntheta
          jacob_space_tmp    = jacobian_space(ir,itheta)
          H_Js_theta(itheta) = H(ir,itheta,iphi) * &
            jacob_space_tmp
          Js_theta(itheta)   = jacob_space_tmp
        end do
        !-> compute \int H Js dtheta
        call compute_omp_theta_integral_CS(geom%Ntheta, &
          geom%dtheta,H_Js_theta,pspline1d_theta, &
          scoeftheta,H_Js_intdtheta(iphi))
        !-> compute \int Js dtheta
        call compute_omp_theta_integral_CS(geom%Ntheta, &
          geom%dtheta,Js_theta,pspline1d_theta, &
          scoeftheta,Js_intdtheta(iphi))
      end do
      !-> compute \int H Js dtheta dphi
      call compute_omp_phi_integral_CS(geom%Nphi,geom%dphi, &
        H_Js_intdtheta,pspline1d_phi,scoefphi, &
        H_Js_intdthetadphi(ir))
      !-> compute \int Js dtheta dphi
      call compute_omp_phi_integral_CS(geom%Nphi,geom%dphi, &
        Js_intdtheta,pspline1d_phi,scoefphi,Js_intdthetadphi(ir))
      !-> compute \int H Js dtheta dphi / \int Js dtheta dphi 
      H_Sflux_avg(ir) = H_Js_intdthetadphi(ir)/Js_intdthetadphi(ir)
    end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
  end subroutine compute_omp_fsurfavg_CS
      
  !***********************************************
  !  INTEGRATION WITH COLLOCATION METHOD
  !***********************************************
  !---------------------------------------------------------
  ! Computes for all function H(x) the integral \int(H(x))dx
  ! with a collocation method defined by the collocation 
  ! coefficient "coeff_intdx"
  !---------------------------------------------------------
  subroutine compute_phase_integral_colloc(Nx,Hx,&
    coeff_intdx,intH_dx)
    integer                     , intent(in)  :: Nx
    real(RKIND), dimension(0:Nx), intent(in)  :: Hx
    real(RKIND), dimension(0:Nx), intent(in)  :: coeff_intdx
    real(RKIND)                 , intent(out) :: intH_dx
    
    integer :: ix
      
    intH_dx = 0._RKIND
    do ix = 0,Nx
      intH_dx = intH_dx + coeff_intdx(ix)*Hx(ix)
    end do
  end subroutine compute_phase_integral_colloc
      
  !--------------------------------------------------
  ! Computes for all function H(r) the integral
  !  I(r) = \int H(r)rdr
  !  with a collocation method
  !--------------------------------------------------
  subroutine compute_intrdr_colloc(Hr,geom,Hr_intrdr)
    real(RKIND), dimension(0:), intent(in)  :: Hr
    type(geometry)            , intent(in)  :: geom
    real(RKIND)               , intent(out) :: Hr_intrdr
      
    integer     :: ir
      
    Hr_intrdr = 0._RKIND
    do ir = 0,geom%Nr
      Hr_intrdr = Hr_intrdr + &
        geom%coeff_intdr(ir)*Hr(ir)*geom%rg(ir)
    end do
  end subroutine compute_intrdr_colloc
      
  !--------------------------------------------------
  ! Computes for all function H(theta) the integral
  !  I = \int H(theta) dtheta with a collocation method
  !--------------------------------------------------
  subroutine compute_theta_integral_colloc(Htheta, &
    geom,Htheta_intdtheta)
    real(RKIND), dimension(0:), intent(in)  :: Htheta
    type(geometry)            , intent(in)  :: geom
    real(RKIND)               , intent(out) :: Htheta_intdtheta
      
    integer :: itheta
      
    Htheta_intdtheta = 0._RKIND
    do itheta = 0,geom%Ntheta
      Htheta_intdtheta = Htheta_intdtheta + &
        geom%coeff_intdtheta(itheta)*Htheta(itheta)
    end do
  end subroutine compute_theta_integral_colloc
      
  !--------------------------------------------------
  ! Computes for all function H(phi) the integral
  !  I = \int H(phi) dphi with a collocation method
  !--------------------------------------------------
  subroutine compute_phi_integral_colloc(Hphi,geom,Hphi_intdphi)
    real(RKIND), dimension(0:), intent(in)  :: Hphi
    type(geometry)            , intent(in)  :: geom
    real(RKIND)               , intent(out) :: Hphi_intdphi
      
    integer :: iphi
      
    Hphi_intdphi = 0._RKIND
    do iphi = 0,geom%Nphi
      Hphi_intdphi = Hphi_intdphi + &
        geom%coeff_intdphi(iphi)*Hphi(iphi)
    end do
  end subroutine compute_phi_integral_colloc
      
  !--------------------------------------------------
  ! Computes for all function H(vpar) the integral
  !  I = \int H(vpar) dvpar with a collocation method
  !--------------------------------------------------
  subroutine compute_vpar_integral_colloc(Hvpar,geom,Hvpar_intdvpar)
    real(RKIND), dimension(0:), intent(in)  :: Hvpar
    type(geometry)            , intent(in)  :: geom
    real(RKIND)               , intent(out) :: Hvpar_intdvpar
      
    integer :: ivpar
      
    Hvpar_intdvpar = 0._RKIND
    do ivpar = 0,geom%Nvpar
      Hvpar_intdvpar = Hvpar_intdvpar + &
        geom%coeff_intdvpar(ivpar)*Hvpar(ivpar)
    end do
  end subroutine compute_vpar_integral_colloc
      
  !--------------------------------------------------
  ! Computes for all function H(r,theta,phi) the integral
  !  I(r) = \int H(r,theta,phi)dthetadphi 
  !  with a collocation method
  !--------------------------------------------------
  subroutine compute_intdthetadphi_colloc(H,geom,H_intdthetadphi)
    real(RKIND), dimension(0:,0:,0:), intent(in)  :: H
    type(geometry)                  , intent(in)  :: geom
    real(RKIND), dimension(0:)      , intent(out) :: H_intdthetadphi
      
    integer     :: ir, itheta, iphi
    real(RKIND) :: H_intdtheta, H_intdthetadphi_tmp
      
    do ir = 0,geom%Nr
      !*** integral according to phi direction ***
      H_intdthetadphi_tmp = 0._RKIND
      do iphi = 0,geom%Nphi
        !*** integral according to theta direction ***
        H_intdtheta = 0._RKIND
        do itheta = 0,geom%Ntheta
          H_intdtheta = H_intdtheta + &
            geom%coeff_intdtheta(itheta)*H(ir,itheta,iphi)
        end do
        H_intdthetadphi_tmp = H_intdthetadphi_tmp + &
          geom%coeff_intdphi(iphi)*H_intdtheta
      end do
      H_intdthetadphi(ir) = H_intdthetadphi_tmp
    end do
  end subroutine compute_intdthetadphi_colloc
      
  !--------------------------------------------------------
  ! Computes for all function H(r,theta,phi) the integral
  !  \int H(r,theta,phi) Js dr dtheta dphi
  ! where Js is the jacobian in space, with a
  ! collocation method
  !--------------------------------------------------------
  subroutine compute_intdvolume_colloc(H,geom,H_intdvolume)
    use coord_system_class, only : jacobian_space
    real(RKIND), &
      dimension(0:,0:,0:), intent(in)  :: H
    type(geometry)       , intent(in)  :: geom
    real(RKIND)          , intent(out) :: H_intdvolume
      
    integer     :: ir, itheta, iphi
    real(RKIND) :: coeff_intdr, coeff_intdtheta, coeff_intdphi
    real(RKIND) :: H_intdthetadphi_tmp
      
    H_intdvolume = 0._RKIND
    do ir = 0,geom%Nr
      coeff_intdr         = geom%coeff_intdr(ir)
      H_intdthetadphi_tmp = 0._RKIND
      !*** integrals in theta and phi direction ***
      do iphi = 0,geom%Nphi
        coeff_intdphi = geom%coeff_intdphi(iphi)
        do itheta = 0,geom%Ntheta
          coeff_intdtheta     = geom%coeff_intdtheta(itheta)
          H_intdthetadphi_tmp = H_intdthetadphi_tmp + &
            H(ir,itheta,iphi)*jacobian_space(ir,itheta) * &
            coeff_intdphi*coeff_intdtheta
        end do
      end do
      !*** integral in r direction ***
      H_intdvolume = H_intdvolume + &
        H_intdthetadphi_tmp*coeff_intdr
    end do
  end subroutine compute_intdvolume_colloc
      
  !---------------------------------------------------------------- 
  ! Compute the (theta,phi)-average of a 3D function 
  !  H(r,theta,phi) as :
  !   H_avg(r) = 1/(Ltheta*Lphi) \int H(r,theta,phi) dtheta dphi 
  ! with a collocation method
  !---------------------------------------------------------------- 
  subroutine compute_thetaphiavg_colloc(H,geom,H_avg)
    use geometry_class
    real(RKIND), dimension(0:,0:,0:), intent(in)    :: H
    type(geometry)                  , intent(in)    :: geom
    real(RKIND), dimension(0:)      , intent(inout) :: H_avg
      
    integer     :: ir, itheta, iphi
    real(RKIND) :: inv_LthetaLphi
    real(RKIND) :: coeff_intdtheta, coeff_intdphi
      
    inv_LthetaLphi = 1._RKIND/(geom%Ltheta*geom%Lphi)
    do ir = 0,geom%Nr
      H_avg(ir) = 0._RKIND
    end do
    do iphi = 0,geom%Nphi
      coeff_intdphi = geom%coeff_intdphi(iphi)
      do itheta = 0,geom%Ntheta
        coeff_intdtheta = geom%coeff_intdtheta(itheta)
        do ir = 0,geom%Nr
          H_avg(ir) =  H_avg(ir) + &
            H(ir,itheta,iphi)*coeff_intdtheta*coeff_intdphi
        end do
      end do
    end do
    do ir = 0,geom%Nr
      H_avg(ir) = H_avg(ir)*inv_LthetaLphi
    end do
  end subroutine compute_thetaphiavg_colloc
      
  !---------------------------------------------------------------- 
  ! Compute the flux surface average of a 3D function 
  !  H(r,theta,phi) as :
  !   H_FSavg(r) =  \int H(r,theta,phi) Js(r,theta) dtheta dphi /
  !                 \int Js(r,theta) dtheta dphi
  ! where Js corresponds to the Jacobian in space, 
  ! with a collocation method
  !---------------------------------------------------------------- 
  subroutine compute_fluxsurfavg_colloc(H,geom,H_Sflux_avg)
    use geometry_class
    use coord_system_class, only : jacobian_space, intdthetadphi_Js 
    real(RKIND), &
      dimension(0:,0:,0:), intent(in)    :: H
    type(geometry)       , intent(in)    :: geom
    real(RKIND), &
      dimension(0:)      , intent(inout) :: H_Sflux_avg
      
    integer     :: ir, itheta, iphi
    real(RKIND) :: intdS_H_tmp, intdthetadphi_Js_tmp
    real(RKIND) :: coeff_intdtheta, coeff_intdphi
    real(RKIND) :: jacob_space_ij
      
    do ir = 0,geom%Nr
      intdS_H_tmp          = 0._RKIND
      intdthetadphi_Js_tmp = intdthetadphi_Js(ir)
      do iphi = 0,geom%Nphi
        coeff_intdphi =  geom%coeff_intdphi(iphi)
        do itheta = 0,geom%Ntheta
          jacob_space_ij  = jacobian_space(ir,itheta)
          coeff_intdtheta = geom%coeff_intdtheta(itheta)
          intdS_H_tmp     = intdS_H_tmp + &
            h(ir,itheta,iphi)*jacob_space_ij * &
            coeff_intdphi*coeff_intdtheta
        end do
      end do
      H_Sflux_avg(ir) = intdS_H_tmp/intdthetadphi_Js_tmp
    end do
  end subroutine compute_fluxsurfavg_colloc
      
  !***************************************************
  !  CHOICE BETWEEN INTEGRATION WITH CUBIC SPLINES
  !  OR COLLOCATION METHOD
  !***************************************************
  !----------------------------------------------------------
  ! Computes for all function H(r,theta,phi) the integral
  !  I(r) = \int(H(r,theta,phi)dthetadphi
  !
  ! by choosing between cubic splines or collocation method
  !----------------------------------------------------------
  subroutine compute_intdthetadphi(H,geom, &
    pspline1d_theta,pspline1d_phi,H_intdthetadphi)
    use globals, only : integration_CS
    real(RKIND), &
      dimension(0:,0:,0:)    , intent(in)    :: H
    type(geometry)           , intent(in)    :: geom
    type(pspline1d), optional, intent(inout) :: pspline1d_theta
    type(pspline1d), optional, intent(inout) :: pspline1d_phi
    real(RKIND), &
      dimension(0:)          , intent(out)   :: H_intdthetadphi
    
    if (integration_CS) then
      call compute_omp_intdthetadphi_CS(H,geom, &
        pspline1d_theta,pspline1d_phi,H_intdthetadphi)
    else
      call compute_intdthetadphi_colloc(H,geom,H_intdthetadphi)
    end if
  end subroutine compute_intdthetadphi
      
  !--------------------------------------------------------
  ! Computes for all function H(r,theta,phi) the integral
  !  \int(H(r,theta,phi) Js dr dtheta dphi
  ! where Js is the jacobian in space,
  !
  ! by choosing between cubic splines or collocation method
  !--------------------------------------------------------
  subroutine compute_intdvolume(H,geom,BCr_left,BCr_right, &
    nspline1d_r,pspline1d_theta,pspline1d_phi,H_intrdrdthdphi)
    use globals, only : integration_CS
    real(RKIND)    , &
          dimension(0:,0:,0:), intent(in)    :: H
    type(geometry)           , intent(in)    :: geom
    integer                  , intent(in)    :: BCr_left, BCr_right
    type(nspline1d), optional, intent(inout) :: nspline1d_r
    type(pspline1d), optional, intent(inout) :: pspline1d_theta
    type(pspline1d), optional, intent(inout) :: pspline1d_phi
    real(RKIND)              , intent(out)   :: H_intrdrdthdphi
      
    if (integration_CS) then
      call compute_omp_intdvolume_CS(H,geom, &
        BCr_left,BCr_right,nspline1d_r,pspline1d_theta, &
        pspline1d_phi,H_intrdrdthdphi)
    else
      call compute_intdvolume_colloc(H,geom,H_intrdrdthdphi)
    end if
  end subroutine compute_intdvolume
      
  !---------------------------------------------------------------- 
  ! Compute the (theta,phi)-average of a 3D function 
  !  H(r,theta,phi) as :
  !   H_avg(r) = 1/(Ltheta*Lphi) \int H(r,theta,phi) dtheta dphi 
  !
  ! by choosing between cubic splines or collocation method
  !---------------------------------------------------------------- 
  subroutine compute_thetaphi_average(H,geom, &
    pspline1d_theta,pspline1d_phi,H_avg)
    use globals, only : integration_CS
    real(RKIND), &
      dimension(0:,0:,0:)    , intent(in)    :: H
    type(geometry)           , intent(in)    :: geom
    type(pspline1d), optional, intent(inout) :: pspline1d_theta
    type(pspline1d), optional, intent(inout) :: pspline1d_phi
    real(RKIND), &
      dimension(0:)          , intent(out)   :: H_avg
    
    if ((integration_CS).and.present(pspline1d_theta)) then
      call compute_omp_thetaphiavg_CS(H,geom, &
        pspline1d_theta,pspline1d_phi,H_avg)
    else
      call compute_thetaphiavg_colloc(H,geom,H_avg)
    end if
  end subroutine compute_thetaphi_average
      
  !---------------------------------------------------------------- 
  ! Compute the flux surface average of a 3D function 
  !  H(r,theta,phi) as :
  !   H_FSavg(r) =  \int H(r,theta,phi) Js(r,theta) dtheta dphi /
  !                 \int Js(r,theta) dtheta dphi
  ! where Js corresponds to the Jacobian in space
  !
  ! by choosing between cubic splines or collocation method
  !---------------------------------------------------------------- 
  subroutine compute_fluxsurf_average(H,geom, &
    pspline1d_theta,pspline1d_phi,H_Sflux_avg)
    use globals, only : integration_CS
    real(RKIND), &
      dimension(0:,0:,0:)    , intent(in)    :: H
    type(geometry)           , intent(in)    :: geom
    type(pspline1d), optional, intent(inout) :: pspline1d_theta
    type(pspline1d), optional, intent(inout) :: pspline1d_phi
    real(RKIND), &
      dimension(0:)          , intent(out)   :: H_Sflux_avg
    
    if ((integration_CS).and.present(pspline1d_theta)) then
      call compute_omp_fsurfavg_CS(H,geom, &
        pspline1d_theta,pspline1d_phi,H_Sflux_avg)
    else
      call compute_fluxsurfavg_colloc(H,geom,H_Sflux_avg)
    end if
  end subroutine compute_fluxsurf_average
end module integration_module
