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
      
!---------------------------------------------
! file : efield.f90
! date : 26/06/2000
!  computation of the electric potential and 
!  the three components of the electric field
!  which is given by the quasi-neutrality
!  equation
! (all the quantities are computed using the
!  cubic spline expression)
!---------------------------------------------
module efield_module
  use prec_const
  use globals, only : Zi, Nr, Ntheta, Nphi
  use mem_alloc_module
  use geometry_class
  use spline1d_class
  implicit none
      
  public :: compute_dPhidr_loc, compute_dPhidtheta_loc, &
    compute_dPhidphi_loc
      
  type(nspline1d), private :: E_nspline1d_r
  type(pspline1d), private :: E_pspline1d_theta
  type(pspline1d), private :: E_pspline1d_phi
     
  !******************************
  contains
  !******************************
  
  
  !-------------------------------------------------
  ! Constructor 
  !-------------------------------------------------
  subroutine new_efield(geom)
    use globals, only : dJ0Phidr, dJ0Phidtheta, dJ0Phidphi
    type(geometry), intent(in) :: geom
      
    call new_spline1d_natural(E_nspline1d_r,geom%Nr,geom%dr)
    call new_spline1d_period(E_pspline1d_theta, &
      geom%Ntheta,geom%dtheta)
    call new_spline1d_period(E_pspline1d_phi,geom%Nphi,geom%dphi)
      
    !*** allocation for the partial derivatives of Phi ***
    call glob_allocate(dJ0Phidr,0,geom%Nr, &
      0,geom%Ntheta,0,geom%Nphi,'dJ0Phidr')
    call glob_allocate(dJ0Phidtheta,0,geom%Nr, &
      0,geom%Ntheta,0,geom%Nphi,'dJ0Phidtheta')
    call glob_allocate(dJ0Phidphi,0,geom%Nr, &
      0,geom%Ntheta,0,geom%Nphi,'dJ0Phidphi')
  end subroutine new_efield
      
  !-------------------------------------------------
  ! Destructor 
  !-------------------------------------------------
  subroutine delete_efield
    use globals, only : dJ0Phidr, dJ0Phidtheta, dJ0Phidphi
    call del_spline1d_natural(E_nspline1d_r)
    call del_spline1d_period(E_pspline1d_theta)
    call del_spline1d_period(E_pspline1d_phi)
    !*** deallocation of the partial derivatives of Phi ***
    call glob_deallocate(dJ0Phidr)
    call glob_deallocate(dJ0Phidtheta)
    call glob_deallocate(dJ0Phidphi)
  end subroutine delete_efield
      
  !-------------------------------------------------
  ! Computation of the radial derivate of the electric
  ! potential d(Phi)/dr by using the spline 
  ! interpolation coefficients 
  !-------------------------------------------------    
  subroutine compute_dPhidr_loc(geom,Phi, &
    begin_dim1,end_dim1,begin_dim2,end_dim2, &
    begin_dim3,end_dim3,dPhidr)
    use globals, only : Phi_BCr_left,Phi_BCr_right
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp1_0Nr, Romp2_0Nr, &
      Romp_scoefr
    type(geometry)       , intent(in) :: geom
    real(RKIND), &
      dimension(0:,0:,0:), intent(in) :: Phi
    integer              , intent(in) :: begin_dim1
    integer              , intent(in) :: end_dim1
    integer              , intent(in) :: begin_dim2
    integer              , intent(in) :: end_dim2
    integer              , intent(in) :: begin_dim3
    integer              , intent(in) :: end_dim3
    real(RKIND), &
      dimension(:,:,:)   , pointer    :: dPhidr
    
    integer                            :: tid, ir, itheta, iphi
    real(RKIND), dimension(:), pointer :: scoefr, rhsr, drhs_dr
      
    !*** dPhi/dr computation ***  
#ifdef _OPENMP
!$OMP PARALLEL private(ir,itheta,iphi,&
!$OMP rhsr,scoefr,drhs_dr,tid) default(shared)
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    rhsr    => Romp1_0Nr(tid)%val
    scoefr  => Romp_scoefr(tid)%val
    drhs_dr => Romp2_0Nr(tid)%val
      
!$OMP DO SCHEDULE(static)
    do iphi = begin_dim3,end_dim3
      do itheta = begin_dim2,end_dim2
        do ir = 0,geom%Nr
          rhsr(ir) = Phi(ir,itheta,iphi)
        end do
        call natural_omp_deriv1d(E_nspline1d_r,geom%Nr,geom%dr, &
          Phi_BCr_left,Phi_BCr_right,rhsr,scoefr,drhs_dr)
        do ir = begin_dim1,end_dim1
          dPhidr(ir,itheta,iphi) = drhs_dr(ir) 
        end do
      enddo
    end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
  end subroutine compute_dPhidr_loc
    
      
  !-------------------------------------------------
  ! Computation of the derivate in theta of the 
  ! electric potential d(Phi)/dtheta by using  
  ! the spline interpolation coefficients
  !-------------------------------------------------    
  subroutine compute_dPhidtheta_loc(geom,Phi, &
    begin_dim1,end_dim1,begin_dim2,end_dim2, &
    begin_dim3,end_dim3,dPhidth)
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp1_0Ntheta, Romp_scoeftheta
    type(geometry)       , intent(in) :: geom
    real(RKIND), &
      dimension(0:,0:,0:), intent(in) :: Phi
    integer              , intent(in) :: begin_dim1
    integer              , intent(in) :: end_dim1
    integer              , intent(in) :: begin_dim2
    integer              , intent(in) :: end_dim2
    integer              , intent(in) :: begin_dim3
    integer              , intent(in) :: end_dim3
    real(RKIND), &
         dimension(:,:,:), pointer    :: dPhidth
    
    integer                            :: ir, itheta, iphi, tid
    real(RKIND), dimension(:), pointer :: scoeftheta, rhstheta
      
    ! *** dPhi/dtheta computation ***  
#ifdef _OPENMP
!$OMP PARALLEL private(ir,itheta,iphi,&
!$OMP rhstheta,scoeftheta,tid) default(shared)
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    scoeftheta => Romp_scoeftheta(tid)%val
    rhstheta   => Romp1_0Ntheta(tid)%val
      
!$OMP DO SCHEDULE(static)
    do iphi = begin_dim3,end_dim3
      do ir = begin_dim1,end_dim1
        do itheta = 0,geom%Ntheta
          rhstheta(itheta) = Phi(ir,itheta,iphi)
        end do
        call period_omp_deriv1d(E_pspline1d_theta, &
          geom%Ntheta,geom%dtheta,rhstheta,scoeftheta)
      
        do itheta = begin_dim2,end_dim2
          dPhidth(ir,itheta,iphi) = rhstheta(itheta) 
        end do
      enddo
    end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
  end subroutine compute_dPhidtheta_loc
    
      
  !-------------------------------------------------
  ! Computation of the derivate in phi of the 
  ! electric potential d(Phi)/dphi by using  
  ! the spline interpolation coefficients
  !-------------------------------------------------    
  subroutine compute_dPhidphi_loc(geom,Phi, &
    begin_dim1,end_dim1,begin_dim2,end_dim2, &
    begin_dim3,end_dim3,dPhidphi)
    use globals        , only : Rarray1_Nphi, Rarray2_Nphi
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    type(geometry)       , intent(in) :: geom
    real(RKIND), &
      dimension(0:,0:,0:), intent(in) :: Phi
    integer              , intent(in) :: begin_dim1
    integer              , intent(in) :: end_dim1
    integer              , intent(in) :: begin_dim2
    integer              , intent(in) :: end_dim2
    integer              , intent(in) :: begin_dim3
    integer              , intent(in) :: end_dim3
    real(RKIND), &
         dimension(:,:,:), pointer    :: dPhidphi
    
    integer                            :: ir, itheta, iphi
    real(RKIND), dimension(:), pointer :: Phiphi_tmp, dPhidphi_tmp
      
    Phiphi_tmp   => Rarray1_Nphi
    dPhidphi_tmp => Rarray2_Nphi
      
    ! *** dPhi/dphi computation ***  
    do itheta = begin_dim2,end_dim2
      do ir = begin_dim1,end_dim1
        do iphi = 0,geom%Nphi
          Phiphi_tmp(iphi) = Phi(ir,itheta,iphi)
        end do
        call period_deriv1d(E_pspline1d_phi,geom%Nphi,geom%dphi, &
          Phiphi_tmp,dPhidphi_tmp)
        do iphi = begin_dim3,end_dim3
          dPhidphi(ir,itheta,iphi) = dPhidphi_tmp(iphi) 
        end do
      enddo
    end do
  end subroutine compute_dPhidphi_loc
      
  !----------------------------------------------------------
  ! Computation of the three partial derivatives of Phi, i.e
  !  dPhi/dr, dPhi/dtheta and dPhi/dphi
  !  (used for the BSL method)
  !----------------------------------------------------------
  subroutine compute_dPhidx(geom,sPhi, &
    SdPhidr,SdPhidtheta, SdPhidphi)
    use globals        , only : Phi_BCr_left, Phi_BCr_right
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp1_0Nr, Romp2_0Nr, &
      Romp1_0Ntheta, Romp1_0Nphi, &
      Romp_scoefr, Romp_scoeftheta, Romp_scoefphi 
    use geometry_class
    type(geometry)       , intent(in)  :: geom
    real(RKIND), &
      dimension(0:,0:,0:), intent(in)  :: sPhi
    real(RKIND), &
      dimension(0:,0:,0:), intent(out) :: SdPhidr
    real(RKIND), &
      dimension(0:,0:,0:), intent(out) :: SdPhidtheta
    real(RKIND), &
      dimension(0:,0:,0:), intent(out) :: SdPhidphi
      
    integer :: tid
    integer :: ir, itheta, iphi
      
    real(RKIND), dimension(:), pointer :: rhsr
    real(RKIND), dimension(:), pointer :: drhs_dr
    real(RKIND), dimension(:), pointer :: rhstheta
    real(RKIND), dimension(:), pointer :: rhsphi
    real(RKIND), dimension(:), pointer :: scoefr
    real(RKIND), dimension(:), pointer :: scoeftheta
    real(RKIND), dimension(:), pointer :: scoefphi
      
#ifdef _OPENMP
!$OMP PARALLEL private(ir,itheta,iphi, &
!$OMP rhsr,drhs_dr,rhstheta,rhsphi, &
!$OMP scoefr,scoeftheta,scoefphi,tid) default(shared)
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    rhsr       => Romp1_0Nr(tid)%val
    drhs_dr    => Romp2_0Nr(tid)%val
    rhstheta   => Romp1_0Ntheta(tid)%val
    rhsphi     => Romp1_0Nphi(tid)%val
    scoefr     => Romp_scoefr(tid)%val
    scoeftheta => Romp_scoeftheta(tid)%val
    scoefphi   => Romp_scoefphi(tid)%val
      
    !*** computation of dPhi/dr ***
!$OMP DO SCHEDULE(static)
    do iphi = 0,geom%Nphi
      do itheta = 0,geom%Ntheta
        do ir = 0,geom%Nr
          rhsr(ir) = sPhi(ir,itheta,iphi)
        end do
        !-> computation of dPhi/dr (saved in rhsr)
        call natural_omp_deriv1d(E_nspline1d_r, &
          geom%Nr,geom%dr,Phi_BCr_left,Phi_BCr_right, &
          rhsr,scoefr,drhs_dr)
        do ir = 0,geom%Nr
          SdPhidr(ir,itheta,iphi) = drhs_dr(ir)
        end do
      enddo
    end do
!$OMP END DO
!$OMP BARRIER
      
    !*** computation of dPhi/dtheta ***
!$OMP DO SCHEDULE(static)
    do iphi = 0,geom%Nphi
      do ir = 0,geom%Nr
        do itheta = 0,geom%Ntheta
          rhstheta(itheta) = sPhi(ir,itheta,iphi)
        end do
        !-> computation of dPhi/dtheta (saved in rhstheta)
        call period_omp_deriv1d(E_pspline1d_theta, &
          geom%Ntheta,geom%dtheta,rhstheta,scoeftheta)
        do itheta = 0,geom%Ntheta
          SdPhidtheta(ir,itheta,iphi) = rhstheta(itheta)
        end do
      enddo
    end do
!$OMP END DO
!$OMP BARRIER
      
!$OMP MASTER
    !*** computation of dPhi/dphi ***
    do itheta = 0,geom%Ntheta
      do ir = 0,geom%Nr
        do iphi = 0,geom%Nphi
          rhsphi(iphi) = sPhi(ir,itheta,iphi)
        end do
        !-> computation of dPhi/dphi (saved in rhsphi)
        call period_omp_deriv1d(E_pspline1d_phi, &
          geom%Nphi,geom%dphi,rhsphi,scoefphi)
        do iphi = 0,geom%Nphi
          SdPhidphi(ir,itheta,iphi) = rhsphi(iphi)
        end do
      enddo
    end do
!$OMP END MASTER
!$OMP BARRIER
!$OMP END PARALLEL
  end subroutine compute_dPhidx
end module efield_module
