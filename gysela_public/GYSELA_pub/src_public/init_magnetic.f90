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
! file : init_magnetic.f90
! date : 21/06/2000
!  Initialization of the magnetic configuration
!   vec_B = B0 R0/R(r,theta)*[zeta(r)*e_theta+e_phi]
!  with zeta(r) = r/(q R0)
!  and R(r,theta) = R0+h(r)*cos(theta) 
!  (Rk : h(r) = r but should be a profile vanishing
!    at the boundaries) 
!-------------------------------------------------------
module init_magnetic_class
  use prec_const
  use mem_alloc_module
  use geometry_class
  use coord_system_class, only : inv_R
  use init_current_class
      
  implicit none
      
  !*** init_magnet type definition ***
  type :: init_magnetic
    !-> covariant coordinates of B
    real(RKIND), dimension(:,:), pointer :: Br
    real(RKIND), dimension(:,:), pointer :: Btheta
    real(RKIND), dimension(:,:), pointer :: Bphi
    !-> contravariant coordinates of B
    real(RKIND), dimension(:,:), pointer :: B_gradr
    real(RKIND), dimension(:,:), pointer :: B_gradtheta
    real(RKIND), dimension(:,:), pointer :: B_gradphi
    !-> norm of B
    real(RKIND), dimension(:,:), pointer :: B_norm
    real(RKIND), dimension(:)  , pointer :: zeta
    real(RKIND), dimension(:)  , pointer :: b0
    !-> Bmin and Bmax
    real(RKIND)                          :: Bmin
    real(RKIND)                          :: Bmax
    !-> partial derivatives of B(r,theta)
    real(RKIND), dimension(:,:), pointer :: dBdr
    real(RKIND), dimension(:,:), pointer :: dBdtheta
  end type init_magnetic
      
  include "Bstar_inline.h"
  !******************************
  contains
  !******************************
#include "Bstar_inline.f90"
     
  !---------------------------------------------------- 
  ! Constructor for magnetic field configuration 
  !----------------------------------------------------
  subroutine new_magnet_config(mthis,geom)
    use globals, only : Bstar_3D, bstar_gradr_3D, &
      bstar_gradtheta_3D, bstar_gradphi_3D
    type(init_magnetic), intent(out) :: mthis
    type(geometry)     , intent(in)  :: geom
    
    !-> allocation for B covariant coordinates and B norm
    call glob_allocate(mthis%B_norm,0,geom%Nr,0, &
      geom%Ntheta,'B_norm')
    call glob_allocate(mthis%Br,0,geom%Nr,0,geom%Ntheta,'Br')
    call glob_allocate(mthis%Btheta,0,geom%Nr, &
      0,geom%Ntheta,'Btheta')
    call glob_allocate(mthis%Bphi,0,geom%Nr, &
      0,geom%Ntheta,'Bphi')
    call glob_allocate(mthis%zeta,0,geom%Nr,'zeta')
    call glob_allocate(mthis%b0,0,geom%Nr,'b0')
    !-> allocation for the partial derivatives of B
    call glob_allocate(mthis%dBdr,0,geom%Nr, &
      0,geom%Ntheta,'dBdr')
    call glob_allocate(mthis%dBdtheta,0,geom%Nr, &
      0,geom%Ntheta,'dBdtheta')
    !-> allocation for the contravariant components of B
    call glob_allocate(mthis%B_gradr,0,geom%Nr, &
      0,geom%Ntheta,'B_gradr')
    call glob_allocate(mthis%B_gradtheta,0,geom%Nr, &
      0,geom%Ntheta,'B_gradtheta')
    call glob_allocate(mthis%B_gradphi,0,geom%Nr, &
      0,geom%Ntheta,'B_gradphi')
    !-> allocation for the computation of Bstar and
    !   the contravariant components of bstar (bstar=Bstar/B)
    call glob_allocate(Bstar_3D,0,geom%Nr,0,geom%Ntheta,&
      0,geom%Nvpar,'Bstar_3D')
    call glob_allocate(bstar_gradr_3D,0,geom%Nr,0,geom%Ntheta,&
      0,geom%Nvpar,'bstar_gradr_3D')
    call glob_allocate(bstar_gradtheta_3D,0,geom%Nr,0,geom%Ntheta,&
      0,geom%Nvpar,'bstar_gradtheta_3D')
    call glob_allocate(bstar_gradphi_3D,0,geom%Nr,0,geom%Ntheta,&
      0,geom%Nvpar,'bstar_gradphi_3D')
  end subroutine new_magnet_config
  
      
  !---------------------------------------------------- 
  ! destructor for magnetic field configuration 
  !---------------------------------------------------- 
  subroutine del_magnet_config(mthis)
    use globals, only : Bstar_3D, bstar_gradr_3D, &
      bstar_gradtheta_3D, bstar_gradphi_3D
    type(init_magnetic), intent(inout) :: mthis
      
    !-> deallocation of B covariant coordinates and B norm
    call glob_deallocate(mthis%B_norm)
    call glob_deallocate(mthis%Br)
    call glob_deallocate(mthis%Btheta)
    call glob_deallocate(mthis%Bphi)
    call glob_deallocate(mthis%zeta)
    call glob_deallocate(mthis%b0)
    !-> deallocation of the partial derivatives of B
    call glob_deallocate(mthis%dBdr)
    call glob_deallocate(mthis%dBdtheta)
    !-> deallocation of the contravariant components of B
    call glob_deallocate(mthis%B_gradr)
    call glob_deallocate(mthis%B_gradtheta)
    call glob_deallocate(mthis%B_gradphi)
    !-> deallocation of the computation of Bstar and
    !   the contravariant components of bstar (bstar=Bstar/B)
    call glob_deallocate(Bstar_3D)
    call glob_deallocate(bstar_gradr_3D)
    call glob_deallocate(bstar_gradtheta_3D)
    call glob_deallocate(bstar_gradphi_3D)
  end subroutine del_magnet_config
      
  !--------------------------------------------------------
  ! - initialization of zeta(r) = r/(q R0)
  !--------------------------------------------------------
  subroutine init_zeta(mthis,geom,init_prof)
    use globals, only : R0
    use init_profile_class
    type(init_magnetic), intent(inout) :: mthis  
    type(geometry)     , intent(in)    :: geom
    type(init_profile) , intent(in)    :: init_prof
      
    integer     :: ir
    real(RKIND) :: zeta_ri
      
    !*** initialisation of zeta(r) = r/(q R0) ***
    do ir = 0,geom%Nr
      zeta_ri        = geom%rg(ir)*init_prof%iota(ir)/R0
      mthis%zeta(ir) = zeta_ri
    end do
  end subroutine init_zeta
      
  !---------------------------------------------------- 
  ! - initialization of the covariant coordinates of
  ! the magnetic field, i.e:
  !   . Br(r,theta)     = 0
  !   . Btheta(r,theta) = B0*R0*r/R*zeta(r)
  !   . Bphi(r,theta)   = B0*R0
  ! because the magnetic field is defined as:
  !   . vec_B = B0R0/R*(zeta*vec_etheta+vec_ephi)
  !       with zeta(r) = r/(q R0)
  !---------------------------------------------------- 
  subroutine init_Bcovariants(mthis,geom)
    use globals, only : R0, B_curvature
    type(init_magnetic), intent(inout) :: mthis  
    type(geometry)     , intent(in)    :: geom
      
    integer     :: ir, itheta
    real(RKIND) :: ri, inv_Rij, zeta_ri
      
    if (B_curvature) then
      do itheta = 0,geom%Ntheta
        do ir = 0,geom%Nr
          ri                      = geom%rg(ir)
          inv_Rij                 = inv_R(ir,itheta)
          zeta_ri                 = mthis%zeta(ir)
          mthis%Br(ir,itheta)     = 0._RKIND
          mthis%Btheta(ir,itheta) = R0*inv_Rij*ri*zeta_ri 
          mthis%Bphi(ir,itheta)   = R0
        end do
      end do
    else
      do itheta = 0,geom%Ntheta
        do ir = 0,geom%Nr
          mthis%Br(ir,itheta)     = 0._RKIND
          mthis%Btheta(ir,itheta) = 0._RKIND
          mthis%Bphi(ir,itheta)   = R(ir,itheta)
        end do
      end do
    end if
  end subroutine init_Bcovariants
      
  !------------------------------------------------
  ! - initialization of the partial derivatives
  !   of the magnetic field 
  !    . dBdr(r,theta)
  !    . dBdtheta(r,theta)
  !------------------------------------------------
  subroutine init_Bderivatives(mthis,geom,init_prof)
    use globals     , only : B_curvature
    use utils_module, only : deriv1
    use init_profile_class
    type(init_magnetic), intent(inout) :: mthis  
    type(geometry)     , intent(in)    :: geom
    type(init_profile) , intent(in)    :: init_prof
    
    integer     :: ir, itheta
    real(RKIND) :: ri, hri, dhri_dr, Rij
    real(RKIND), dimension(0:geom%Nr)     :: Br_tmp
    real(RKIND), dimension(0:geom%Nr)     :: dBdr_tmp
    real(RKIND), dimension(0:geom%Ntheta) :: Btheta_tmp
    real(RKIND), dimension(0:geom%Ntheta) :: dBdtheta_tmp
      
    !*** Initialization of B ***
    if (B_curvature) then
      !-> computation of dB/dr
      do itheta = 0,geom%Ntheta
        do ir = 0,geom%Nr
          Br_tmp(ir) = mthis%B_norm(ir,itheta)
        end do
        call deriv1(Br_tmp,dBdr_tmp,geom%Nr,geom%dr,0)
        do ir = 0,geom%Nr
          mthis%dBdr(ir,itheta) = dBdr_tmp(ir)
        end do
      end do
      !-> computation of dB/dtheta
      do ir = 0,geom%Nr
        do itheta = 0,geom%Ntheta
          Btheta_tmp(itheta) = mthis%B_norm(ir,itheta)
        end do
        call deriv1(Btheta_tmp,dBdtheta_tmp,&
          geom%Ntheta,geom%dtheta,1)
        do itheta = 0,geom%Ntheta
          mthis%dBdtheta(ir,itheta) = dBdtheta_tmp(itheta)
        end do
      end do
    else
      mthis%dBdr     = 0._RKIND
      mthis%dBdtheta = 0._RKIND
    end if
  end subroutine init_Bderivatives
      
  !---------------------------------------------------- 
  ! Initialization of the magnetic configuration
  !  - initialisation of covariant components,
  !  - computation of the contravariant components,
  !  - computation of the norm,
  !  - computation of the partial derivatives and
  !
  ! Computation of jacobian in the velocity space
  !----------------------------------------------------
  subroutine init_magnetic_config(mthis,geom, &
    coord_sys,init_prof,init_curr)
    use globals, only : memory_test
    use init_profile_class
    use coord_system_class
    type(init_magnetic), intent(inout) :: mthis
    type(geometry)     , intent(in)    :: geom
    type(coord_system) , intent(in)    :: coord_sys 
    type(init_profile) , intent(in)    :: init_prof
    type(init_current) , intent(in)    :: init_curr
    
    real(RKIND) :: Ti_max, mu_max, err_threshold
      
    !*** allocation of the arrays ***
    call new_magnet_config(mthis,geom)
      
    if (.not. memory_test) then
      !*** initialisation of the magnetic field B ***
      !-> initialisation of the covariant components of B
      call init_zeta(mthis,geom,init_prof)
      call init_Bcovariants(mthis,geom)
      !-> computation of the norm of B
      call compute_norm(coord_sys,mthis%Br, &
        mthis%Btheta,mthis%Bphi,mthis%B_norm)
      !-> computation of the contravariant components of the
      !->  unit vector vec_b (with vec_b= vec_B/B)
      call compute_contravariant_vector(coord_sys, &
        mthis%Br,mthis%Btheta,mthis%Bphi, &
        mthis%B_gradr,mthis%B_gradtheta,mthis%B_gradphi)
      !-> computation of the partial derivatives of B
      call init_Bderivatives(mthis,geom,init_prof)
      !-> computation of the min and the max of B
      mthis%Bmin = minval(mthis%B_norm)
      mthis%Bmax = maxval(mthis%B_norm)
      !-> cheking of mu_max 
      if (geom%Nmu.ne.0) then
        err_threshold = 1.e-5_RKIND
        Ti_max = init_prof%Ti(0)
        mu_max = -Ti_max*mthis%Bmin*log(err_threshold)
        if (geom%mug(geom%Nmu).lt.mu_max) then
          if (pglobal_id.eq.0) then
            print*,'---> REMARK : max(mu) = ', &
              geom%mug(geom%Nmu),' < ',mu_max
          end if
        end if
      end if
    end if      
  end subroutine init_magnetic_config
      
  !*********************************************************
  ! Computation of quantities associated to Bstar
  !*********************************************************
  !---------------------------------------------------------
  ! Computation of the scalar product mu0*vec_J.vec_b where
  !  the scalar product is defined as 
  !    vec_J.vec_b = 1/B*(J_gradx1*B_gradx1 + 
  !                  J_gradx2*B_gradx2 + J_gradx3*B_gradx3) 
  !
  ! Rk: This scalar product appears in the 
  !     definition of Bstar as:
  !  Bstar = B_norm + mi*vpar*mu0*vec_J.vec_b/(e*B_norm)
  !---------------------------------------------------------
  subroutine compute_scalprod_mu0Jb(geom,mthis,init_curr, &
    ir,itheta,Sscalprod_mu0Jb)
    type(geometry)     , intent(in)  :: geom
    type(init_magnetic), intent(in)  :: mthis
    type(init_current) , intent(in)  :: init_curr
    integer            , intent(in)  :: ir
    integer            , intent(in)  :: itheta
    real(RKIND)        , intent(out) :: Sscalprod_mu0Jb
    
    real(RKIND) :: Bnorm_ij
    
    !*** compute the scalar product of mu0*vec_J ***
    !***  and vec_b                              ***
    !-> B norm
    Bnorm_ij = mthis%B_norm(ir,itheta) 
    !-> mu0*vec_J.vec_B
    Sscalprod_mu0Jb = &
      mthis%B_gradr(ir,itheta) * &
      init_curr%mu0_J_gradr(ir,itheta) + &
      mthis%B_gradtheta(ir,itheta) * &
      init_curr%mu0_J_gradtheta(ir,itheta) + &
      mthis%B_gradphi(ir,itheta) * &
      init_curr%mu0_J_gradphi(ir,itheta)
    !-> mu0*vec_J.vec_b = mu0*vec_J.vec_B/B
    Sscalprod_mu0Jb = Sscalprod_mu0Jb/Bnorm_ij
  end subroutine compute_scalprod_mu0Jb
      
  !-------------------------------------------------------
  ! Computation of Bstar 
  !   at the mesh point (ir,itheta,ivpar), as :
  !  Bstar = B_norm + mi*vpar*mu0*vec_J.vec_b/(e*B_norm)
  !-------------------------------------------------------
  subroutine compute_Bstar(geom,mthis,init_curr, &
    ir,itheta,ivpar,SBstar)
    use globals, only : Zi
    type(geometry)     , intent(in)  :: geom
    type(init_magnetic), intent(in)  :: mthis
    type(init_current) , intent(in)  :: init_curr
    integer            , intent(in)  :: ir
    integer            , intent(in)  :: itheta
    integer            , intent(in)  :: ivpar
    real(RKIND)        , intent(out) :: SBstar
      
    real(RKIND) :: vparl, Bnorm_ij
    real(RKIND) :: scalprod_mu0Jb_tmp
      
    vparl = geom%vparg(ivpar)
    !-> B norm
    Bnorm_ij = mthis%B_norm(ir,itheta) 
    !-> computation of the scalar product of vec_J and vec_b
    !->  as vec_J.vec_b = 1/B*(J_gradx1*B_gradx1 + 
    !->  J_gradx2*B_gradx2 + J_gradx3*B_gradx3) 
    call compute_scalprod_mu0Jb(geom,mthis,init_curr, &
      ir,itheta,scalprod_mu0Jb_tmp)
    !-> computation of Bstar
    SBstar = Bnorm_ij + Zi*vparl*scalprod_mu0Jb_tmp/Bnorm_ij
  end subroutine compute_Bstar
      
  !---------------------------------------------------- 
  ! Computation of the jacobian in the velocity 
  !  space, which is equal to :
  !    . 1 in the 4D case 
  !    . 2*pi*Bstar(r,theta,vpar) in the 5D case
  !----------------------------------------------------   
  subroutine compute_jacobian_velocity(geom,mthis, &
    init_curr,ir,itheta,ivpar,Sjacobian_velocity)
    use globals, only : mumin
    type(geometry)     , intent(in)  :: geom
    type(init_magnetic), intent(in)  :: mthis
    type(init_current) , intent(in)  :: init_curr
    integer            , intent(in)  :: ir
    integer            , intent(in)  :: itheta
    integer            , intent(in)  :: ivpar
    real(RKIND)        , intent(out) :: Sjacobian_velocity
      
    real(RKIND) :: Bstar_tmp
      
    if ((geom%Nmu.ne.0).or.(mumin.ne.0._RKIND)) then
#ifdef NOPRECOMPUTE 
        call compute_Bstar(geom,mthis,init_curr, &
          ir,itheta,ivpar,Bstar_tmp)
#else
        call precomputed_Bstar(ir,itheta,ivpar,Bstar_tmp)
#endif
      Sjacobian_velocity = TWOPI*Bstar_tmp
    else
      !-> case Nmu=0 and mumin=0
      Sjacobian_velocity = 1._RKIND
    end if
  end subroutine compute_jacobian_velocity
      
  !-------------------------------------------------------
  !  Computes the contravariant components of bstar
  !  (i.e bstar_gradxi) are defined as :
  !    bstar_gradxi = 1/Bstar* [B_gradxi + 
  !                   mi*vpar*mu0*J_gradxi/(e*B_norm)]
  !-------------------------------------------------------
  subroutine compute_bstar_contravariant(geom, &
    mthis,init_curr,ir,itheta,ivpar, &
    Sbstar_gradx1,Sbstar_gradx2,Sbstar_gradx3)
    use globals, only : Zi
    type(geometry)       , intent(in)  :: geom
    type(init_magnetic)  , intent(in)  :: mthis
    type(init_current)   , intent(in)  :: init_curr
    integer              , intent(in)  :: ir
    integer              , intent(in)  :: itheta
    integer              , intent(in)  :: ivpar
    real(RKIND), optional, intent(out) :: Sbstar_gradx1
    real(RKIND), optional, intent(out) :: Sbstar_gradx2
    real(RKIND), optional, intent(out) :: Sbstar_gradx3
    
    real(RKIND) :: vparl
    real(RKIND) :: Bij, Bstar_ijl
      
    vparl = geom%vparg(ivpar)
    Bij   = mthis%B_norm(ir,itheta)
    call compute_Bstar(geom,mthis,init_curr, &
      ir,itheta,ivpar,Bstar_ijl)
    !-> bstar_gradxi = 1/Bstar* [ b_gradxi + 
    !                  mi*vpar*mu0*J_gradxi/(e*B_norm) ]
    if (present(Sbstar_gradx1)) then
      Sbstar_gradx1 = ( mthis%B_gradr(ir,itheta) + &
        Zi*vparl*init_curr%mu0_J_gradr(ir,itheta)/Bij ) / &
        Bstar_ijl
    end if
    if (present(Sbstar_gradx2)) then
      Sbstar_gradx2 = ( mthis%B_gradtheta(ir,itheta) + &
        Zi*vparl*init_curr%mu0_J_gradtheta(ir,itheta) / &
        Bij ) / Bstar_ijl
    end if
    if (present(Sbstar_gradx3)) then
      Sbstar_gradx3 = ( mthis%B_gradphi(ir,itheta) + &
        Zi*vparl*init_curr%mu0_J_gradphi(ir,itheta)/Bij ) / &
        Bstar_ijl
    end if
  end subroutine compute_bstar_contravariant
      
  !-------------------------------------------------------
  !  Computes and save in 3D arrays:
  !    - Bstar and
  !    - the contravariant components of bstar
  ! (where bstar=Bstar/B)
  !-------------------------------------------------------
  subroutine init_precompute_Bstar(init_magnet,init_curr,geom, &
    coord_sys,init_prof)
    use globals, only : Bstar_3D, bstar_gradr_3D, &
      bstar_gradtheta_3D, bstar_gradphi_3D
    use init_profile_class
    type(init_magnetic)  , intent(in)    :: init_magnet
    type(init_current)   , intent(in)    :: init_curr
    type(geometry)       , intent(in)    :: geom
    type(coord_system)   , intent(in)    :: coord_sys 
    type(init_profile)   , intent(in)    :: init_prof
    
    integer     :: ir, itheta, ivpar
    real(RKIND) :: ri, bstar_ijl
    real(RKIND) :: bstar_gradr_tmp, bstar_gradth_tmp
    real(RKIND) :: bstar_gradphi_tmp
      
    !*** initialization of the covariant coordinates ***
    if (.not.memory_test) then
      do ivpar = 0, geom%Nvpar
        do itheta = 0, geom%Ntheta
          do ir = 0, geom%Nr
            call compute_Bstar(geom,init_magnet,init_curr, &
              ir,itheta,ivpar,Bstar_ijl)
            call compute_bstar_contravariant(geom, &
              init_magnet,init_curr,ir,itheta,ivpar, &
              Sbstar_gradx1=bstar_gradr_tmp, &
              Sbstar_gradx2=bstar_gradth_tmp, &
              Sbstar_gradx3=bstar_gradphi_tmp)
            Bstar_3D(ir,itheta,ivpar)           = Bstar_ijl
            bstar_gradr_3D(ir,itheta,ivpar)     = bstar_gradr_tmp
            bstar_gradtheta_3D(ir,itheta,ivpar) = bstar_gradth_tmp
            bstar_gradphi_3D(ir,itheta,ivpar)   = bstar_gradphi_tmp
          end do
        end do
      end do
    end if
  end subroutine init_precompute_bstar
end module init_magnetic_class
