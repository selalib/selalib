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
      
!------------------------------------------------------
! file : velocities_inline.f90
! date : 15/02/2010
! Rk : all the quantities are computed using the
!  cubic spline expression
! ( Module used in the case of the Backward 
!   Semi-Lagrangian (BSL) scheme )
!------------------------------------------------------
      
!--------------------------------------------------------
! Computation of the radial contravariant component
!  of the curvature drift velocity, i.e vD.gradr as:
!    vD_gradr = K_curv*vg*[B,r]
!  with 
!    vg =  Zi(vpar^2 + mu B)/(Bstar*B)
!  and [B,r] the Poisson Bracket is defined as:
!    [B,r] = -1/(Js B)*Bphi*dB/dtheta
!  where Js is the jacobian in space
!
! Rk : K_curv=1 if B_curvature=.true. and equal 
!       to 0 otherwise. 
!--------------------------------------------------------
subroutine compute_vDgradr(geom,init_magnet, &
  init_curr,begin_dim1,end_dim1,begin_dim2,end_dim2, &
  begin_dim4,end_dim4,SvD_gradr)
  use globals           , only : mu_id, Zi, R0, K_curv
  use coord_system_class, only : jacobian_space
  use geometry_class
  use init_profile_class
  use init_magnetic_class
  type(geometry)     , intent(in)    :: geom
  type(init_magnetic), intent(in)    :: init_magnet
  type(init_current) , intent(in)    :: init_curr
  integer            , intent(in)    :: begin_dim1
  integer            , intent(in)    :: end_dim1
  integer            , intent(in)    :: begin_dim2
  integer            , intent(in)    :: end_dim2
  integer            , intent(in)    :: begin_dim4
  integer            , intent(in)    :: end_dim4
  real(RKIND), &
    dimension(0:)    , intent(inout) :: SvD_gradr
      
  integer     :: ir, itheta, ivpar
  integer     :: icount
  real(RKIND) :: vpar, vpar2, mu, vg
  real(RKIND) :: Bij, Bstar_ijl
  real(RKIND) :: dBdtheta_tmp, Bphi_tmp, J_tmp
  real(RKIND) :: PoissBrack_B_r
      
  mu     = geom%mug(mu_id)
  icount = 0
  do ivpar = begin_dim4,end_dim4
    vpar  = geom%vparg(ivpar)
    vpar2 = vpar*vpar
    do itheta = begin_dim2,end_dim2
      do ir = begin_dim1,end_dim1
        Bij = init_magnet%B_norm(ir,itheta)
#ifdef NOPRECOMPUTE 
          call compute_Bstar(geom,init_magnet,init_curr, &
            ir,itheta,ivpar,Bstar_ijl)
#else
          call precomputed_Bstar(ir,itheta,ivpar,Bstar_ijl)
#endif
        !-> vg = Zi (vpar^2 + mu B)/(B*Bstar)
        vg  = real(Zi)*(vpar2 + mu*Bij)/(Bij*Bstar_ijl)
        !-> dB/dtheta and Bphi
        dBdtheta_tmp = init_magnet%dBdtheta(ir,itheta)
        Bphi_tmp     = init_magnet%Bphi(ir,itheta)
        !-> [B,r] = -1/(Js B)*Bphi*dB/dtheta
        J_tmp          = jacobian_space(ir,itheta)
        PoissBrack_B_r = -1._RKIND/(J_tmp*Bij) * &
          Bphi_tmp*dBdtheta_tmp
        !-> vD.gradr =  vg*[B,r]
        SvD_gradr(icount) = K_curv*vg*PoissBrack_B_r
        icount = icount + 1
      end do
    end do
  end do
end subroutine compute_vDgradr
      
!--------------------------------------------------------
! Computation of the poloidal contravariant component
!  of the curvature drift velocity, i.e vD.gradtheta as:
!    vD_gradtheta = K_curv*vg*[B,theta]
!  with 
!    vg =  Zi(vpar^2 + mu B)/(Bstar*B)
!  and [B,theta] the Poisson Bracket is defined as:
!    [B,theta] = 1/(Js B)*Bphi*dB/dr
!  where Js is the jacobian in space
!
! Rk : K_curv=1 if B_curvature=.true. and equal 
!       to 0 otherwise. 
!--------------------------------------------------------
subroutine compute_vDgradtheta(geom,init_magnet,init_curr, &
  begin_dim1,end_dim1,begin_dim2,end_dim2, &
  begin_dim4,end_dim4,SvD_gradtheta)
  use globals, only : mu_id, Zi, R0, K_curv
  use geometry_class
  use init_profile_class
  use init_magnetic_class
  type(geometry)     , intent(in)    :: geom
  type(init_magnetic), intent(in)    :: init_magnet
  type(init_current) , intent(in)    :: init_curr
  integer            , intent(in)    :: begin_dim1
  integer            , intent(in)    :: end_dim1
  integer            , intent(in)    :: begin_dim2
  integer            , intent(in)    :: end_dim2
  integer            , intent(in)    :: begin_dim4
  integer            , intent(in)    :: end_dim4
  real(RKIND), &
    dimension(0:)    , intent(inout) :: SvD_gradtheta
      
  integer     :: ir, itheta, ivpar
  integer     :: icount
  real(RKIND) :: vpar, vpar2, mu, vg
  real(RKIND) :: Bij, Bstar_ijl
  real(RKIND) :: dBdr_tmp, Bphi_tmp, J_tmp
  real(RKIND) :: PoissBrack_B_theta
      
  mu     = geom%mug(mu_id)
  icount = 0
  do ivpar = begin_dim4,end_dim4
    vpar  = geom%vparg(ivpar)
    vpar2 = vpar*vpar
    do itheta = begin_dim2,end_dim2
      do ir = begin_dim1,end_dim1
        Bij = init_magnet%B_norm(ir,itheta)
#ifdef NOPRECOMPUTE 
          call compute_Bstar(geom,init_magnet,init_curr, &
            ir,itheta,ivpar,Bstar_ijl)
#else
          call precomputed_Bstar(ir,itheta,ivpar,Bstar_ijl)
#endif
        !-> vg = Zi (vpar^2 + mu B)/(B*Bstar)
        vg  = real(Zi)*(vpar2 + mu*Bij)/(Bij*Bstar_ijl)
        !-> dB/dr and Bphi
        dBdr_tmp = init_magnet%dBdr(ir,itheta)
        Bphi_tmp = init_magnet%Bphi(ir,itheta)
        !-> [B,theta] = 1/(Js B)*Bphi*dB/dr
        J_tmp              = jacobian_space(ir,itheta)
        PoissBrack_B_theta = 1._RKIND/(J_tmp*Bij) * &
          Bphi_tmp*dBdr_tmp
        !-> vD.gradtheta =  vg*[B,theta]
        SvD_gradtheta(icount) = K_curv*vg*PoissBrack_B_theta
        icount = icount + 1
      end do
    end do
  end do
end subroutine compute_vDgradtheta
      
!--------------------------------------------------------
! Computation of the toroidal contravariant component
!  of the curvature drift velocity, i.e vD.gradphi as:
!    vD_gradphi = K_curv*vg*[B,phi]
!  with 
!    vg =  Zi(vpar^2 + mu B)/(Bstar*B)
!  and [B,phi] the Poisson Bracket is defined as:
!    [B,phi] = 1/(Js B)*(Br*dB/dtheta-Btheta*dB/dr)
!  where Js is the jacobian in space
!
! Rk : K_curv=1 if B_curvature=.true. and equal 
!       to 0 otherwise. 
!--------------------------------------------------------
subroutine compute_vDgradphi(geom,init_prof,init_magnet, &
  init_curr,begin_dim1,end_dim1,begin_dim2,end_dim2, &
  begin_dim4,end_dim4,SvD_gradphi)
  use globals, only : mu_id, Zi, R0, K_curv
  use geometry_class
  use init_profile_class
  use init_magnetic_class
  type(geometry)     , intent(in)    :: geom
  type(init_profile) , intent(in)    :: init_prof
  type(init_magnetic), intent(in)    :: init_magnet
  type(init_current) , intent(in)    :: init_curr
  integer            , intent(in)    :: begin_dim1
  integer            , intent(in)    :: end_dim1
  integer            , intent(in)    :: begin_dim2
  integer            , intent(in)    :: end_dim2
  integer            , intent(in)    :: begin_dim4
  integer            , intent(in)    :: end_dim4
  real(RKIND), &
    dimension(0:)    , intent(inout) :: SvD_gradphi
      
  integer     :: ir, itheta, ivpar
  integer     :: icount
  real(RKIND) :: vpar, vpar2, mu, vg
  real(RKIND) :: Bij, Bstar_ijl
  real(RKIND) :: dBdr_tmp, dBdtheta_tmp
  real(RKIND) :: Br_tmp, Btheta_tmp, J_tmp
  real(RKIND) :: PoissBrack_B_phi
      
  mu     = geom%mug(mu_id)
  icount = 0
  do ivpar = begin_dim4,end_dim4
    vpar  = geom%vparg(ivpar)
    vpar2 = vpar*vpar
    do itheta = begin_dim2,end_dim2
      do ir = begin_dim1,end_dim1
        Bij = init_magnet%B_norm(ir,itheta)
#ifdef NOPRECOMPUTE 
          call compute_Bstar(geom,init_magnet,init_curr, &
            ir,itheta,ivpar,Bstar_ijl)
#else
          call precomputed_Bstar(ir,itheta,ivpar,Bstar_ijl)
#endif
        !-> vg = Zi (vpar^2 + mu B)/(B*Bstar)
        vg  = real(Zi)*(vpar2 + mu*Bij)/(Bij*Bstar_ijl)
        !-> dB/dr and dB/dtheta
        dBdr_tmp     = init_magnet%dBdr(ir,itheta)
        dBdtheta_tmp = init_magnet%dBdtheta(ir,itheta)
        !-> Br and Btheta
        Br_tmp     = init_magnet%Br(ir,itheta)
        Btheta_tmp = init_magnet%Btheta(ir,itheta)
        !-> [B,phi] = 1/(Js B)*(Br*dB/dtheta-Btheta*dB/dr)
        J_tmp            = jacobian_space(ir,itheta)
        PoissBrack_B_phi = 1._RKIND/(J_tmp*Bij) * &
          (Br_tmp*dBdtheta_tmp-Btheta_tmp*dBdr_tmp)
        !-> vD.gradphi =  vg*[B,phi]
        SvD_gradphi(icount) = K_curv*vg*PoissBrack_B_phi
        icount = icount + 1
      end do
    end do
  end do
end subroutine compute_vDgradphi
      
!-------------------------------------------------------------
! Computation of the parallel gradient of Phi
!  (used for dvpar/dt advection)
!     gradpar_J0Phi = bstar.grad(J0.Phi)
!                   = bstar.gradx1*d(J0.Phi)/dx1 + 
!                     bstar.gradx2*d(J0.Phi)/dx2 + 
!                     bstar.gradx3*d(J0.Phi)/dx3 
!   where the contravariant components of bstar
!    (i.e bstar_gradxi) are defined as :
!    bstar_gradxi = 1/Bstar* [B_gradxi + 
!                   mi*vpar*mu0*J_gradxi/(e*B_norm)]
!-------------------------------------------------------------
subroutine compute_gradpar_J0Phi(geom, &
  init_magnet,init_curr,begin_dim1,end_dim1, &
  begin_dim2,end_dim2,begin_dim3,end_dim3, &
  begin_dim4,end_dim4,Sgradpar_J0Phi)
  use globals, only : &
    dJ0Phidr, dJ0Phidtheta, dJ0Phidphi
  use geometry_class
  use init_magnetic_class
  use init_current_class
  type(geometry)     , intent(in)    :: geom
  type(init_magnetic), intent(in)    :: init_magnet
  type(init_current) , intent(in)    :: init_curr
  integer            , intent(in)    :: begin_dim1
  integer            , intent(in)    :: end_dim1
  integer            , intent(in)    :: begin_dim2
  integer            , intent(in)    :: end_dim2
  integer            , intent(in)    :: begin_dim3
  integer            , intent(in)    :: end_dim3
  integer            , intent(in)    :: begin_dim4
  integer            , intent(in)    :: end_dim4
  real(RKIND), &
    dimension(0:)    , intent(inout) :: Sgradpar_J0Phi
      
  integer     :: ir, itheta, iphi, ivpar
  integer     :: icount
  real(RKIND) :: bstar_gradr_tmp, dJ0Phidr_tmp
  real(RKIND) :: bstar_gradth_tmp, dJ0Phidtheta_tmp
  real(RKIND) :: bstar_gradphi_tmp, dJ0Phidphi_tmp
      
  icount = 0
  do ivpar = begin_dim4,end_dim4
    do iphi = begin_dim3,end_dim3
      do itheta = begin_dim2,end_dim2
        do ir = begin_dim1,end_dim1
          !-> bstar_gradxi = 1/Bstar* [ B_gradxi + 
          !                  mi*vpar*mu0*J_gradxi/(e*B_norm) ]
#ifdef NOPRECOMPUTE 
            call compute_bstar_contravariant(geom, &
              init_magnet,init_curr,ir,itheta,ivpar, &
              Sbstar_gradx1=bstar_gradr_tmp, &
              Sbstar_gradx2=bstar_gradth_tmp, &
              Sbstar_gradx3=bstar_gradphi_tmp)
#else
            call precomputed_bstar_gradr(ir,itheta,ivpar, &
              bstar_gradr_tmp)
            call precomputed_bstar_gradtheta(ir,itheta,ivpar, &
              bstar_gradth_tmp)
            call precomputed_bstar_gradphi(ir,itheta,ivpar, &
              bstar_gradphi_tmp)
#endif
          !-> d(J0.Phi)/dr, d(J0.Phi)/dtheta, d(J0.Phi)/dphi
          dJ0Phidr_tmp            = dJ0Phidr(ir,itheta,iphi)
          dJ0Phidtheta_tmp        = dJ0Phidtheta(ir,itheta,iphi)
          dJ0Phidphi_tmp          = dJ0Phidphi(ir,itheta,iphi)
          Sgradpar_J0Phi(icount)  = &
            dJ0Phidr_tmp*bstar_gradr_tmp + &
            dJ0Phidtheta_tmp*bstar_gradth_tmp + &
            dJ0Phidphi_tmp*bstar_gradphi_tmp
          icount = icount+1
        end do
      enddo
    end do
  end do
end subroutine compute_gradpar_J0Phi
      
!--------------------------------------------------------------
!  Computation of vExB_gradr by using the expression:
!     vExB_gradr = 1/Bstar * [Phi,r]
!  where the Poisson bracket [Phi,r] corresponds to
!    [Phi,r] = 1/(Js B)*( Btheta*dPhi/dphi
!                 - Bphi*dPhi/dtheta )
!  where Js is the jacobian in space
! 
!  Rk : This subroutine can be called for computing the
!   ExB drift of the guiding centers or of the particles :
!    -> J0.Phi (for guiding center) in the advections
!    -> Phi (for particles) in the computation of the fluxes
!--------------------------------------------------------------
subroutine compute_vExBgradr(geom,init_magnet,init_curr, &
  dPhidtheta_3D,dPhidphi_3D, &
  begin_dim1,end_dim1,begin_dim2,end_dim2, &
  begin_dim3,end_dim3,begin_dim4,end_dim4,SvExB_gradr)
  use geometry_class
  use init_magnetic_class
  use init_current_class
  type(geometry)     , intent(in)    :: geom
  type(init_magnetic), intent(in)    :: init_magnet
  type(init_current) , intent(in)    :: init_curr
  real(RKIND), &
     dimension(:,:,:), pointer       :: dPhidtheta_3D
  real(RKIND), &
     dimension(:,:,:), pointer       :: dPhidphi_3D
  integer            , intent(in)    :: begin_dim1
  integer            , intent(in)    :: end_dim1
  integer            , intent(in)    :: begin_dim2
  integer            , intent(in)    :: end_dim2
  integer            , intent(in)    :: begin_dim3
  integer            , intent(in)    :: end_dim3
  integer            , intent(in)    :: begin_dim4
  integer            , intent(in)    :: end_dim4
  real(RKIND), &
    dimension(0:)    , intent(inout) :: SvExB_gradr
      
  integer     :: ir, itheta, iphi, ivpar
  integer     :: icount
  real(RKIND) :: Bij, Bstar_ijl
  real(RKIND) :: J_tmp
  real(RKIND) :: dPhidtheta_tmp, dPhidphi_tmp
  real(RKIND) :: Btheta_tmp, Bphi_tmp
  real(RKIND) :: PoissBrack_Phi_r
      
  !*** computation of vExB_gradr as:                ***
  !***   vExB_gradr = 1/Bstar * [Phi,r]             ***
  !*** with [Phi,r] = 1/(Js B)*( Btheta*dPhi/dphi   ***
  !***                - Bphi*dPhi/dtheta )          ***
  icount = 0
  do ivpar = begin_dim4,end_dim4
    do iphi = begin_dim3,end_dim3
      do itheta = begin_dim2,end_dim2
        do ir = begin_dim1,end_dim1
          Bij = init_magnet%B_norm(ir,itheta)
#ifdef NOPRECOMPUTE 
            call compute_Bstar(geom,init_magnet,init_curr, &
              ir,itheta,ivpar,Bstar_ijl)
#else
            call precomputed_Bstar(ir,itheta,ivpar,Bstar_ijl)
#endif
          !-> derivatives of Phi
          dPhidtheta_tmp = dPhidtheta_3D(ir,itheta,iphi)
          dPhidphi_tmp   = dPhidphi_3D(ir,itheta,iphi)
          !-> Btheta and Bphi
          Btheta_tmp = init_magnet%Btheta(ir,itheta)
          Bphi_tmp   = init_magnet%Bphi(ir,itheta)
          !-> [Phi,r] = 1/(Js B) * 
          !     ( Btheta*dPhi/dphi- Bphi*dPhi/dtheta ) 
          J_tmp            = jacobian_space(ir,itheta)
          PoissBrack_Phi_r = 1._RKIND/(J_tmp*Bij) * &
            (Btheta_tmp*dPhidphi_tmp-Bphi_tmp*dPhidtheta_tmp)
          !-> vEXB.gradr = 1/Bstar * [Phi,r]
          SvExB_gradr(icount) = PoissBrack_Phi_r/Bstar_ijl
          icount              = icount+1
        end do
      enddo
    end do
  end do
end subroutine compute_vExBgradr
      
!-----------------------------------------------------------
!  Computation of vExB_gradtheta
!     vExB_gradtheta = 1/Bstar * [Phi,theta]
!  where the Poisson bracket [Phi,theta] corresponds to
!    [Phi,theta] = 1/(Js B)*( -Br*dPhi/dphi
!                  + Bphi*dPhi/dr )
!  where Js is the jacobian in space
!
!  Rk : This subroutine can be called for computing the
!   ExB drift of the guiding centers or of the particles :
!    -> J0.Phi (for guiding center) in the advections
!    -> Phi (for particles) in the computation of the fluxes
!-----------------------------------------------------------
subroutine compute_vExBgradtheta(geom, &
  init_magnet,init_curr,dPhidr_3D,dPhidphi_3D, &
  begin_dim1,end_dim1,begin_dim2,end_dim2, &
  begin_dim3,end_dim3,begin_dim4,end_dim4,SvExB_gradtheta)
  use geometry_class
  use init_magnetic_class
  use init_current_class
  type(geometry)     , intent(in)    :: geom
  type(init_magnetic), intent(in)    :: init_magnet
  type(init_current) , intent(in)    :: init_curr
  real(RKIND), &
     dimension(:,:,:), pointer       :: dPhidr_3D
  real(RKIND), &
     dimension(:,:,:), pointer       :: dPhidphi_3D
  integer            , intent(in)    :: begin_dim1
  integer            , intent(in)    :: end_dim1
  integer            , intent(in)    :: begin_dim2
  integer            , intent(in)    :: end_dim2
  integer            , intent(in)    :: begin_dim3
  integer            , intent(in)    :: end_dim3
  integer            , intent(in)    :: begin_dim4
  integer            , intent(in)    :: end_dim4
  real(RKIND), &
    dimension(0:)    , intent(inout) :: SvExB_gradtheta
      
  integer     :: ir, itheta, iphi, ivpar
  integer     :: icount
  real(RKIND) :: Bij, Bstar_ijl
  real(RKIND) :: J_tmp
  real(RKIND) :: dPhidr_tmp, dPhidphi_tmp
  real(RKIND) :: Br_tmp, Bphi_tmp
  real(RKIND) :: PoissBrack_Phi_theta
      
  !*** computation of vExB_gradtheta as:         ***
  !***   vExB_gradtheta = 1/Bstar * [Phi,theta]  ***
  !*** with                                      ***    
  !***   [Phi,r] = 1/(Js B)*( -Br*dPhi/dphi      ***
  !***             + Bphi*dPhi/dr )              ***
  icount = 0
  do ivpar = begin_dim4,end_dim4
    do iphi = begin_dim3,end_dim3
      do itheta = begin_dim2,end_dim2
        do ir = begin_dim1,end_dim1
          Bij = init_magnet%B_norm(ir,itheta)
#ifdef NOPRECOMPUTE 
            call compute_Bstar(geom,init_magnet,init_curr, &
              ir,itheta,ivpar,Bstar_ijl)
#else
            call precomputed_Bstar(ir,itheta,ivpar,Bstar_ijl)
#endif
          !-> derivatives of Phi
          dPhidr_tmp   = dPhidr_3D(ir,itheta,iphi)
          dPhidphi_tmp = dPhidphi_3D(ir,itheta,iphi)
          !-> Br and Bphi
          Br_tmp   = init_magnet%Br(ir,itheta)
          Bphi_tmp = init_magnet%Bphi(ir,itheta)
          !-> [Phi,theta] = 1/(Js B) *
          !      ( -Br*dPhi/dphi+ Bphi*dPhi/dr )
          J_tmp                = jacobian_space(ir,itheta)
          PoissBrack_Phi_theta = 1._RKIND/(J_tmp*Bij) * &
            (-Br_tmp*dPhidphi_tmp + Bphi_tmp*dPhidr_tmp)
          !-> vEXB.gradtheta = 1/Bstar * [Phi,theta]
          SvExB_gradtheta(icount) = PoissBrack_Phi_theta / &
            Bstar_ijl
          icount                  = icount+1
        end do
      enddo
    end do
  end do
end subroutine compute_vExBgradtheta
      
!--------------------------------------------------------
!  Computation of vExB_gradphi
!     vExB_gradphi = 1/Bstar * [Phi,phi]
!  where the Poisson bracket [Phi,phi] corresponds to
!    [Phi,phi] = 1/(Js B)*( Br*dPhi/dtheta
!                  - Btheta*dPhi/dr )
!  where Js is the jacobian in space
!
!  Rk : This subroutine can be called for computing the
!   ExB drift of the guiding centers or of the particles :
!    -> J0.Phi (for guiding center) in the advections
!    -> Phi (for particles) in the computation of the fluxes
!--------------------------------------------------------
subroutine compute_vExBgradphi(geom, &
  init_magnet,init_curr,dPhidr_3D,dPhidtheta_3D, &
  begin_dim1,end_dim1,begin_dim2,end_dim2, &
  begin_dim3,end_dim3,begin_dim4,end_dim4,SvExB_gradphi)
  use geometry_class
  use init_magnetic_class
  use init_current_class
      
  type(geometry)     , intent(in)    :: geom
  type(init_magnetic), intent(in)    :: init_magnet
  type(init_current) , intent(in)    :: init_curr
  real(RKIND), &
     dimension(:,:,:), pointer       :: dPhidr_3D
  real(RKIND), &
     dimension(:,:,:), pointer       :: dPhidtheta_3D
  integer            , intent(in)    :: begin_dim1
  integer            , intent(in)    :: end_dim1
  integer            , intent(in)    :: begin_dim2
  integer            , intent(in)    :: end_dim2
  integer            , intent(in)    :: begin_dim3
  integer            , intent(in)    :: end_dim3
  integer            , intent(in)    :: begin_dim4
  integer            , intent(in)    :: end_dim4
  real(RKIND), &
    dimension(0:)    , intent(inout) :: SvExB_gradphi
      
  integer     :: ir, itheta, iphi, ivpar
  integer     :: icount
  real(RKIND) :: Bij, Bstar_ijl
  real(RKIND) :: J_tmp
  real(RKIND) :: dPhidr_tmp, dPhidtheta_tmp
  real(RKIND) :: Br_tmp, Btheta_tmp
  real(RKIND) :: PoissBrack_Phi_phi
      
  !*** computation of vExB_gradphi as:              ***
  !***   vExB_gradphi = 1/Bstar * [Phi,phi]         ***
  !*** with [Phi,phi] = 1/(Js B)*( Br*dPhi/dtheta   ***
  !***                   - Btheta*dPhi/dr )         ***
  icount = 0
  do ivpar = begin_dim4,end_dim4
    do iphi = begin_dim3,end_dim3
      do itheta = begin_dim2,end_dim2
        do ir = begin_dim1,end_dim1
          Bij = init_magnet%B_norm(ir,itheta)
#ifdef NOPRECOMPUTE 
            call compute_Bstar(geom,init_magnet,init_curr, &
              ir,itheta,ivpar,Bstar_ijl)
#else
            call precomputed_Bstar(ir,itheta,ivpar,Bstar_ijl)
#endif
          !-> derivatives of Phi
          dPhidr_tmp     = dPhidr_3D(ir,itheta,iphi)
          dPhidtheta_tmp = dPhidtheta_3D(ir,itheta,iphi)
          !-> Br and Btheta
          Br_tmp     = init_magnet%Br(ir,itheta)
          Btheta_tmp = init_magnet%Btheta(ir,itheta)
          !-> [Phi,phi] = 1/(Js B)*( Br*dPhi/dtheta
          !                - Btheta*dPhi/dr 
          J_tmp                = jacobian_space(ir,itheta)
          PoissBrack_Phi_phi = 1._RKIND/(J_tmp*Bij) * &
            (Br_tmp*dPhidtheta_tmp - Btheta_tmp*dPhidr_tmp)
          !-> vEXB.gradphi = 1/Bstar * [Phi,phi]
          SvExB_gradphi(icount) = PoissBrack_Phi_phi/Bstar_ijl
          icount                = icount+1
        end do
      enddo
    end do
  end do
end subroutine compute_vExBgradphi
      
!-----------------------------------------------------
! Computation of the parallel projection of B
!    vec_bstar.gradB = bstar_gradx1*partial_x1 B + 
!                      bstar_gradx2*partial_x2 B 
!  where vec_bstar = vec_Bstar/Bstar
!  where the contravariant components of bstar
!  (i.e bstar_gradxi) are defined as :
!    bstar_gradxi = 1/Bstar* [B_gradxi + 
!                   mi*vpar*mu0*J_gradxi/(e*B_norm)]
!-----------------------------------------------------
subroutine compute_gradpar_B(geom,mthis,init_curr, &
  begin_dim1,end_dim1,begin_dim2,end_dim2, &
  begin_dim4,end_dim4,Sgradpar_B)
  use globals, only : Zi
  use geometry_class
  use init_magnetic_class
  use init_current_class
  type(geometry)     , intent(in)    :: geom
  type(init_magnetic), intent(in)    :: mthis
  type(init_current) , intent(in)    :: init_curr
  integer            , intent(in)    :: begin_dim1
  integer            , intent(in)    :: end_dim1
  integer            , intent(in)    :: begin_dim2
  integer            , intent(in)    :: end_dim2
  integer            , intent(in)    :: begin_dim4
  integer            , intent(in)    :: end_dim4
  real(RKIND), &
    dimension(0:)    , intent(inout) :: Sgradpar_B
  
  integer     :: ir, itheta, ivpar  
  integer     :: icount
  real(RKIND) :: bstar_gradr_tmp, bstar_gradtheta_tmp
  real(RKIND) :: dBdr_tmp, dBdtheta_tmp
      
  icount = 0
  do ivpar = begin_dim4,end_dim4
    do itheta = begin_dim2,end_dim2
      do ir = begin_dim1,end_dim1
        !-> bstar_gradxi = 1/Bstar* [ b_gradxi + 
        !               mi*vpar*mu0*J_gradxi/(e*B_norm) ]
#ifdef NOPRECOMPUTE 
          call compute_bstar_contravariant(geom, &
            mthis,init_curr,ir,itheta,ivpar, &
            Sbstar_gradx1=bstar_gradr_tmp, &
            Sbstar_gradx2=bstar_gradtheta_tmp)
#else
          call precomputed_bstar_gradr(ir,itheta,ivpar, &
            bstar_gradr_tmp)
          call precomputed_bstar_gradtheta(ir,itheta,ivpar, &
            bstar_gradtheta_tmp)
#endif
        !-> dB/dr and dB/dtheta
        dBdr_tmp         = mthis%dBdr(ir,itheta)
        dBdtheta_tmp     = mthis%dBdtheta(ir,itheta)
        !->  vec_bstar.gradB = bstar_gradr*dBdr +
        !                      bstar_gradtheta*dBdtheta  
        Sgradpar_B(icount) = bstar_gradr_tmp*dBdr_tmp + &
          bstar_gradtheta_tmp*dBdtheta_tmp
        icount = icount + 1
      end do
    end do
  end do
end subroutine compute_gradpar_B
