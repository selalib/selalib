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
! file : Bstar_inline.f90
! date : 25/03/2010
!------------------------------------------------------
      
!-------------------------------------------------------
! Computation of Bstar 
!   at the mesh point (ir,itheta,ivpar), as :
!  Bstar = B_norm + mi*vpar*mu0*vec_J.vec_b/(e*B_norm)
!-------------------------------------------------------
subroutine precomputed_Bstar(ir,itheta,ivpar,SBstar)
  use prec_const
  use globals, only : Bstar_3D
  integer            , intent(in)  :: ir
  integer            , intent(in)  :: itheta
  integer            , intent(in)  :: ivpar
  real(RKIND)        , intent(out) :: SBstar
      
  SBstar = Bstar_3D(ir,itheta,ivpar)
end subroutine precomputed_Bstar
      
!-------------------------------------------------------
!  Computes the toroidal contravariant component of 
!  bstar defined as :
!    bstar_gradr = 1/Bstar* [B_gradr + 
!                    mi*vpar*mu0*J_gradr/(e*B_norm)]
!-------------------------------------------------------
subroutine precomputed_bstar_gradr(ir,itheta,&
  ivpar,Sbstar_gradr) 
  use prec_const
  use globals, only : bstar_gradr_3D
  integer         , intent(in)  :: ir
  integer         , intent(in)  :: itheta
  integer         , intent(in)  :: ivpar
  real(RKIND)     , intent(out) :: Sbstar_gradr
  
  Sbstar_gradr = bstar_gradr_3D(ir,itheta,ivpar)
end subroutine precomputed_bstar_gradr
      
!-------------------------------------------------------
!  Computes the toroidal contravariant component of 
!  bstar defined as :
!    bstar_gradtheta = 1/Bstar* [B_gradtheta + 
!                    mi*vpar*mu0*J_gradtheta/(e*B_norm)]
!-------------------------------------------------------
subroutine precomputed_bstar_gradtheta(ir,itheta,&
  ivpar,Sbstar_gradtheta) 
  use prec_const
  use globals, only : bstar_gradtheta_3D
  integer         , intent(in)  :: ir
  integer         , intent(in)  :: itheta
  integer         , intent(in)  :: ivpar
  real(RKIND)     , intent(out) :: Sbstar_gradtheta
  
  Sbstar_gradtheta = bstar_gradtheta_3D(ir,itheta,ivpar)
end subroutine precomputed_bstar_gradtheta
      
!-------------------------------------------------------
!  Computes the toroidal contravariant component of 
!  bstar defined as :
!    bstar_gradphi = 1/Bstar* [B_gradphi + 
!                    mi*vpar*mu0*J_gradphi/(e*B_norm)]
!-------------------------------------------------------
subroutine precomputed_bstar_gradphi(ir,itheta,&
  ivpar,Sbstar_gradphi) 
  use prec_const
  use globals, only : bstar_gradphi_3D
  integer         , intent(in)  :: ir
  integer         , intent(in)  :: itheta
  integer         , intent(in)  :: ivpar
  real(RKIND)     , intent(out) :: Sbstar_gradphi
  
  Sbstar_gradphi = bstar_gradphi_3D(ir,itheta,ivpar)
end subroutine precomputed_bstar_gradphi
