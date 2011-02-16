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
      
!---------------------------------------------------------
! file : advec1D_BSL.f90
! date : 31/10/2000
!   1D advection in phi and vparallel directions used in
!   a Backward Semi-Lagrangian method
!---------------------------------------------------------
module advec1D_BSL_module
  use globals, only : Nr, Ntheta, Nphi, Nvpar, Nbthread
  use mem_alloc_module
  use integration_module
  use clock_module
  use geometry_class
  use init_profile_class
  use init_magnetic_class
  use fdistribu5d_class
  use interpolation_module
  
  implicit none
      
  include "Bstar_inline.h"
  include "velocities_inline.h"
  !******************************
  contains
  !******************************
#include "Bstar_inline.f90"
#include "velocities_inline.f90"
      
  !------------------------------------------------------------
  ! splitting in phi direction 
  !  - Solving of : df/dt+vpar/R*df/dphi=0
  !     . dphi = vpar bstar_gradphi + vExB_gradphi + vD_gradphi
  !     . f**(r,theta,f,vpar,mu) = f*(r,theta,phi-dphi,vpar,mu)
  !  with bstar_gradphi = 1/Bstar* [B_gradphi + 
  !                   mi*vpar*mu0*J_gradphi/(e*B_norm)]
  !------------------------------------------------------------
  subroutine advec1D_BSL_phi(geom,dt,init_prof, &
    init_magnet,init_curr,f,maxalp_phi)
    use globals          , only : istart, iend, jstart, jend, &
      dJ0Phidr, dJ0Phidtheta
#ifdef _OPENMP
    use OMPutils_module  , only : omp_get_thread_num
#endif
    use OMPutils_module  , only : Romp_scoefphi, Romp1_0Nphi, &
      Romp2_0Nphi, Romp3_0Nphi, Romp1_0Nvpar
!R3 #include "r3_info.h" !R3
    type(geometry)     , intent(in)    :: geom
    real(RKIND)        , intent(in)    :: dt
    type(init_profile) , intent(in)    :: init_prof
    type(init_magnetic), intent(in)    :: init_magnet
    type(init_current) , intent(in)    :: init_curr
    type(fdistribu5d)  , intent(inout) :: f
    real(RKIND)        , intent(out)   :: maxalp_phi
      
    integer     :: ir, itheta, iphi, ivpar
    real(RKIND) :: phimin, phimax, phistar, disp_phi
    real(RKIND) :: vparl
    real(RKIND) :: bstar_gradphi_tmp
    real(RKIND) :: vExB_gradphi_tmp, vD_gradphi_tmp
    real(RKIND) :: finterpol, nbLphi
    integer     :: tid
    
    real(RKIND), dimension(:), pointer :: rhsphi
    real(RKIND), dimension(:), pointer :: phig
    real(RKIND), dimension(:), pointer :: scoefphi
      
    real(RKIND), dimension(:), pointer :: vExB_gradphi_1D
    real(RKIND), dimension(:), pointer :: vD_gradphi_1D
      
    real(RKIND), dimension(1:Nbthread) :: Romp_maxalp_phi
      
!R3 call r3_info_begin (r3_info_index_0, 'CPU_advec1d_sl_phi') !R3
    call clck_time(bclock_advec1Dphi)
      
#ifdef _OPENMP
!$OMP PARALLEL private(ir,itheta,iphi,ivpar, &
!$OMP   phimin,phimax,phistar,disp_phi,vparl, &
!$OMP   bstar_gradphi_tmp,vExB_gradphi_tmp,vD_gradphi_tmp, &
!$OMP   vExB_gradphi_1D, vD_gradphi_1D, nbLphi, &
!$OMP   finterpol,rhsphi,scoefphi,phig,tid) default(shared)
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
      
    Romp_maxalp_phi(tid) = 0._RKIND
    scoefphi             => Romp_scoefphi(tid)%val
    phig                 => Romp1_0Nphi(tid)%val
    rhsphi               => Romp2_0Nphi(tid)%val
    vExB_gradphi_1D      => Romp3_0Nphi(tid)%val
    vD_gradphi_1D        => Romp1_0Nvpar(tid)%val
      
    !-> initialisation of phig
    do iphi = 0,Nphi
       phig(iphi) = geom%phig(iphi)
    end do
      
    phimin = geom%phig(0)
    phimax = geom%phig(Nphi)
      
!$OMP DO SCHEDULE(static)
    do itheta = jstart,jend 
      do ir = f%istart_modif,f%iend_modif
        !-> computation of vD_gradphi
        call compute_vDgradphi(geom,init_prof,init_magnet, &
          init_curr,ir,ir,itheta,itheta,0,f%n4,vD_gradphi_1D)
        do ivpar = 0,f%n4
          vparl          = geom%vparg(ivpar)
          vD_gradphi_tmp = vD_gradphi_1D(ivpar) 
          do iphi = 0,f%n3-1
            rhsphi(iphi) = f%values(ir,itheta,iphi,ivpar)
          end do
          rhsphi(Nphi) = rhsphi(0)
          !*** cubic spline computation ***
          call period_omp_spline_coef(f%pspline1d_phi, &
            rhsphi,scoefphi)
          !*** splitting in phi direction ***
          !-> computation of bstar.gradphi
#ifdef NOPRECOMPUTE 
            call compute_bstar_contravariant(geom, &
              init_magnet,init_curr,ir,itheta,ivpar, &
              Sbstar_gradx3=bstar_gradphi_tmp)
#else
            call precomputed_bstar_gradphi(ir,itheta,ivpar, &
              bstar_gradphi_tmp)
#endif
          !-> computation of vExB_gradphi
          call compute_vExBgradphi(geom,init_magnet,init_curr, &
            dJ0Phidr,dJ0Phidtheta,ir,ir,itheta,itheta, &
            0,f%n3,ivpar,ivpar,vExB_gradphi_1D)
          do iphi = 0,f%n3-1
            vExB_gradphi_tmp = vExB_gradphi_1D(iphi) 
            disp_phi         = dt * (vparl*bstar_gradphi_tmp + &
              vExB_gradphi_tmp + vD_gradphi_tmp)
            Romp_maxalp_phi(tid) = max(Romp_maxalp_phi(tid), &
              abs(disp_phi)/geom%dphi)
            phistar    = geom%phig(iphi) - disp_phi	  
            nbLphi     = 0._RKIND
            if (phistar.ge.phimax) then
              nbLphi  = abs(floor((phistar-phimin)/geom%Lphi))
              phistar = phistar - nbLphi*geom%Lphi  
            elseif (phistar.lt.phimin) then
              nbLphi  = -abs(floor((phistar-phimin)/geom%Lphi))
              phistar = phistar - nbLphi*geom%Lphi  
            endif
            phistar = MIN(MAX(phistar,phimin),phimax)
            !*** f interpolation ***
            call interpol1d_phi(geom,scoefphi, &
              phig,phistar,finterpol)
            f%values(ir,itheta,iphi,ivpar) = finterpol
          end do
        end do
      end do
    end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
      
    maxalp_phi = maxval(Romp_maxalp_phi(1:Nbthread))
      
    call clck_time(eclock_advec1Dphi)
    call clck_diff(bclock_advec1Dphi,eclock_advec1Dphi, &
      global_time_advec1Dphi)
!R3 call r3_info_end (r3_info_index_0) !R3
  end subroutine advec1D_BSL_phi
      
  !------------------------------------------------------------
  ! splitting in v parallel direction
  !  - Solving of : df/dt+(dvpar/dt)*df/dvpar=0 as:
  !     . f*(r,theta,f,vpar,mu) = fn(r,theta,phi,vpar-dpar,mu)
  !    with
  !     . dpar = dt*[ -Zi*gradpar(J0.Phi) - mu*gradpar B + 
  !              K_curv*vpar/B*(vExB.grad_B) ]
  !    where
  !        -> gradpar_J0Phi = vec_bstar.grad_Phi
  !        -> gradpar_B     = vec_bstar.grad_B
  !        -> vExB.grad_B   = vExB_gradr*dB/dr + 
  !                         vExB_gradtheta*dB/dtheta
  ! Rk : K_curv=1 if B_curvature=.true. and equal 
  !       to 0 otherwise. 
  !------------------------------------------------------------
  subroutine advec1D_BSL_vpar(geom,dt,init_prof, &
    init_magnet,init_curr,f,maxalp_vpar)
    use globals        , only : Zi, K_curv, &
      istart, iend, jstart, jend, &
      dJ0Phidr, dJ0Phidtheta, dJ0Phidphi
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp_scoefvpar, &
      Romp1_0Nvpar, Romp2_0Nvpar, Romp3_0Nvpar, &
      Romp4_0Nvpar, Romp5_0Nvpar, Romp6_0Nvpar
!R3 #include "r3_info.h" !R3
    type(geometry)       , intent(in)    :: geom
    real(RKIND)          , intent(in)    :: dt
    type(init_profile)   , intent(in)    :: init_prof
    type(init_magnetic)  , intent(in)    :: init_magnet
    type(init_current)   , intent(in)    :: init_curr
    type(fdistribu5d)    , intent(inout) :: f
    real(RKIND)          , intent(out)   :: maxalp_vpar
      
    integer     :: ir, itheta, iphi, ivpar
    real(RKIND) :: mu_value, vstar, disp_vpar, vparl
    real(RKIND) :: inv_Bij
    real(RKIND) :: finterpol
    real(RKIND) :: dBdr_tmp, dBdtheta_tmp
    real(RKIND) :: gradpar_J0Phi_tmp, gradpar_B_tmp
    real(RKIND) :: vExB_gradB_tmp
    integer     :: tid
      
    real(RKIND), dimension(:), pointer :: rhsvpar
    real(RKIND), dimension(:), pointer :: scoefvpar
    real(RKIND), dimension(:), pointer :: vparg
      
    !-> components of ExB drift velocity
    real(RKIND), dimension(:), pointer :: vExB_gradr_1D
    real(RKIND), dimension(:), pointer :: vExB_gradtheta_1D
    !-> gradpar.B
    real(RKIND), dimension(:), pointer :: gradpar_B_1D
    !-> gradpar.Phi
    real(RKIND), dimension(:), pointer :: gradpar_J0Phi_1D
      
    real(RKIND), dimension(1:Nbthread) :: Romp_maxalp_vpar
      
!R3 call r3_info_begin (r3_info_index_0, 'CPU_advec1d_sl_vpar') !R3
    call clck_time(bclock_advec1Dvpar)
      
    !*** initialisation of mu value ***
    mu_value = f%mu_value
      
#ifdef _OPENMP
!$OMP PARALLEL private(ir,itheta,iphi,ivpar,vstar,disp_vpar,vparl, &
!$OMP inv_Bij,finterpol,dBdr_tmp,dBdtheta_tmp, &
!$OMP gradpar_J0Phi_tmp,gradpar_B_tmp,vExB_gradB_tmp, &
!$OMP vExB_gradr_1D,vExB_gradtheta_1D, &
!$OMP gradpar_B_1D,gradpar_J0Phi_1D, &
!$OMP rhsvpar,scoefvpar,vparg,tid) default(shared)
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
      
    Romp_maxalp_vpar(tid) = 0._RKIND
    scoefvpar             => Romp_scoefvpar(tid)%val
    vparg                 => Romp1_0Nvpar(tid)%val
    rhsvpar               => Romp2_0Nvpar(tid)%val
      
    vExB_gradr_1D     => Romp3_0Nvpar(tid)%val
    vExB_gradtheta_1D => Romp4_0Nvpar(tid)%val
    gradpar_B_1D      => Romp5_0Nvpar(tid)%val
    gradpar_J0Phi_1D  => Romp6_0Nvpar(tid)%val
      
    do ivpar = 0,Nvpar
       vparg(ivpar) = geom%vparg(ivpar)
    end do
      
!$OMP DO SCHEDULE(static)
    do iphi = 0,Nphi-1
      !*** vpar advection ***
      do itheta = jstart,jend
        do ir = f%istart_modif,f%iend_modif 
          inv_Bij      = 1._RKIND/init_magnet%B_norm(ir,itheta)
          dBdr_tmp     = init_magnet%dBdr(ir,itheta)
          dBdtheta_tmp = init_magnet%dBdtheta(ir,itheta)
          !-> computation of vExB contravariant components
          call compute_vExBgradr(geom,init_magnet,init_curr, &
            dJ0Phidtheta,dJ0Phidphi,ir,ir,itheta,itheta, &
            iphi,iphi,0,Nvpar,vExB_gradr_1D)
          call compute_vExBgradtheta(geom, &
            init_magnet,init_curr,dJ0Phidr,dJ0Phidphi,  &
            ir,ir,itheta,itheta,iphi,iphi,0,Nvpar, &
            vExB_gradtheta_1D)
          !-> computation of gradpar.B
          call compute_gradpar_B(geom,init_magnet,init_curr, &
            ir,ir,itheta,itheta,0,Nvpar,gradpar_B_1D)
          !-> computation of gradpar.Phi
          call compute_gradpar_J0Phi(geom, &
            init_magnet,init_curr,ir,ir,itheta,itheta, &
            iphi,iphi,0,Nvpar,gradpar_J0Phi_1D)
          !-> cubic spline computation 
          do ivpar = 0,Nvpar
            rhsvpar(ivpar) = f%values(ir,itheta,iphi,ivpar)
          end do
          call natural_omp_spline_coef(f%nspline1d_vpar(tid), &
            rhsvpar,f%BCvpar_left,f%BCvpar_right,scoefvpar)
          !*** splitting in vpar direction ***          
          do ivpar = 0,Nvpar
            vparl             = geom%vparg(ivpar)
            gradpar_B_tmp     = gradpar_B_1D(ivpar)
            gradpar_J0Phi_tmp = gradpar_J0Phi_1D(ivpar)
            !-> computation of vExB.gradB
            vExB_gradB_tmp = vExB_gradr_1D(ivpar)*dBdr_tmp + &
              vExB_gradtheta_1D(ivpar)*dBdtheta_tmp
            !-> computation of vpar displacement
            disp_vpar = dt*(- mu_value*gradpar_B_tmp - &
              real(Zi)*gradpar_J0Phi_tmp + &
              K_curv*vparl*inv_Bij*vExB_gradB_tmp)
            Romp_maxalp_vpar(tid) = max(Romp_maxalp_vpar(tid), &
              abs(disp_vpar/geom%dvpar))
            vstar       = geom%vparg(ivpar) - disp_vpar       
            vstar       = min(max(geom%vparg(0),vstar), &
              geom%vparg(Nvpar))   
            !-> f interpolation 
            call interpol1d_vpar(geom,scoefvpar, &
              vparg,vstar,finterpol)
            f%values(ir,itheta,iphi,ivpar) = finterpol
          end do
        end do
      end do
    end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
      
    maxalp_vpar = maxval(Romp_maxalp_vpar(1:Nbthread))
      
    call clck_time(eclock_advec1Dvpar)
    call clck_diff(bclock_advec1Dvpar,eclock_advec1Dvpar, &
      global_time_advec1Dvpar)
!R3 call r3_info_end (r3_info_index_0) !R3
  end subroutine advec1D_BSL_vpar
end module advec1D_BSL_module
