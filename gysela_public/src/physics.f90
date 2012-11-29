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
      
!----------------------------------------------------------------- 
! file : physics.f90
! date : 16/04/2010
!  Compute all the physics quantities 
!   (radial profiles) into two main subroutines, one dedicated
!   to the quantities associated to the guiding center
!   and the other one to the particles.
!   
! - Guiding-center quantities:
!    -> for the poloidal velocities
!      . nvpoloGCr_turb (due to ExB velocity)
!      . nvpoloGCr_neo (due to curvature)
!      . nvpoloGCr_vpar  (due to parallel velocity)
!    -> for toroidal angular momentum conservation
!      . LphiGCr         : toroidal angular momentum
!      . TpolGCr         : polarization term
!      . GammaGCr_turb   : guiding-center turbulent flux
!                          (due to ExB velocity)
!      . GammaGCr_neo    : guiding-center neoclassical flux
!      . RSphiGCr_turb   : toroidal Reynolds stress due to ExB
!      . RSphiGCr_neo    : neoclassical toroidal Reynolds stress
!    -> QGCr_turb        : guiding-center heat flux due to ExB
!    -> QGCr_neo         : guiding-center heat flux due to 
!                          curvature drift
!    -> niGCr            : guiding-center density
!    -> pressGCr         : guiding-center pressure
!    -> dpressGCr_dt     : time derivative of the guiding-center
!                          pressure: (P(tn+1)-P(tn))/dt
!    -> pressGCr_perp    : guiding-center perpendicular pressure
!    -> pressGCr_par     : guiding-center parallel pressure
!    -> dLphiGCr_dt      : time derivative of the toroidal angular
!                           momentum: (Lphi(tn+1)-Lphi(tn))/dt
!    -> QGCr_perp_turb   : guiding-center perpendicular heat flux
!                          due to ExB
!    -> QGCr_perp_neo    : guiding-center perpendicular neoclassical 
!                          heat flux
!    -> QGCr_par_turb    : guiding-center parallel heat flux due to ExB
!    -> QGCr_par_neo     : guiding-center parallel neoclassical heat flux 
!    -> QGCr_dtvpar      : guiding-center dtvpar term 
!                          (in parallel pressure balance)
!    -> dpGCr_dt_perp    : time derivative of the perpendicular pressure
!    -> dpGCr_dt_par     : time derivative of the parallel pressure
!    -> entropy    
!    -> L2 norm    
!    -> kinetic energy
!
! - Particles quantities: 
!    -> nir         : ion density
!    -> nuparr      : parallel velocity
!    -> pressr      : pressure
!    -> Gammar_turb : turbulent flux (due to ExB velocity)
!    -> Gammar_neo  : neoclassical flux (due to curvature drift)
!    -> Qr_turb     : turbulent heat flux (due to ExB velocity)
!    -> Qr_neo      : neoclassical heat flux 
!                     (due to curvature drift)
!    -> RSthetar    : poloidal Reynold stress
!    -> RSparr      : parallel Reynold stress
!    -> nvpolor_mag : poloidal velocity due to the magnetisation
!    -> potential energy
!
! Rk: All the previous radial profiles correspond to
!     flux surface average, i.e:
!       <.>_FS = \int . Js dtheta dphi /
!                \int Js dtheta dphi 
!      with Js the jacobian in space
!----------------------------------------------------------------- 
module physics_module
  use prec_const
  implicit none
      
  !-------------------------------------
  !  PRIVATE VARIABLES
  !------------------------------------- 
  !-> variables for the guiding-center and particle moments
  real(RKIND), dimension(:,:), pointer, private :: Fmoments_loc
  real(RKIND), dimension(:,:), pointer, private :: GC_Fmoments
  real(RKIND), dimension(:,:), pointer, private :: part_Fmoments
      
  !-------------------------------------
  !  PUBLIC VARIABLES
  !------------------------------------- 
  public
  !--> for guiding-center radial profiles
  real(RKIND), dimension(:), pointer, public :: nvpoloGCr_turb_diag
  real(RKIND), dimension(:), pointer, public :: nvpoloGCr_neo_diag
  real(RKIND), dimension(:), pointer, public :: nvpoloGCr_vpar_diag
  real(RKIND), dimension(:), pointer, public :: LphiGCr_diag
  real(RKIND), dimension(:), pointer, public :: TpolGCr_diag
  real(RKIND), dimension(:), pointer, public :: GammaGCr_turb_diag
  real(RKIND), dimension(:), pointer, public :: GammaGCr_neo_diag
  real(RKIND), dimension(:), pointer, public :: RSphiGCr_turb_diag
  real(RKIND), dimension(:), pointer, public :: RSphiGCr_neo_diag
  real(RKIND), dimension(:), pointer, public :: QGCr_turb_diag
  real(RKIND), dimension(:), pointer, public :: QGCr_neo_diag
  real(RKIND), dimension(:), pointer, public :: niGCr_diag
  real(RKIND), dimension(:), pointer, public :: pressGCr_diag
  real(RKIND), dimension(:), pointer, public :: pressGCr_perp_diag
  real(RKIND), dimension(:), pointer, public :: pressGCr_par_diag
  real(RKIND), dimension(:), pointer, public :: dpressGCr_dt_diag
  real(RKIND), dimension(:), pointer, public :: dLphiGCr_dt_diag
  real(RKIND), dimension(:), pointer, public :: QGCr_perp_turb_diag
  real(RKIND), dimension(:), pointer, public :: QGCr_perp_neo_diag
  real(RKIND), dimension(:), pointer, public :: QGCr_par_turb_diag
  real(RKIND), dimension(:), pointer, public :: QGCr_par_neo_diag
  real(RKIND), dimension(:), pointer, public :: QGCr_dtvpar_diag
  real(RKIND), dimension(:), pointer, public :: dpGCr_dt_perp_diag
  real(RKIND), dimension(:), pointer, public :: dpGCr_dt_par_diag
      
  !--> for particle radial profiles
  real(RKIND), dimension(:), pointer, public :: nir_diag
  real(RKIND), dimension(:), pointer, public :: nuparr_diag
  real(RKIND), dimension(:), pointer, public :: pressr_diag
  real(RKIND), dimension(:), pointer, public :: Gammar_turb_diag
  real(RKIND), dimension(:), pointer, public :: Gammar_neo_diag
  real(RKIND), dimension(:), pointer, public :: Qr_turb_diag
  real(RKIND), dimension(:), pointer, public :: Qr_neo_diag
  real(RKIND), dimension(:), pointer, public :: RSthetar_diag
  real(RKIND), dimension(:), pointer, public :: RSparr_diag
  real(RKIND), dimension(:), pointer, public :: nvpolor_mag_diag
      
  !--> for entropy and L2norm 
  real(RKIND), public :: entropy_diag
  real(RKIND), public :: L2norm_diag
  !--> for kinetic and potential energy
  real(RKIND), public :: Enkin_diag
  real(RKIND), public :: Enpot_diag
      
  include "Bstar_inline.h"
  include "velocities_inline.h"
  !******************************
  contains
  !******************************
#include "Bstar_inline.f90"
#include "velocities_inline.f90"
      
  !--------------------------------------------------------
  ! Constructor
  !--------------------------------------------------------
  subroutine new_physics(geom)
    use globals, only : nbmoments_GC, nbmoments_part
    use geometry_class
    use mem_alloc_module
    type(geometry), intent(in) :: geom
      
    integer :: max_nbmoments
      
    max_nbmoments = max(nbmoments_GC,nbmoments_part)
      
    !*** allocation for fluid moments computation ***
    call glob_allocate(Fmoments_loc,0,geom%Nr,1,max_nbmoments)
    call glob_allocate(GC_Fmoments,0,geom%Nr,1,nbmoments_GC)
    call glob_allocate(part_Fmoments,0,geom%Nr,1,nbmoments_part)
      
    !*** allocation of arrays for output saving of ***
    !***  guiding-center radial profiles           ***
    call glob_allocate(nvpoloGCr_turb_diag,0,geom%Nr, &
      'nvpoloGCr_turb_diag')
    call glob_allocate(nvpoloGCr_neo_diag,0,geom%Nr, &
      'nvpoloGCr_neo_diag')
    call glob_allocate(nvpoloGCr_vpar_diag,0,geom%Nr, &
      'nvpoloGCr_vpar_diag')
    call glob_allocate(LphiGCr_diag,0,geom%Nr,'LphiGCr_diag')
    call glob_allocate(TpolGCr_diag,0,geom%Nr,'TpolGCr_diag')
    call glob_allocate(GammaGCr_turb_diag,0,geom%Nr, &
      'GammaGCr_turb_diag')
    call glob_allocate(GammaGCr_neo_diag,0,geom%Nr, &
      'GammaGCr_neo_diag')
    call glob_allocate(RSphiGCr_turb_diag,0,geom%Nr, &
      'RSphiGCr_turb_diag')
    call glob_allocate(RSphiGCr_neo_diag,0,geom%Nr, &
      'RSphiGCr_neo_diag')
    call glob_allocate(QGCr_turb_diag,0,geom%Nr,'QGCr_turb_diag')
    call glob_allocate(QGCr_neo_diag,0,geom%Nr,'QGCr_neo_diag')
    call glob_allocate(niGCr_diag,0,geom%Nr,'niGCr_diag')
    call glob_allocate(pressGCr_diag,0,geom%Nr,'pressGCr_diag')
    call glob_allocate(pressGCr_perp_diag,0,geom%Nr, &
      'pressGCr_perp_diag')
    call glob_allocate(pressGCr_par_diag,0,geom%Nr, &
      'pressGCr_par_diag')
    call glob_allocate(dpressGCr_dt_diag,0,geom%Nr, &
      'dpressGCr_dt_diag')
    call glob_allocate(dLphiGCr_dt_diag,0,geom%Nr, &
      'dLphiGCr_dt_diag')
    call glob_allocate(QGCr_perp_turb_diag,0,geom%Nr, &
      'QGCr_perp_turb_diag')
    call glob_allocate(QGCr_perp_neo_diag,0,geom%Nr, &
      'QGCr_perp_neo_diag')
    call glob_allocate(QGCr_par_turb_diag,0,geom%Nr, &
      'QGCr_par_turb_diag')
    call glob_allocate(QGCr_par_neo_diag,0,geom%Nr, &
      'QGCr_par_neo_diag')
    call glob_allocate(QGCr_dtvpar_diag,0,geom%Nr,'QGCr_dtvpar_diag')
    call glob_allocate(dpGCr_dt_par_diag,0,geom%Nr, &
      'dpGCr_dt_par_diag')
    call glob_allocate(dpGCr_dt_perp_diag,0,geom%Nr, &
      'dpGCr_dt_perp_diag')
      
    !*** allocation of arrays for output saving of ***
    !***  particle radial profiles                 ***
    call glob_allocate(nir_diag,0,geom%Nr,'nir_diag')
    call glob_allocate(nuparr_diag,0,geom%Nr,'nuparr_diag')
    call glob_allocate(pressr_diag,0,geom%Nr,'pressr_diag')
    call glob_allocate(Gammar_turb_diag,0,geom%Nr, &
      'Gammar_turb_diag')
    call glob_allocate(Gammar_neo_diag,0,geom%Nr,'Gammar_neo_diag')
    call glob_allocate(Qr_turb_diag,0,geom%Nr,'Qr_turb_diag')
    call glob_allocate(Qr_neo_diag,0,geom%Nr,'Qr_neo_diag')
    call glob_allocate(RSthetar_diag,0,geom%Nr,'RSthetar_diag')
    call glob_allocate(RSparr_diag,0,geom%Nr,'RSparr_diag')
    call glob_allocate(nvpolor_mag_diag,0,geom%Nr, &
      'nvpolor_mag_diag')
  end subroutine new_physics
      
  !--------------------------------------------------------
  ! Destructor
  !--------------------------------------------------------
  subroutine delete_physics()
    use mem_alloc_module
      
    !*** allocation for fluid moments computation ***
    call glob_deallocate(Fmoments_loc)
    call glob_deallocate(GC_Fmoments)
    call glob_deallocate(part_Fmoments)
      
    !*** deallocation of arrays for output saving of ***
    !***  guiding-center radial profiles             ***
    call glob_deallocate(nvpoloGCr_turb_diag)
    call glob_deallocate(nvpoloGCr_neo_diag)
    call glob_deallocate(nvpoloGCr_vpar_diag)
    call glob_deallocate(LphiGCr_diag)
    call glob_deallocate(TpolGCr_diag)
    call glob_deallocate(GammaGCr_turb_diag)
    call glob_deallocate(GammaGCr_neo_diag)
    call glob_deallocate(RSphiGCr_turb_diag)
    call glob_deallocate(RSphiGCr_neo_diag)
    call glob_deallocate(QGCr_turb_diag)
    call glob_deallocate(QGCr_neo_diag)
    call glob_deallocate(niGCr_diag)
    call glob_deallocate(pressGCr_diag)
    call glob_deallocate(pressGCr_perp_diag)
    call glob_deallocate(pressGCr_par_diag)
    call glob_deallocate(dpressGCr_dt_diag)
    call glob_deallocate(dLphiGCr_dt_diag)
    call glob_deallocate(QGCr_perp_turb_diag)
    call glob_deallocate(QGCr_perp_neo_diag)
    call glob_deallocate(QGCr_par_turb_diag)
    call glob_deallocate(QGCr_par_neo_diag)
    call glob_deallocate(QGCr_dtvpar_diag)
    call glob_deallocate(dpGCr_dt_par_diag)
    call glob_deallocate(dpGCr_dt_perp_diag)
      
    !*** deallocation of arrays for output saving of ***
    !***  particle radial profiles                   ***
    call glob_deallocate(nir_diag)
    call glob_deallocate(nuparr_diag)
    call glob_deallocate(pressr_diag)
    call glob_deallocate(Gammar_turb_diag)
    call glob_deallocate(Gammar_neo_diag)
    call glob_deallocate(Qr_turb_diag)
    call glob_deallocate(Qr_neo_diag)
    call glob_deallocate(RSthetar_diag)
    call glob_deallocate(RSparr_diag)
    call glob_deallocate(nvpolor_mag_diag)
  end subroutine delete_physics
      
  !----------------------------------------------------------------
  ! Computation of the fluid moments associated to the 
  !  guiding-centers, corresponding to 
  !   1) nvpoloGC_turb(r)        = <\int f*vE_GC.etheta d3v>_FS
  !   2) nvpoloGC_neo(r)         = <\int f*vD.etheta d3v>_FS
  !   3) nvpoloGC_vpar(r)        = <\int f*vpar bstar.etheta d3v>_FS
  ! where vE_GC is the guiding center ExB velocity
  !   4) LphiGC(r)               = <\int f*vpar R^2*B*gradphi/B d3v>_FS
  !   5) TpolGC(r)               = <\int f*dJ0Phi/dphi d3v>_FS
  !   6) GammaGC_turb(r)         = <\int f*vE_GC.gradr d3v>_FS
  !   7) GammaGC_neo(r)          = <\int f*vD_GC.gradr d3v>_FS
  !   8) RSphiGC_turb(r)         = <\int f*vE_GC.gradr*vpar* 
  !                                R^2*B*gradphi/B d3v>_FS
  !   9) RSphiGC_neo(r)          = <\int f*vD_GC.gradr*vpar*
  !                                R^2*B*gradphi/B d3v>_FS
  !   10) QGC_turb(r)            = <\int f*E*vE_GC.gradr d3v>_FS
  !   11) QGC_neo(r)             = <\int f*E*vD_GC.gradr d3v>_FS 
  !   12) niGC(r)                = <\int f d3v>_FS 
  !   13) pressGC(r)             = <\int f*E d3v>_FS 
  !   14) pressGC_perp(r)        = <\int f*mu*B d3v>_FS 
  !   15) pressGC_par(r)         = <\int f*v_par^2 d3v>_FS 
  !   16) dpressGC_dt(r)         = (pressGC(tn+1)-pressGC(tn))/dt
  !   17) dLphiGC_dt(r)          = (LphiGC(tn+1)-LphiGC(tn))/dt
  !   18) QGCr_perp_turb(r)      = <\int f*mu*B*(vE_GC+vD_GC).gradr d3v>_FS
  !   19) QGCr_perp_neo(r)       = <\int f*mu*B*(vE_GC+vD_GC).gradr d3v>_FS
  !   20) QGCr_par_turb(r)       = <\int f*0.5*vpar*vpar*
  !                               (vE_GC+vD_GC).gradr d3v>_FS
  !   21) QGCr_par_neo(r)        = <\int f*0.5*vpar*vpar*
  !                               (vE_GC+vD_GC).gradr d3v>_FS
  !   22) QGCr_dtvpar(r)         = <\int f*v_par*dv_par/dt>
  !   23) dpGCr_dt_perp(r)       = (pressGC_perp(tn+1)-pressGC_perp(tn))/dt
  !   24) dpGCr_dt_par(r)        = (pressGC_par(tn+1)-pressGC_par(tn))/dt
  !   25) entropy    
  !   26) L2 norm    
  !   27) kinetic energy
  !
  ! Rk1: <.>_FS corresponds to the flux surface average,
  !      i.e <.>_FS = \int . Js dtheta dphi /
  !                   \int Js dtheta dphi 
  !      with Js the jacobian in space
  !
  ! Rk2: the poloidal velocity nvpolo corresponds to the
  !     sum of the three previous moments
  !     nvpolo(r) = nvpoloGC_turb(r) + nvpoloGC_neo(r) + 
  !                 nvpoloGC_vpar(r)
  !       NB: vD contains the magnetisation velocity
  !----------------------------------------------------------------
  subroutine compute_GC_moments(geom,coord_sys, &
    init_prof,init_magnet,init_curr, &
    fn,fnp1,Sfmu_eq,SdJ0Phidr,SdJ0Phidtheta, &
    SdJ0Phidphi,dt_diag,iter_time)
    use globals, only : nbmoments_GC, &
      istart, iend, jstart, jend, mu_id, &
      Nbproc_tot, Nbthread, Rarray1_Nr, Rarray2_Nr, &
      diag_targ, ipow, &
      reduceddt, reducedbegin, reducedend, K_curv, Zi
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp1_0Nr, Romp2_0Nr, &
      Romp3_0Nr, Romp4_0Nr, Romp5_0Nr, Romp6_0Nr, &
      Romp1_1nbmom_0Nr
    use clock_module
    use init_profile_class
    use init_magnetic_class
    use init_current_class
    use fdistribu5d_class
    implicit none
    include "mpiperso.h"    
    type(geometry)               , intent(in) :: geom
    type(coord_system)           , intent(in) :: coord_sys
    type(init_profile)           , intent(in) :: init_prof
    type(init_magnetic)          , intent(in) :: init_magnet
    type(init_current)           , intent(in) :: init_curr
    type(fdistribu5d)            , intent(in) :: fn
    type(fdistribu5d)            , intent(in) :: fnp1
    real(RKIND), dimension(:,:,:), pointer    :: Sfmu_eq
    real(RKIND), dimension(:,:,:), pointer    :: SdJ0Phidr 
    real(RKIND), dimension(:,:,:), pointer    :: SdJ0Phidtheta 
    real(RKIND), dimension(:,:,:), pointer    :: SdJ0Phidphi
    real(RKIND)                  , intent(in) :: dt_diag
    real(RKIND)                  , intent(in) :: iter_time
      
    integer     :: imom, tid, ithread, ierr
    integer     :: ir, itheta, iphi, ivpar, imu
    integer     :: nbmoments, nbelements, nbadapt
    real(RKIND) :: sqrt_g22, vpar, mu_value
    real(RKIND) :: Bij, energy, fnval, fnp1val
    !-> variables for the poloidal velocities
    real(RKIND) :: bstar_gradtheta_tmp
    real(RKIND) :: vpar_etheta, vExB_etheta, vD_etheta
    real(RKIND), dimension(:), pointer :: vExB_gradtheta_1D
    real(RKIND), dimension(:), pointer :: vD_gradtheta_1D
    !-> variables for entropy and L2norm
    real(RKIND) :: entropy, L2norm
    real(RKIND), dimension(1:Nbthread) :: entropy_thread 
    real(RKIND), dimension(1:Nbthread) :: L2norm_thread
    !-> variables for the kinetic energy
    real(RKIND) :: fmfeq_E, Enkin
    real(RKIND), dimension(1:Nbthread) :: Enkin_thread 
    !-> variables for toroidal angular momentum
    real(RKIND) :: Rij, inv_Bij, B_gphi
    real(RKIND) :: vExB_gradr, vD_gradr
    real(RKIND), dimension(:), pointer :: vExB_gradr_1D
    real(RKIND), dimension(:), pointer :: vD_gradr_1D
    !-> variables for parallel and perp pressure
    real(RKIND) :: gradpar_B_tmp, gradpar_J0Phi_tmp
    real(RKIND) :: dBdr_tmp, dBdtheta_tmp, dvpar
    real(RKIND) :: vExB_gradB_tmp
    real(RKIND), dimension(:), pointer :: gradpar_B_1D
    real(RKIND), dimension(:), pointer :: gradpar_J0Phi_1D
      
    !-> variables for the flux surface average
    real(RKIND) :: jacob_space_tmp, jacob_vel_tmp
    real(RKIND) :: coeff_intdr_tmp
    real(RKIND) :: coeff_intdtheta_tmp, coeff_intdphi_tmp
    real(RKIND) :: coeff_intdvpar_tmp, coeff_intdmu_tmp
    real(RKIND) :: coeff_int_volume
    real(RKIND) :: intdthetadphi_Js_tmp
    real(RKIND) :: coeff_pressr
      
    ! -> variables for time perfomance analysis
    integer(TIMEPREC) :: tdeb, tfin
    real(RKIND)       :: tcomm, tcomp1, tcomp2
      
    !-> BECAREFUL: if you want to add the computation of
    !    a guiding-center moment, do not forget to change
    !    nbmoments_GC in 'globals.f90' file
    nbmoments        = nbmoments_GC
    mu_value         = fnp1%mu_value
    coeff_intdmu_tmp = geom%coeff_intdmu(mu_id)
      
    !*** initialisation to 0 of the radial arrays ***
    do imom = 1,nbmoments
      do ir =  0,geom%Nr
        Fmoments_loc(ir,imom) = 0._RKIND
        GC_Fmoments(ir,imom)  = 0._RKIND
      end do
    end do
      
    call clck_time(tdeb)
#ifdef _OPENMP
!$OMP PARALLEL private(tid,ir,itheta,iphi,ivpar,imom, &
!$OMP Bij,sqrt_g22,vpar,energy,fnval,fnp1val,fmfeq_E, &
!$OMP dBdr_tmp, dBdtheta_tmp, gradpar_B_tmp, gradpar_J0Phi_tmp, &
!$OMP coeff_intdr_tmp,coeff_intdtheta_tmp,&
!$OMP coeff_intdphi_tmp,coeff_intdvpar_tmp, &
!$OMP coeff_int_volume,jacob_space_tmp,jacob_vel_tmp, &
!$OMP bstar_gradtheta_tmp,vExB_gradtheta_1D,vD_gradtheta_1D, &
!$OMP gradpar_J0Phi_1D, gradpar_B_1D, &
!$OMP inv_Bij, Rij, B_gphi, vExB_gradr, vD_gradr, &
!$OMP vExB_gradr_1D, vD_gradr_1D, vExB_gradB_tmp, &
!$OMP vpar_etheta,vD_etheta,vExB_etheta,dvpar) &
!$OMP default(shared)
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
      
    vExB_gradtheta_1D => Romp1_0Nr(tid)%val
    vD_gradtheta_1D   => Romp2_0Nr(tid)%val
    vExB_gradr_1D     => Romp3_0Nr(tid)%val
    vD_gradr_1D       => Romp4_0Nr(tid)%val
    gradpar_B_1D      => Romp5_0Nr(tid)%val
    gradpar_J0Phi_1D  => Romp6_0Nr(tid)%val
      
    !*** initialisation of the values to 0 ***
    do imom = 1,nbmoments
      do ir = 0,geom%Nr
        Romp1_1nbmom_0Nr(tid)%val(imom,ir) = 0._RKIND
      end do
    end do
    entropy_thread(tid) = 0._RKIND
    L2norm_thread(tid)  = 0._RKIND
    Enkin_thread(tid)   = 0._RKIND
      
!$OMP DO SCHEDULE(STATIC)
    do ivpar = 0,geom%Nvpar
      vpar               = geom%vparg(ivpar)
      coeff_intdvpar_tmp = geom%coeff_intdvpar(ivpar)
      do iphi = 0,geom%Nphi-1
        coeff_intdphi_tmp = geom%coeff_intdphi(iphi)
        do itheta = jstart,jend
          coeff_intdtheta_tmp = geom%coeff_intdtheta(itheta)
          !*** computation of the poloidal velocities    ***
          !-> computation of vExB_gradtheta(0:Nr)
          call compute_vExBgradtheta(geom,init_magnet,init_curr, &
            SdJ0Phidr,SdJ0Phidphi,istart,iend,itheta,itheta, &
            iphi,iphi,ivpar,ivpar,vExB_gradtheta_1D)
          !-> computation of vD_gradtheta(0:Nr)
          call compute_vDgradtheta(geom,init_magnet,init_curr, &
            istart,iend,itheta,itheta,ivpar,ivpar,vD_gradtheta_1D)
          !-> computation of vExB_gradr(0:Nr)
          call compute_vExBgradr(geom,init_magnet,init_curr, &
            SdJ0Phidtheta,SdJ0Phidphi,istart,iend,itheta,itheta, &
            iphi,iphi,ivpar,ivpar,vExB_gradr_1D)
          !-> computation of vD_gradr(0:Nr)
          call compute_vDgradr(geom,init_magnet,init_curr, &
            istart,iend,itheta,itheta,ivpar,ivpar,vD_gradr_1D)
          !-> computation of gradpar.B(0:Nr)
          call compute_gradpar_B(geom,init_magnet,init_curr, &
            istart,iend,itheta,itheta,ivpar,ivpar,gradpar_B_1D)
          !-> computation of gradpar.Phi(0:Nr)
          call compute_gradpar_J0Phi(geom,init_magnet,init_curr, &
            istart,iend,itheta,itheta,iphi,iphi,ivpar,ivpar, &
            gradpar_J0Phi_1D)
          do ir = istart,iend
            Bij          = init_magnet%B_norm(ir,itheta)
            dBdr_tmp     = init_magnet%dBdr(ir,itheta)
            dBdtheta_tmp = init_magnet%dBdtheta(ir,itheta)
            energy       = 0.5_RKIND*vpar*vpar+mu_value*Bij
            fnval        = fn%values(ir,itheta,iphi,ivpar)
            fnp1val      = fnp1%values(ir,itheta,iphi,ivpar)
            fmfeq_E      = (fnp1val-Sfmu_eq(ir,itheta,ivpar))*energy
            !-> jacobian_velocity = 2*pi*Bstar
            call compute_jacobian_velocity(geom,init_magnet, &
              init_curr,ir,itheta,ivpar,jacob_vel_tmp)            
            !-> jacobian in space
            jacob_space_tmp = jacobian_space(ir,itheta)
            !-> compute the elements of volume in space and velocity
            coeff_int_volume = jacob_space_tmp * &
              coeff_intdtheta_tmp*coeff_intdphi_tmp * &
              jacob_vel_tmp*coeff_intdvpar_tmp*coeff_intdmu_tmp 
            !-> compute bstar.gradtheta
#ifdef NOPRECOMPUTE 
            call compute_bstar_contravariant(geom, &
              init_magnet,init_curr,ir,itheta,ivpar, &
              Sbstar_gradx2=bstar_gradtheta_tmp)
#else
            call precomputed_bstar_gradtheta(ir,itheta,ivpar, &
              bstar_gradtheta_tmp)
#endif
            !*** computation of the projection of the      ***
            !***  poloidal velocities on the unit poloidal ***
            !***  vector                                   ***
            sqrt_g22    = sqrt(coord_sys%g22(ir,itheta))
            !-> vpar.etheta = vpar*r*bstar_gradtheta 
            vpar_etheta = vpar*sqrt_g22*bstar_gradtheta_tmp
            !-> vD.etheta = r*vD_gradtheta 
            vD_etheta   = sqrt_g22*vD_gradtheta_1D(ir-istart)
            !-> vE.etheta = r*vExB_gradtheta
            vExB_etheta = sqrt_g22*vExB_gradtheta_1D(ir-istart)
            !*** Computations for toroidal angular momentum ***
            inv_Bij     = 1._RKIND/Bij
            B_gphi      = init_magnet%B_gradphi(ir,itheta)
            Rij         = R(ir,itheta)
            !-> vD.gradr
            vD_gradr    = vD_gradr_1D(ir-istart)
            !-> vExB.gradr
            vExB_gradr  = vExB_gradr_1D(ir-istart)
            !-> computations for dvpar
            gradpar_B_tmp     = gradpar_B_1D(ir-istart)
            gradpar_J0Phi_tmp = gradpar_J0Phi_1D(ir-istart)
            !-> computation of vExB.gradB
            vExB_gradB_tmp = vExB_gradr_1D(ir-istart)*dBdr_tmp + &
              vExB_gradtheta_1D(ir-istart)*dBdtheta_tmp
            !-> computation of vpar displacement
            dvpar = - mu_value*gradpar_B_tmp - &
              real(Zi)*gradpar_J0Phi_tmp + &
              K_curv*vpar*inv_Bij*vExB_gradB_tmp
            !*** computation of the integrals in velocity and ***
            !***  in theta and phi                            ***
            !-> compute nvpoloGCr_turb
            Romp1_1nbmom_0Nr(tid)%val(1,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(1,ir) + &
              fnp1val*vExB_etheta*coeff_int_volume
            !-> compute nvpoloGCr_neo
            Romp1_1nbmom_0Nr(tid)%val(2,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(2,ir) + &
              fnp1val*vD_etheta*coeff_int_volume
            !-> compute nvpoloGCr_vpar
            Romp1_1nbmom_0Nr(tid)%val(3,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(3,ir) + &
              fnp1val*vpar_etheta*coeff_int_volume
            !-> compute LphiGCr
            Romp1_1nbmom_0Nr(tid)%val(4,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(4,ir) + &
              fnp1val*vpar*Rij*Rij*B_gphi*inv_Bij * &
              coeff_int_volume
            !-> compute TpolGCr
            Romp1_1nbmom_0Nr(tid)%val(5,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(5,ir) + &
              fnp1val* SdJ0Phidphi(ir,itheta,iphi) * &
              coeff_int_volume
            !-> compute GammaGCr_turb
            Romp1_1nbmom_0Nr(tid)%val(6,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(6,ir) + &
              fnp1val*vExB_gradr*coeff_int_volume
            !-> compute GammaGCr_neo
            Romp1_1nbmom_0Nr(tid)%val(7,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(7,ir) + &
              fnp1val*vD_gradr*coeff_int_volume
            !-> compute RSphiGCr_turb
            Romp1_1nbmom_0Nr(tid)%val(8,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(8,ir) + &
              fnp1val*vpar*Rij*Rij*B_gphi*inv_Bij * &
              vExB_gradr*coeff_int_volume
            !-> compute RSphiGCr_neo
            Romp1_1nbmom_0Nr(tid)%val(9,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(9,ir) + &
              fnp1val*vpar*Rij*Rij*B_gphi*inv_Bij * &
              vD_gradr*coeff_int_volume
            !-> compute QGCr_turb
            Romp1_1nbmom_0Nr(tid)%val(10,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(10,ir) + &
              fnp1val*vExB_gradr*energy* &
              coeff_int_volume
            !-> compute QGCr_neo
            Romp1_1nbmom_0Nr(tid)%val(11,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(11,ir) + &
              fnp1val*vD_gradr*energy* &
              coeff_int_volume
            !-> compute niGCr
            Romp1_1nbmom_0Nr(tid)%val(12,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(12,ir) + &
              fnp1val*coeff_int_volume
            !-> compute pressGCr_tnp1
            Romp1_1nbmom_0Nr(tid)%val(13,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(13,ir) + &
              fnp1val*energy*coeff_int_volume
            !-> compute pressGCr_perp
            Romp1_1nbmom_0Nr(tid)%val(14,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(14,ir) + &
              fnp1val*mu_value*Bij*coeff_int_volume
            !-> compute pressGCr_par
            Romp1_1nbmom_0Nr(tid)%val(15,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(15,ir) + &
              fnp1val*vpar*vpar*coeff_int_volume
            !-> compute pressGCr_tn (for dpressGCr_dt)
            Romp1_1nbmom_0Nr(tid)%val(16,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(16,ir) + &
              fnval*energy*coeff_int_volume
            !-> compute LphiGCr_tn (for dLphiGCr_dt)
            Romp1_1nbmom_0Nr(tid)%val(17,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(17,ir) + &
              fnval*vpar*Rij*Rij*B_gphi*inv_Bij * &
              coeff_int_volume
            !-> compute QGCr_perp_turb
            Romp1_1nbmom_0Nr(tid)%val(18,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(18,ir) + &
              fnp1val*vExB_gradr*mu_value*Bij* &
              coeff_int_volume
            !-> compute QGCr_neo_turb
            Romp1_1nbmom_0Nr(tid)%val(19,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(19,ir) + &
              fnp1val*vD_gradr*mu_value*Bij* &
              coeff_int_volume
            !-> compute QGCr_par_turb
            Romp1_1nbmom_0Nr(tid)%val(20,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(20,ir) + &
              fnp1val*vExB_gradr*vpar*vpar* &
              0.5_RKIND*coeff_int_volume
            !-> compute QGCr_par_neo
            Romp1_1nbmom_0Nr(tid)%val(21,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(21,ir) + &
              fnp1val*vD_gradr*vpar*vpar* &
              0.5_RKIND*coeff_int_volume
            !-> compute Par_press_dtvpar
            Romp1_1nbmom_0Nr(tid)%val(22,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(22,ir) + &
              fnp1val*vpar*dvpar* &
              coeff_int_volume
            !-> compute press_perp_tn (for dpGCr_dt_perp)
            Romp1_1nbmom_0Nr(tid)%val(23,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(23,ir) + &
              fnval*mu_value*Bij*coeff_int_volume
            !-> compute press_par_tn (for dpGCr_dt_par)
            Romp1_1nbmom_0Nr(tid)%val(24,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(24,ir) + &
              fnval*vpar*vpar*coeff_int_volume
            !-> compute entropy
            coeff_intdr_tmp = geom%coeff_intdr(ir)
            if (fnp1val .ne. 0._RKIND) then
              entropy_thread(tid) = entropy_thread(tid) - &
                fnp1val*log(abs(fnp1val))*coeff_int_volume * &
                coeff_intdr_tmp
            end if
            !-> compute L2norm
            L2norm_thread(tid)  = L2norm_thread(tid)  + &
              fnp1val*fnp1val*coeff_int_volume * &
              coeff_intdr_tmp
            !-> compute Enkin
            Enkin_thread(tid)  = Enkin_thread(tid)  + &
              fmfeq_E*coeff_int_volume * &
              coeff_intdr_tmp
          end do ! end of do ir
        end do   ! end of do itheta
      end do     ! end of do iphi
    end do       ! end of do ivpar
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
    !*** fill an unique array containing all the       ***
    !***  radial profiles of the fluid moments         ***
    !*** Rq: each processor fills the part istart:iend ***
    !***     and all the others values are put equal   ***
    !***     to 0                                      ***
    entropy = 0._RKIND
    L2norm  = 0._RKIND
    Enkin   = 0._RKIND
    do ithread = 1,Nbthread
      do ir = istart,iend
        do imom = 1,nbmoments
          Fmoments_loc(ir,imom) = Fmoments_loc(ir,imom) + & 
            Romp1_1nbmom_0Nr(ithread)%val(imom,ir)
        end do
      end do
      entropy = entropy + entropy_thread(ithread)
      L2norm  = L2norm + L2norm_thread(ithread) 
      Enkin   = Enkin + Enkin_thread(ithread) 
    end do
      
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcomp1)
      
    !*** MPI_REDUCE ****
    call clck_time(tdeb)
    if (Nbproc_tot.ne.1) then
      nbelements = nbmoments*(geom%Nr+1)
      call MPI_REDUCE(Fmoments_loc,GC_Fmoments, &
        nbelements,MPI_REAL8,MPI_SUM,diag_targ(2), &
        MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(entropy,entropy_diag,1, &
        MPI_REAL8,MPI_SUM,diag_targ(1),MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(L2norm,L2norm_diag,1, &
        MPI_REAL8,MPI_SUM,diag_targ(1),MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(Enkin,Enkin_diag,1, &
        MPI_REAL8,MPI_SUM,diag_targ(1),MPI_COMM_WORLD,ierr)
    else
      do imom = 1,nbmoments
        do ir = 0,geom%Nr
          GC_Fmoments(ir,imom) = Fmoments_loc(ir,imom)
        end do
      end do
      entropy_diag = entropy
      L2norm_diag  = L2norm
      Enkin_diag   = Enkin
    end if
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcomm)
      
    !*** computation of the flux surface average ***
    call clck_time(tdeb)
    !-> for the computation of the pressure and its time derivative
    coeff_pressr = 1._RKIND/(0.5_RKIND+float(ipow))
    if ((reduceddt .ne. 0._RKIND) .and. &
      (iter_time .ge. reducedbegin) .and. &
      (iter_time .le. reducedend)) then
      nbadapt    = MAX(1, nint(dt_diag / reduceddt))
    else
      nbadapt    = 1
    end if
    do ir = 0,geom%Nr
      intdthetadphi_Js_tmp    = intdthetadphi_Js(ir)
      nvpoloGCr_turb_diag(ir) = &
        GC_Fmoments(ir,1)/intdthetadphi_Js_tmp
      nvpoloGCr_neo_diag(ir)  = &
        GC_Fmoments(ir,2)/intdthetadphi_Js_tmp
      nvpoloGCr_vpar_diag(ir) = &
        GC_Fmoments(ir,3)/intdthetadphi_Js_tmp
      LphiGCr_diag(ir)        = &
        GC_Fmoments(ir,4)/intdthetadphi_Js_tmp
      TpolGCr_diag(ir)        = &
        GC_Fmoments(ir,5)/intdthetadphi_Js_tmp
      GammaGCr_turb_diag(ir)  = &
        GC_Fmoments(ir,6)/intdthetadphi_Js_tmp
      GammaGCr_neo_diag(ir)   = &
        GC_Fmoments(ir,7)/intdthetadphi_Js_tmp
      RSphiGCr_turb_diag(ir)  = &
        GC_Fmoments(ir,8)/intdthetadphi_Js_tmp
      RSphiGCr_neo_diag(ir)   = &
        GC_Fmoments(ir,9)/intdthetadphi_Js_tmp
      QGCr_turb_diag(ir)      = &
        GC_Fmoments(ir,10)/intdthetadphi_Js_tmp
      QGCr_neo_diag(ir)       = &
        GC_Fmoments(ir,11)/intdthetadphi_Js_tmp
      niGCr_diag(ir)          = &
        GC_Fmoments(ir,12)/intdthetadphi_Js_tmp
      pressGCr_diag(ir)       = coeff_pressr* &
        GC_Fmoments(ir,13)/intdthetadphi_Js_tmp
      pressGCr_perp_diag(ir)  = &
        GC_Fmoments(ir,14)/intdthetadphi_Js_tmp
      pressGCr_par_diag(ir)   = &
        GC_Fmoments(ir,15)/intdthetadphi_Js_tmp
      dpressGCr_dt_diag(ir)   = coeff_pressr* &
        (GC_Fmoments(ir,13)-GC_Fmoments(ir,16))/ &
        (intdthetadphi_Js_tmp*dt_diag) *nbadapt
      dLphiGCr_dt_diag(ir)    = &
        (GC_Fmoments(ir,4)-GC_Fmoments(ir,17))/ &
        (intdthetadphi_Js_tmp*dt_diag) *nbadapt
      QGCr_perp_turb_diag(ir) = &
        GC_Fmoments(ir,18)/intdthetadphi_Js_tmp
      QGCr_perp_neo_diag(ir) = &
        GC_Fmoments(ir,19)/intdthetadphi_Js_tmp
      QGCr_par_turb_diag(ir)       = &
        GC_Fmoments(ir,20)/intdthetadphi_Js_tmp
      QGCr_par_neo_diag(ir)       = &
        GC_Fmoments(ir,21)/intdthetadphi_Js_tmp
      QGCr_dtvpar_diag(ir)    = &
        GC_Fmoments(ir,22)/intdthetadphi_Js_tmp
      dpGCr_dt_perp_diag(ir)    = &
        (GC_Fmoments(ir,14)-GC_Fmoments(ir,23))/ &
        (intdthetadphi_Js_tmp*dt_diag) *nbadapt
      dpGCr_dt_par_diag(ir)    = &
        (GC_Fmoments(ir,15)-GC_Fmoments(ir,24))/ &
        (intdthetadphi_Js_tmp*dt_diag) *nbadapt
    end do
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcomp2)
      
#ifdef TIMER    
    write(6,'(I4,A,F13.7,A,F13.7,A,F13.7)') pglobal_id, &
      " Temps GC moments: comp ", tcomp1," comm ", tcomm, &
      " comp2 ", tcomp2
#endif
  end subroutine compute_GC_moments
      
  !----------------------------------------------------------------
  ! Computation of the fluid moments associated to the 
  !  particles, corresponding to 
  !   1) Ni(r)           = <\int J0.f d3v>_FS
  !   2) nVpar(r)        = <\int J0.f*vpar d3v>_FS
  !   3) P(r)            = <\int J0.f*E d3v>_FS
  !   4) Gamma_turb(r)   = <\int J0.f*vE.er d3v>_FS
  !   5) Gamma_neo(r)    = <\int J0.f*vD.er d3v>_FS
  !   6) Q_turb(r)       = <\int J0.f*E*vE.er d3v>_FS
  !   7) Q_neo(r)        = <\int J0.f*E*vD.er d3v>_FS 
  !   8) RStheta(r)      = <\int J0.(f*vExB.gradr) d3v>_FS
  !   9) RSpar(r)        = <\int J0.f vtilde_r vpar d3v>_FS 
  !   10) nvpolor_mag(r) = -<rot(M).etheta>_FS
  !   11) potential energy
  !
  ! Rk: <.>_FS corresponds to the flux surface average,
  !      i.e <.>_FS = \int . Js dtheta dphi /
  !                   \int Js dtheta dphi 
  !      with Js the jacobian in space
  !
  ! Rk2: Compute specifically the magnetisation component in nvpolo
  !    -> nvpolor_mag = nvpolo - nvpolor_GuidingCentre
  !  . nvpolor_mag = -rot(M).etheta                   
  !               = -r*rot(M).gradtheta              
  !     with rot(M).gradtheta = -1/Js d/dr(M_phi)     
  !      where Js is the jacobian in space and        
  !      M_phi is the toroidal covariant component of 
  !       M = -\int (\int J0.f vec_b Jv dvpar)*mu dmu 
  !     i.e M_phi = - b_phi *                         
  !          \int (\int J0.f vec_b Jv dvpar)*mu dmu     
  !----------------------------------------------------------------
  subroutine compute_part_moments(geom,coord_sys, &
    init_prof,init_magnet,init_curr,J0,fnp1,sJ0ftmp,sJ0fmu_eq,sPhi)
    use globals, only : nbmoments_part, &
      istart, iend, jstart, jend, &
      kstart, kend, lstart, lend, mu_id, Nbproc_tot, &
      Nbthread, ipow, Zi, Rarray1_Nr, Rarray1_PNrPNthetaNphi, &
      Rarray2_PNrPNthetaNphi, Rarray3_PNrPNthetaNphi, &
      diag_targ, Rarray1_NrNthetaNphi, Rarray2_NrNthetaNphi
#ifdef _OPENMP
    use OMPutils_module, only : omp_get_thread_num
#endif
    use OMPutils_module, only : Romp1_1Nrp1_1Nthetap1, &
      Romp2_1Nrp1_1Nthetap1, Comp1_1Nrp1_1Ntheta, Romp1_0Nr, &
      Romp2_0Nr, Romp3_0Nr, Romp4_0Nr, Romp5_0Nr, &
      Romp1_1nbmom_0Nr, omp1dvector
    use clock_module
    use coord_system_class
    use init_profile_class
    use init_magnetic_class
    use init_current_class
    use fdistribu5d_class
    use gyroaverage_class
    use efield_module, only : compute_dPhidr_loc, &
      compute_dPhidtheta_loc, compute_dPhidphi_loc
    use utils_module , only : deriv1
    implicit none
    include "mpiperso.h"    
    type(geometry)                  , intent(in)    :: geom
    type(coord_system)              , intent(in)    :: coord_sys
    type(init_profile)              , intent(in)    :: init_prof
    type(init_magnetic)             , intent(in)    :: init_magnet
    type(init_current)              , intent(in)    :: init_curr
    type(J0operator)                , intent(inout) :: J0
    type(fdistribu5d)               , intent(in)    :: fnp1
    type(fdistribu5d)               , intent(inout) :: sJ0ftmp
    real(RKIND), dimension(:,:,:)   , pointer       :: sJ0fmu_eq
    real(RKIND), dimension(0:,0:,0:), intent(in)    :: sPhi
      
    integer     :: imom, tid, ithread, ierr
    integer     :: ir, itheta, iphi, ivpar
    integer     :: nbmoments, nbelements
    real(RKIND) :: sqrt_g22, vpar, mu_value
    real(RKIND) :: Bij, energy
    real(RKIND) :: coeff_pressr
    real(RKIND) :: J0f, J0fmu_eq, J0fv, J0fE, fval
    !-> variables for the flux surface average
    real(RKIND) :: jacob_space_tmp, jacob_vel_tmp
    real(RKIND) :: coeff_intdr_tmp, coeff_intdtheta_tmp
    real(RKIND) :: coeff_intdphi_tmp, coeff_intdvpar_tmp
    real(RKIND) :: coeff_intdmu_tmp
    real(RKIND) :: coeff_int_velocities, coeff_int_volume
    real(RKIND) :: intdthetadphi_Js_tmp
    !-> variables for the Reynold Stress
    real(RKIND) :: RStheta_tmp
    type(omp1dvector), dimension(:), pointer :: RStheta
    !-> variables for the velocities
    real(RKIND) :: vpar_er, vExB_er, vD_er
    real(RKIND) :: vExB_etheta, vD_etheta
    real(RKIND), dimension(:,:,:), pointer :: dPhidr
    real(RKIND), dimension(:,:,:), pointer :: dPhidtheta
    real(RKIND), dimension(:,:,:), pointer :: dPhidphi
    real(RKIND), dimension(:)    , pointer :: vExB_gradr_1D
    real(RKIND), dimension(:)    , pointer :: vExB_gradtheta_1D
    real(RKIND), dimension(:)    , pointer :: vD_gradr_1D
    real(RKIND), dimension(:)    , pointer :: vD_gradtheta_1D
    !-> variables for the theta component of the magnetisation
    real(RKIND)                        :: bphi  
    real(RKIND), dimension(:), pointer :: intMphi_dthetadphi
    !-> variables for potential energy
    real(RKIND) :: Enpot
    real(RKIND), dimension(1:Nbthread) :: Enpot_thread 
    !-> variables for J0.f
    complex(CKIND), dimension(:,:), pointer :: Acomp
    real(RKIND)   , dimension(:,:), pointer :: J0f_rtheta
    real(RKIND)   , dimension(:,:), pointer :: J0_fvExBgradr
    ! -> variables for time perfomance analysis
    integer(TIMEPREC) :: tdeb, tfin
    real(RKIND)       :: tdPhidx, tgyro, tcomm, tcomp1, tcomp2
      
    !*** computation of the moments associated to the particles ***
    !-> BECAREFUL: if you want to add the computation of
    !    a particle moment, do not forget to change
    !    nbmoments_part in 'globals.f90' file
    nbmoments          = nbmoments_part
    intMphi_dthetadphi => Rarray1_Nr
      
    !*** initialisation to 0 of the radial arrays ***
    do imom = 1,nbmoments
      do ir =  0,geom%Nr
        Fmoments_loc(ir,imom)  = 0._RKIND
        part_Fmoments(ir,imom) = 0._RKIND
      end do
    end do
      
    mu_value          = fnp1%mu_value
    coeff_intdmu_tmp  = geom%coeff_intdmu(mu_id)
      
    !*** computation of the derivatives of Phi
    call clck_time(tdeb)
    call ppbarrier_timer()
    dPhidtheta => Rarray1_NrNthetaNphi
    dPhidphi   => Rarray2_NrNthetaNphi
    call compute_dPhidtheta_loc(geom,sPhi, &
      0,geom%Nr,0,geom%Ntheta,kstart,kend,dPhidtheta)
    call compute_dPhidphi_loc(geom,sPhi, &
       0,geom%Nr,0,geom%Ntheta,0,geom%Nphi,dPhidphi)
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tdPhidx)
      
    !*** computation of J0.f ***
    call clck_time(tdeb)
    
    !-> transposition of f to have all the informations in (r,theta)
    call pptranspose_forward(fnp1%values)
      
    RStheta => Romp5_0Nr
#ifdef _OPENMP
!$OMP PARALLEL private(tid,vExB_gradr_1D, &
!$OMP Acomp,J0f_rtheta,J0_fvExBgradr, &
!$OMP ir,itheta,iphi,ivpar, &
!$OMP coeff_intdvpar_tmp,coeff_intdphi_tmp,fval, &
!$OMP coeff_intdtheta_tmp,jacob_space_tmp, &
!$OMP coeff_int_volume,jacob_vel_tmp) &
!$OMP default(shared)
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
    vExB_gradr_1D => Romp1_0Nr(tid)%val
    Acomp         => Comp1_1Nrp1_1Ntheta(tid)%val
    J0f_rtheta    => Romp1_1Nrp1_1Nthetap1(tid)%val
    J0_fvExBgradr => Romp2_1Nrp1_1Nthetap1(tid)%val
    do ir = 0,fnp1%n1
      RStheta(tid)%val(ir) = 0._RKIND
    end do
!$OMP DO SCHEDULE(STATIC)
    do ivpar = lstart,lend
      coeff_intdvpar_tmp  = geom%coeff_intdvpar(ivpar)
      do iphi = kstart,kend
        coeff_intdphi_tmp   = geom%coeff_intdphi(iphi)
        do itheta = 0,fnp1%n2
          !-> computation of vExB_gradr(0:Nr)
          call compute_vExBgradr(geom,init_magnet,init_curr, &
            dPhidtheta,dPhidphi,0,geom%Nr,itheta,itheta, &
            iphi,iphi,ivpar,ivpar,vExB_gradr_1D)
          do ir = 0,fnp1%n1
            fval = f4D_transp(ir,itheta,iphi,ivpar)
            J0f_rtheta(ir+1,itheta+1)    = fval
            J0_fvExBgradr(ir+1,itheta+1) = fval*vExB_gradr_1D(ir)
          end do
        end do
        !--> computation of J0.f
        call omp_compute_gyrofunction_2D(J0,geom, &
          mu_id,Acomp,J0f_rtheta)
        !--> computation of J0.(f.vExB_gradr)
        call omp_compute_gyrofunction_2D(J0,geom, &
          mu_id,Acomp,J0_fvExBgradr)
        
        !-> \int J0.(f*vExB.gradr) d3v (for RStheta)
        do itheta = 0,geom%Ntheta
          coeff_intdtheta_tmp = geom%coeff_intdtheta(itheta)
          do ir = 0,geom%Nr
            !-> jacobian_velocity = 2*pi*Bstar
            call compute_jacobian_velocity(geom,init_magnet, &
              init_curr,ir,itheta,ivpar,jacob_vel_tmp)            
            !-> jacobian in space
            jacob_space_tmp = jacobian_space(ir,itheta)
            !-> compute the elements of volume in space and velocity
            coeff_int_volume = jacob_space_tmp * &
              coeff_intdtheta_tmp*coeff_intdphi_tmp * &
              jacob_vel_tmp*coeff_intdvpar_tmp*coeff_intdmu_tmp 
            RStheta(tid)%val(ir) =  RStheta(tid)%val(ir) + &
              J0_fvExBgradr(ir+1,itheta+1) * coeff_int_volume
          end do
        end do
        do itheta = 0,fnp1%n2
          do ir = 0,fnp1%n1
            f4D_transp(ir,itheta,iphi,ivpar) = &
              J0f_rtheta(ir+1,itheta+1)
          end do
        end do        
      end do
    end do
!$OMP END DO
!$OMP END PARALLEL
    !->  f4D_transp -> J0ftmp%values
    call pptranspose_backward(sJ0ftmp%values)
    
    do ir = 0, geom%Nr
      RStheta_tmp = 0._RKIND
      do ithread = 1,Nbthread
        RStheta_tmp = RStheta_tmp + &
          RStheta(ithread)%val(ir)
      end do
      Fmoments_loc(ir,8) = RStheta_tmp
    end do
      
    call ppbarrier_timer()
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tgyro)
      
    !*** computation of the derivatives of Phi
    call clck_time(tdeb)
    call ppbarrier_timer()
    dPhidr     => Rarray1_PNrPNthetaNphi
    dPhidtheta => Rarray2_PNrPNthetaNphi
    dPhidphi   => Rarray3_PNrPNthetaNphi
    call compute_dPhidr_loc(geom,sPhi, &
      istart,iend,jstart,jend,0,geom%Nphi,dPhidr)
    call compute_dPhidtheta_loc(geom,sPhi, &
      istart,iend,jstart,jend,0,geom%Nphi,dPhidtheta)
    call compute_dPhidphi_loc(geom,sPhi, &
      istart,iend,jstart,jend,0,geom%Nphi,dPhidphi)
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tdPhidx)
      
    call clck_time(tdeb)
      
#ifdef _OPENMP
!$OMP PARALLEL private(tid,ir,itheta,iphi,ivpar,imom, &
!$OMP Bij,sqrt_g22,vpar,energy,J0f,J0fmu_eq,J0fv,J0fE, &
!$OMP coeff_intdr_tmp,coeff_intdtheta_tmp, &
!$OMP coeff_intdphi_tmp,coeff_intdvpar_tmp, &
!$OMP coeff_int_velocities,coeff_int_volume, &
!$OMP jacob_space_tmp,jacob_vel_tmp, &
!$OMP vExB_gradr_1D,vExB_gradtheta_1D, &
!$OMP vD_gradr_1D,vD_gradtheta_1D, bphi, &
!$OMP vExB_er,vExB_etheta,vD_er,vD_etheta) &
!$OMP default(shared)
!$OMP BARRIER
    tid = 1+omp_get_thread_num()
#else
    tid = 1
#endif
      
    vExB_gradr_1D     => Romp1_0Nr(tid)%val
    vExB_gradtheta_1D => Romp2_0Nr(tid)%val
    vD_gradr_1D       => Romp3_0Nr(tid)%val
    vD_gradtheta_1D   => Romp4_0Nr(tid)%val
      
    !*** initialisation of the values to 0 ***
    do imom = 1,nbmoments
      do ir = 0,geom%Nr
        Romp1_1nbmom_0Nr(tid)%val(imom,ir) = 0._RKIND
      end do
    end do
    Enpot_thread(tid) = 0._RKIND
      
!$OMP DO SCHEDULE(STATIC)
    do ivpar = 0,geom%Nvpar
      vpar                = geom%vparg(ivpar)
      coeff_intdvpar_tmp  = geom%coeff_intdvpar(ivpar)
      do iphi = 0,geom%Nphi-1
        coeff_intdphi_tmp = geom%coeff_intdphi(iphi)
        do itheta = jstart,jend
          coeff_intdtheta_tmp = geom%coeff_intdtheta(itheta)
          !*** computation of the radial velocities    ***
          !-> computation of vExB_gradr(0:Nr)
          call compute_vExBgradr(geom,init_magnet,init_curr, &
            dPhidtheta,dPhidphi,istart,iend,itheta,itheta, &
            iphi,iphi,ivpar,ivpar,vExB_gradr_1D)
          !-> computation of vExB_gradtheta(0:Nr)
          call compute_vExBgradtheta(geom, &
            init_magnet,init_curr, &
            dPhidr,dPhidphi,istart,iend,itheta,itheta, &
            iphi,iphi,ivpar,ivpar,vExB_gradtheta_1D)
          !-> computation of vD_gradr(0:Nr)
          call compute_vDgradr(geom,init_magnet,init_curr, &
            istart,iend,itheta,itheta,ivpar,ivpar,vD_gradr_1D)
          !-> computation of vD_gradtheta(0:Nr)
          call compute_vDgradtheta(geom,init_magnet,init_curr, &
            istart,iend,itheta,itheta,ivpar,ivpar,vD_gradtheta_1D)
          do ir = istart,iend
            coeff_intdr_tmp = geom%coeff_intdr(ir)
            Bij             = init_magnet%B_norm(ir,itheta)
            energy          = 0.5_RKIND*vpar*vpar+mu_value*Bij
            J0f             = sJ0ftmp%values(ir,itheta,iphi,ivpar)
            J0fmu_eq        = sJ0fmu_eq(ir,itheta,ivpar)
            J0fv            = J0f*vpar
            J0fE            = J0f*energy
            !-> jacobian_velocity = 2*pi*Bstar
            call compute_jacobian_velocity(geom,init_magnet, &
              init_curr,ir,itheta,ivpar,jacob_vel_tmp)            
            !-> jacobian in space
            jacob_space_tmp = jacobian_space(ir,itheta)
            !-> compute the elements of velocity
            coeff_int_velocities = jacob_vel_tmp * &
              coeff_intdvpar_tmp*coeff_intdmu_tmp 
            !-> compute the elements in the total phase space
            coeff_int_volume = jacob_space_tmp * &
              coeff_intdtheta_tmp*coeff_intdphi_tmp * &
              coeff_int_velocities
            !*** computation of the projection of the   ***
            !***  radial velocities on the unit radial  ***
            !***  vector                                ***
            sqrt_g22 = sqrt(coord_sys%g22(ir,itheta)) 
            !-> vE.er = vExB_gradr
            vExB_er = vExB_gradr_1D(ir-istart)
            !-> vE.etheta = sqrt(g_{22})*vExB_gradtheta
            vExB_etheta = sqrt_g22*vExB_gradtheta_1D(ir-istart)
            !-> vD.er = vD_gradr 
            vD_er = vD_gradr_1D(ir-istart)
            !-> vD.etheta = sqrt(g_{22})*vD_gradtheta
            vD_etheta = sqrt_g22*vD_gradtheta_1D(ir-istart)
            
            !*** computation of b_phi = B_phi/B ***
            bphi = init_magnet%Bphi(ir,itheta)/Bij
      
            !*** computation of the integrals in velocity and ***
            !***  in theta and phi                            ***
            !-> \int J0.f d3v  (for Ni)
            Romp1_1nbmom_0Nr(tid)%val(1,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(1,ir) + &
              J0f*coeff_int_volume
            !-> \int J0.f*vpar d3v (for nVpar)
            Romp1_1nbmom_0Nr(tid)%val(2,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(2,ir) + &
              J0fv*coeff_int_volume
            !-> \int J0.f*E d3v (for P)
            Romp1_1nbmom_0Nr(tid)%val(3,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(3,ir) + &
              J0fE*coeff_int_volume
            !-> \int J0.f*vE.er d3v (for Gamma_turb)
            Romp1_1nbmom_0Nr(tid)%val(4,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(4,ir) + &
              J0f*vExB_er*coeff_int_volume
            !-> \int J0.f*vD.er d3v (for Gamma_neo)
            Romp1_1nbmom_0Nr(tid)%val(5,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(5,ir) + &
              J0f*vD_er*coeff_int_volume
            !-> \int J0.f*E*vE.er d3v (for Q_turb)
            Romp1_1nbmom_0Nr(tid)%val(6,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(6,ir) + &
              J0fE*vExB_er*coeff_int_volume
            !-> \int J0.f*E*vD.er d3v (for Q_neo)
            Romp1_1nbmom_0Nr(tid)%val(7,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(7,ir) + &
              J0fE*vD_er*coeff_int_volume
            !-> Romp1_1nbmom_0Nr(tid)%val(8,ir) is already computed
            !-> \int J0.f vtilde_r vpar d3v (for RSpar)
            Romp1_1nbmom_0Nr(tid)%val(9,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(9,ir) + &
              vExB_er*J0fv*coeff_int_volume
            !-> \int M_phi dtheta dphi (for nvpolor_mag)
            !->  where M_phi = - \int (\int J0.f b_phi 
            !->                 Jv dvpar)*mu dmu            
            Romp1_1nbmom_0Nr(tid)%val(10,ir) = &
              Romp1_1nbmom_0Nr(tid)%val(10,ir) - &
              real(Zi)*J0f*mu_value*bphi * &
              coeff_int_velocities * &
              coeff_intdtheta_tmp*coeff_intdphi_tmp
            !-> compute \int JO.f Phi d3v
            Enpot_thread(tid) = Enpot_thread(tid) + &
              0.5_RKIND*(J0f-J0fmu_eq)*sPhi(ir,itheta,iphi) * &
               coeff_int_volume*coeff_intdr_tmp            
          end do ! end of do ir
        end do   ! end of do itheta
      end do     ! end of do iphi
    end do       ! end of do ivpar
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
      
    !*** fill an unique array containing all the       ***
    !***  radial profiles of the fluid moments         ***
    !*** Rq: each processor fills the part istart:iend ***
    !***     and all the others values are put equal   ***
    !***     to 0                                      ***
    Enpot = 0._RKIND
    do ithread = 1,Nbthread
      do ir = istart,iend
        do imom = 1,nbmoments
          Fmoments_loc(ir,imom) = &
            Fmoments_loc(ir,imom) + & 
            Romp1_1nbmom_0Nr(ithread)%val(imom,ir)
        end do
      end do
      Enpot = Enpot + Enpot_thread(ithread)
    end do
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcomp1)
      
    !*** MPI_REDUCE ****
    call clck_time(tdeb)
    if (Nbproc_tot.ne.1) then
      nbelements = nbmoments*(geom%Nr+1)
      call MPI_REDUCE(Fmoments_loc,part_Fmoments, &
        nbelements,MPI_REAL8,MPI_SUM,diag_targ(3), &
        MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(Enpot,Enpot_diag,1, &
        MPI_REAL8,MPI_SUM,diag_targ(1), &
        MPI_COMM_WORLD,ierr)
    else
      do imom = 1,nbmoments
        do ir = 0,geom%Nr
          part_Fmoments(ir,imom) = Fmoments_loc(ir,imom)
        end do
      end do
      Enpot_diag = Enpot
    end if
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcomm)
      
    !*** computation of the flux surface average ***
    call clck_time(tdeb)
    coeff_pressr = 1._RKIND/(0.5_RKIND+float(ipow))
    
    do ir = 0,geom%Nr
      intdthetadphi_Js_tmp   = intdthetadphi_Js(ir)
      !-> compute Ni(r) = <\int J0.f d3v>_FS
      nir_diag(ir)           = & 
        part_Fmoments(ir,1)/intdthetadphi_Js_tmp
      !-> compute nVpar(r) = <\int J0.f*vpar d3v>_FS
      nuparr_diag(ir)        = &
        part_Fmoments(ir,2)/intdthetadphi_Js_tmp
      !-> compute P(r) = <\int J0.f*E d3v>_FS
      pressr_diag(ir)        = coeff_pressr * &
        part_Fmoments(ir,3)/intdthetadphi_Js_tmp
      !-> compute Gamma_turb(r) = <\int J0.f*vE.er d3v>_FS
      Gammar_turb_diag(ir)   = &
        part_Fmoments(ir,4)/intdthetadphi_Js_tmp
      !-> compute Gamma_neo(r) = <\int J0.f*vD.er d3v>_FS
      Gammar_neo_diag(ir)    = &
        part_Fmoments(ir,5)/intdthetadphi_Js_tmp
      !-> compute Q_turb(r) = <\int J0.f*E*vE.er d3v>_FS
      Qr_turb_diag(ir)       = &
        part_Fmoments(ir,6)/intdthetadphi_Js_tmp
      !-> compute Q_neo(r) = <\int J0.f*E*vD.er d3v>_FS 
      Qr_neo_diag(ir)        = &
        part_Fmoments(ir,7)/intdthetadphi_Js_tmp
      !-> compute RStheta(r) = <\int J0.(f*vExB.gradr) d3v>_FS
      RSthetar_diag(ir)      = &
        part_Fmoments(ir,8)/intdthetadphi_Js_tmp
      !-> compute RSpar(r) = <\int J0.f vtilde_r vpar d3v>_FS
      RSparr_diag(ir)        = &
        part_Fmoments(ir,9)/intdthetadphi_Js_tmp
      !-> compute \int M_phi dtheta dphi 
      intMphi_dthetadphi(ir) = part_Fmoments(ir,10)
    end do
 
    !*** Computation of the flux surface average of the  ***
    !***  theta component of the magnetisation, i.e:     ***
    !***  nvpolor_mag = -<rot(M).etheta>_FS              ***
    !***              = -r*<rot(M).gradtheta>_FS         ***
    !***              = r d/dr(\int M_phi dtheta dphi) / ***
    !***                     \int Js dtheta dphi         ***
    !***  because rot(M).gradtheta = -1/Js d/dr(M_phi)   ***
    !***   with Js is the jacobian in space and          ***
    !-> A(r) = d/dr (\int M_phi dtheta dphi)
    call deriv1(intMphi_dthetadphi,nvpolor_mag_diag, &
      geom%Nr,geom%dr,0)
    !-> nvpolor_mag(r) = r*A(r) / \int Js dtheta dphi
    do ir = 0,geom%Nr
      intdthetadphi_Js_tmp = intdthetadphi_Js(ir)
      !---> ATTENTION (no dependence of sqrt_g22 in theta ?)
      sqrt_g22             = sqrt(coord_sys%g22(ir,0))
      nvpolor_mag_diag(ir) = -sqrt_g22 * &
        nvpolor_mag_diag(ir)/intdthetadphi_Js_tmp
    end do
    call ppbarrier_timer()
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tcomp2)
      
#ifdef TIMER    
    write(6,'(I4,A,F13.7,A,F13.7,A,F13.7,A,F13.7,A,F13.7)') &
      pglobal_id," Temps part moments: dPhidx ",tdPhidx, &
      " J0f ",tgyro, " comp1 ",tcomp1, &
      " comm ", tcomm, " comp2 ", tcomp2
#endif
  end subroutine compute_part_moments
end module physics_module
