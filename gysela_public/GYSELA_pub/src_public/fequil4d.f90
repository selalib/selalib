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
! file : fequil4d.f90
! date : 19/05/2006
! - computation of the initial contribution of the 
!  equilibrium part of the distribution function :
!   feq(r,theta,vpar) for the value of mu equal to mu_id
!-------------------------------------------------------
module fequil4d_module
  use prec_const
  use MPIutils_module
  use mem_alloc_module
  use spline1d_class
  use geometry_class
  use init_profile_class
  use init_magnetic_class
      
  implicit none      
      
  integer        , private :: istart_fmu, iend_fmu
  integer        , private :: jstart_fmu, jend_fmu
  integer        , private :: f4d_BCr_left, f4d_BCr_right
  integer        , private :: f4d_BCvpar_left, f4d_BCvpar_right
  type(nspline1d), private :: nspline1d_r
  type(pspline1d), private :: pspline1d_theta
  type(nspline1d), private :: nspline1d_vpar
      
  include "Bstar_inline.h"
  include "velocities_inline.h"
  !******************************
  contains
  !******************************
#include "Bstar_inline.f90"
#include "velocities_inline.f90"
      
  !*********************************************************
  ! EQUILIBRIUM PART OF THE DISTRIBUTION FUNCTION
  !*********************************************************
  !----------------------------------------------------------
  ! allocate the arrays associated to the equilibrium part
  !----------------------------------------------------------
  subroutine allocate_equilibrium_part(geom,Sfmu_eq,SJ0fmu_eq, &
    M2_equil,nGieq_rtheta, neq_r,dneq_dr,Tieq_r)
    use globals, only : BC_Hermite, istart, iend, &
      jstart, jend, hffilterfref, hffilterfreq
    type(geometry)               , intent(in) :: geom
    real(RKIND), dimension(:,:,:), pointer    :: Sfmu_eq
    real(RKIND), dimension(:,:,:), pointer    :: SJ0fmu_eq
    real(RKIND), dimension(:,:)  , pointer    :: M2_equil
    real(RKIND), dimension(:,:)  , pointer    :: nGieq_rtheta 
    real(RKIND), dimension(:)    , pointer    :: neq_r
    real(RKIND), dimension(:)    , pointer    :: dneq_dr
    real(RKIND), dimension(:)    , pointer    :: Tieq_r
      
!baoter
    if ((hffilterfref.eq.1).and.(hffilterfreq.gt.0)) then
      istart_fmu = 0
      iend_fmu   = geom%Nr
      jstart_fmu = 0
      jend_fmu   = geom%Ntheta
    else
      istart_fmu = istart
      iend_fmu   = iend
      jstart_fmu = jstart
      jend_fmu   = jend
    end if
!eaoter
    f4d_BCr_left     = BC_Hermite
    f4d_BCr_right    = BC_Hermite
    f4d_BCvpar_left  = BC_Hermite
    f4d_BCvpar_right = BC_Hermite
    call new_spline1d_natural(nspline1d_r,geom%Nr,geom%dr)
    call new_spline1d_period(pspline1d_theta, &
      geom%Ntheta,geom%dtheta)
    call new_spline1d_natural(nspline1d_vpar,geom%Nvpar,geom%dvpar)
    call glob_allocate(Sfmu_eq,istart_fmu,iend_fmu, &
      jstart_fmu,jend_fmu,0,geom%Nvpar,'Sfmu_eq')
    call glob_allocate(SJ0fmu_eq,istart,iend, &
      jstart,jend,0,geom%Nvpar,'SJ0fmu_eq')
    call glob_allocate(M2_equil,0,geom%Nr,0,geom%Ntheta,'M2_equil')
    call glob_allocate(nGieq_rtheta,0,geom%Nr, &
      0,geom%Ntheta,'nGieq_rtheta')
    call glob_allocate(neq_r,0,geom%Nr,'neq_r')
    call glob_allocate(dneq_dr,0,geom%Nr,'dneq_dr')
    call glob_allocate(Tieq_r,0,geom%Nr,'Tieq_r')
  end subroutine allocate_equilibrium_part
      
  !---------------------------------------------------------
  ! delete the arrays associated to the equilibrium part
  !---------------------------------------------------------
  subroutine deallocate_equilibrium_part(Sfmu_eq,SJ0fmu_eq, &
    M2_equil,nGieq_rtheta,neq_r,dneq_dr,Tieq_r)
    real(RKIND), dimension(:,:,:), pointer :: Sfmu_eq 
    real(RKIND), dimension(:,:,:), pointer :: SJ0fmu_eq
    real(RKIND), dimension(:,:)  , pointer :: M2_equil
    real(RKIND), dimension(:,:)  , pointer :: nGieq_rtheta 
    real(RKIND), dimension(:)    , pointer :: neq_r 
    real(RKIND), dimension(:)    , pointer :: dneq_dr
    real(RKIND), dimension(:)    , pointer :: Tieq_r
      
    call del_spline1d_natural(nspline1d_r)
    call del_spline1d_period(pspline1d_theta)
    call del_spline1d_natural(nspline1d_vpar)
      
    call glob_deallocate(Sfmu_eq)
    call glob_deallocate(SJ0fmu_eq)
    call glob_deallocate(M2_equil)
    call glob_deallocate(nGieq_rtheta)
    call glob_deallocate(neq_r)
    call glob_deallocate(dneq_dr)
    call glob_deallocate(Tieq_r)
  end subroutine deallocate_equilibrium_part
      
  !-----------------------------------------------------------
  ! Computation of rbar = rp-(q_rp/rp)*[psi(r)-psi(rp)]
  !    -(q_rp/rp)*[vpar*R-vpar0*R0*H(energy-mu*Bmax)] 
  !  where . rp     = r(peak)
  !        . psi(r) = -\int_0^r r/q(r) dr
  !        . vpar0  = sign(vpar) * 
  !                   sqrt(vpar^2+2*mu*[B(r,theta)-Bmin])
  !        . Bmin   = min(B)
  ! with H(energy-mu*Bmax) is close to the heaviside 
  !  ( => tanh for having a smoothness transition between 
  !     trapped and passing particles)
  !  Rk : H is equal to 0 for the trapped particles
  !-----------------------------------------------------------
  subroutine compute_rbar(geom,ir,itheta,ivpar,imu, &
    init_prof,init_magnet,rbar,psibar)
    use globals, only : R0, a, q0, deltaq, alphaq, &
      feq_choice, vpar0_option
    type(geometry)     , intent(in)  :: geom
    integer            , intent(in)  :: ir, itheta, ivpar, imu
    type(init_profile) , intent(in)  :: init_prof
    type(init_magnetic), intent(in)  :: init_magnet
    real(RKIND)        , intent(out) :: rbar, psibar
      
    real(RKIND) :: rp, q_rp, psi_rp
    real(RKIND) :: vparl, Bmin, Bmax, Bij, vpar0, mum
    real(RKIND) :: energy, part_type, H_passing, v_tmp, sign_vpar
      
    rp        = init_prof%rpeak
    q_rp      = init_prof%q_rp 
    psi_rp    = init_prof%psi_rp
    ! -> vpar0 = sign(vpar)*sqrt(vpar^2+2*mu*[B(r,theta)-Bmin])
    Bmin      = init_magnet%Bmin
    Bmax      = init_magnet%Bmax
    Bij       = init_magnet%B_norm(ir,itheta)
    vparl     = geom%vparg(ivpar)
    mum       = geom%mug(imu)
    energy    = 0.5_RKIND*vparl**2 + mum*Bij
    sign_vpar = sign(1._RKIND,vparl)
    
    !*** computation of vpar0 ***
    select case (vpar0_option)
    case(0)
      ! ->  vpar0 = 0 
      vpar0     = 0._RKIND
      H_passing = 0._RKIND
    case(1)  
      ! -> vpar0 = sign(vpar)*sqrt(2*[E-mu*Bmin])*H(E-mu*Bmax)
      part_type = energy-mum*Bmax
      if (part_type.le.0._RKIND) then
        vpar0     = 0._RKIND
        H_passing = 0._RKIND
      else
        v_tmp     = sqrt(2._RKIND*(energy-mum*Bmin))
        vpar0     = sign_vpar*v_tmp
        H_passing = 1._RKIND
      end if
    case(2)  
      ! -> vpar0 = sign(vpar)*sqrt(2*[E-mu*Bmax])*H(E-mu*Bmax)
      part_type = energy-mum*Bmax
      if (part_type.le.0._RKIND) then
        vpar0     = 0._RKIND
        H_passing = 0._RKIND
      else
        v_tmp     = sqrt(2._RKIND*(energy-mum*Bmax))
        vpar0     = sign_vpar*v_tmp
        H_passing = 1._RKIND
      end if
    case(3)
      ! -> vpar0 = sign(vpar)*sqrt(2*[E-mu*B0])*H(E-mu*Bmax)
      part_type = energy-mum*Bmax
      if (part_type.le.0._RKIND) then
        vpar0     = 0._RKIND
        H_passing = 0._RKIND
      else
        v_tmp     = sqrt(2._RKIND*(energy-mum*1._RKIND))
        vpar0     = sign_vpar*v_tmp
        H_passing = (tanh((part_type-3._RKIND*geom%dvpar) / &
          (geom%dvpar)) + 1._RKIND) * 0.5_RKIND 
      end if
    case(4)
      ! -> vpar0 = sign(vpar)*sqrt(2*[E-mu*B0])*H(E-mu*B0)
      part_type = energy-mum*1._RKIND
      if (part_type.le.0._RKIND) then
        vpar0     = 0._RKIND
        H_passing = 0._RKIND
      else
        v_tmp     = sqrt(2._RKIND*(energy-mum*1._RKIND))
        vpar0     = sign_vpar*v_tmp
        H_passing = 1._RKIND
      end if
    case default
      if (pglobal_id.eq.0) then
        print*, 'option vpar0_option = ', &
          vpar0_option, ' not treated'
      end if
    end select
      
    psibar = init_prof%psi(ir)+R0*(vparl/Bij-H_passing*vpar0)
    if (alphaq.ne.2._RKIND) then
      ! -> rbar = rp-(q_rp/rp)*[psi(r)-psi(rp)- 
      !    R0*vpar/(B(ir,itheta)-vpar0*R0)] 
      rbar = rp-(q_rp/rp) * &
        (init_prof%psi(ir)-psi_rp+R0*vparl/Bij-H_passing*vpar0*R0)
    else
      ! in the case of a parobolic profile for q:
      !  q=q0+deltaq*(r/a)^2 the expression of r can be 
      !  deduced directly from psi(r) so the correction is 
      !  made on psi: psibar = psi + corr and then 
      !  rbar is obtained with psibar
      ! -> psibar = psi-R0*vpar/(q*B(ir,itheta))-vpar0*R0/q)] 
      ! -> rbar = sqrt(q0a^2/deltaq * 
      !           (exp(-2._RKIND*deltaq*psibar/a^2)-1))
      if (psibar.lt.0) then
        rbar = sqrt(q0*a**2/deltaq * &
          (exp(-2._RKIND*deltaq*psibar/a**2)-1))      
      else
        !-> in this case rbar < rmin so is put equal to rmin 
        rbar = geom%rg(0)
      end if
    end if
  end subroutine compute_rbar
      
  !------------------------------------------------------------
  ! Computation of the equilibrium part feq of the 
  !   distribution function at the point (ir,itheta,ivpar,imu) :
  !    feq(rbar,theta,vpar,mu) = n0(rbar)/(2*PI*Ti(rbar))**(3/2)
  !                *exp(-1/Ti(rbar)*(1/2*vpar**2+mu*B(r,theta)))
  !  with rbar = (q_rp/rp)*[psi(r)-R0*vpar/(B(ir,itheta)] 
  !  where rp = r(peak)
  !------------------------------------------------------------
  subroutine compute_feq_value(geom,ir,itheta,ivpar,imu, &
    init_prof,init_magnet,feq)
    use globals, only : Rarray1_Nr, &
      rbar_in_feq, feq_choice, &
      min_err_rbar, max_err_rbar, ipow, use_psibar
    use utils_module
    type(geometry)     , intent(in)  :: geom
    integer            , intent(in)  :: ir, itheta, ivpar, imu
    type(init_profile) , intent(in)  :: init_prof
    type(init_magnetic), intent(in)  :: init_magnet
    real(RKIND)        , intent(out) :: feq
      
    integer     :: ipos, ibase
    integer     :: ipsimin, ipsimax, ir_loc
    real(RKIND) :: denom, Ti_r, Ti_rp, n0_rp
    real(RKIND) :: rbar, Ti_rbar, n0_rbar
    real(RKIND) :: psibar, psibar_minus_psi, &
      test_inside, deltar, psimin, psimax
    real(RKIND) :: energy
    real(RKIND), dimension(-1:2) :: sbase
      
    if (.not.rbar_in_feq) then
      Ti_r   = init_prof%Ti(ir)
      denom  = sqrt(2._RKIND*PI*Ti_r)*(2._RKIND*PI*Ti_r)**ipow
      energy = 0.5_RKIND*geom%vparg(ivpar)**2 + &
        geom%mug(imu)*init_magnet%B_norm(ir,itheta)
      feq    = init_prof%n0(ir)/denom * exp(-energy/Ti_r)
    else
      !*** computation of rbar ***
      call compute_rbar(geom,ir,itheta,ivpar,imu, &
        init_prof,init_magnet,rbar,psibar)
      min_err_rbar = min(min_err_rbar,abs(rbar-geom%rg(ir)))
      max_err_rbar = max(max_err_rbar,abs(rbar-geom%rg(ir)))
      if (use_psibar) then
        psimin  = init_prof%psimin
        psimax  = init_prof%psimax
        ipsimin = init_prof%ipsimin
        ipsimax = init_prof%ipsimax
        test_inside = psibar-psimin
        if (test_inside.le.0._RKIND) then
          deltar  = -test_inside/(geom%rg(ipsimin) * &
            init_prof%iota(ipsimin))
          Ti_rbar = init_prof%Ti(ipsimin)+deltar * &
            init_prof%dTidr(ipsimin)
          n0_rbar = init_prof%n0(ipsimin)+deltar * &
            init_prof%dn0dr(ipsimin)
          Ti_rbar = max(Ti_rbar,0._RKIND)
          n0_rbar = max(n0_rbar,0._RKIND)
        else
          test_inside = psimax-psibar
          if (test_inside.le.0._RKIND) then
            deltar  = test_inside/(geom%rg(ipsimax) * &
              init_prof%iota(ipsimax))
            if (abs(deltar).le.geom%rg(0)) then
              Ti_rbar = init_prof%Ti(ipsimax) + &
                deltar*init_prof%dTidr(ipsimax)
              n0_rbar = init_prof%n0(ipsimax) + &
                deltar*init_prof%dn0dr(ipsimax)
            else if (abs(deltar).le.2._RKIND*geom%rg(0)) then
              Ti_rbar = init_prof%Ti(ipsimax)- &
                (2._RKIND*geom%rg(0)-abs(deltar)) * &
                init_prof%dTidr(ipsimax)
              n0_rbar = init_prof%n0(ipsimax)- &
                (2._RKIND*geom%rg(0)-abs(deltar)) * &
                init_prof%dn0dr(ipsimax)
            else
              ipos = floor((abs(deltar)-2._RKIND*geom%rg(0)) / &
                geom%dr)
              if (ipos.ge.geom%Nr) then
                Ti_rbar = init_prof%Ti(ipsimin)
                n0_rbar = init_prof%n0(ipsimin)
              else
                Ti_rbar = init_prof%Ti(ipos)
                n0_rbar = init_prof%n0(ipos)
              end if
            end if
          else
            !*** Looking for the nearest psi value
            !*** interpolation (2nd order Taylor expansion)
            test_inside = abs(psibar-psimin)
            do ipos = 0,geom%Nr
              psibar_minus_psi = psibar-init_prof%psi(ipos)
              if (abs(psibar_minus_psi).le.test_inside) then
                deltar  = -psibar_minus_psi / &
                  (geom%rg(ipos)*init_prof%iota(ipos))
                Ti_rbar = init_prof%Ti(ipos) + &
                  deltar*init_prof%dTidr(ipos) + &
                  deltar*deltar*init_prof%d2Tidr2(ipos)/2._RKIND
                n0_rbar = init_prof%n0(ipos) + &
                  deltar*init_prof%dn0dr(ipos) + &
                  deltar*deltar*init_prof%d2n0dr2(ipos)/2._RKIND
                test_inside = abs(psibar_minus_psi)
              end if
            end do
          end if
        end if
      else
        !*** computation of Ti(rbar) and n0(rbar) ***
        !***  by cubic spline interpolation       ***
        ! ->  ipos = position in rg of the first knot 
        !      at the left of rbar
        ipos = floor( (rbar-geom%rg(0))/geom%dr)
        if (ipos.le.0) then
          Ti_rbar = init_prof%Ti(0)
          n0_rbar = init_prof%n0(0)
        else
          if (ipos.gt.geom%Nr) then
            Ti_rbar = init_prof%Ti(geom%Nr)
            n0_rbar = init_prof%n0(geom%Nr)
          else
            ! -> cubic spline interpolation
            if (ipos.eq.geom%Nr) ipos = geom%Nr-1
            call spline_basis(geom%rg(ipos),rbar, &
              geom%rg(ipos+1),geom%dr,sbase)
            Ti_rbar = 0._RKIND
            n0_rbar = 0._RKIND
            do ibase = -1,2
              Ti_rbar = Ti_rbar + &
                init_prof%Ti_coef1d(ipos+ibase)*sbase(ibase)
              n0_rbar = n0_rbar + &
                init_prof%n0_coef1d(ipos+ibase)*sbase(ibase)
            enddo
          end if
        end if
      end if
      !*** computation of feq depending on rbar ***
      denom  = sqrt(2._RKIND*PI*Ti_rbar) * &
        (2._RKIND*PI*Ti_rbar)**ipow
      energy = 0.5_RKIND*geom%vparg(ivpar)**2 + &
        geom%mug(imu)*init_magnet%B_norm(ir,itheta) 
      Ti_rp  = init_prof%Ti_rp
      n0_rp  = init_prof%n0_rp
      select case (feq_choice)
      case (1,5)
        feq = n0_rbar/denom * exp(-energy/Ti_rbar)
      case (2,4,6)
        feq = exp(-(rbar/init_prof%rpeak)**2)
      case (3)
        feq = exp(-energy)
      case (7)
        denom  = sqrt(2._RKIND*PI*Ti_rp)*(2._RKIND*PI*Ti_rp)**ipow
        feq = (n0_rbar+(energy/Ti_rp-0.5_RKIND-float(ipow)) &
          *n0_rbar*(Ti_rbar-Ti_rp)/Ti_rp)/denom*exp(-energy/Ti_rp)
      case (8)
        denom  = sqrt(2._RKIND*PI*Ti_rp)*(2._RKIND*PI*Ti_rp)**ipow
        feq    = n0_rbar*exp(-energy/Ti_rp)/denom
      end select
    end if
    if ( abs(feq).le.1.e-30_RKIND ) feq = 1.e-30_RKIND
  end subroutine compute_feq_value
      
  !-------------------------------------------------------------
  ! computation of the equilibrium terms (global variables):
  !    . feq(r,theta,vpar,mu) the equilibrium part of f
  !    . J0.feq the gyroaveraged feq
  !    . nGi_eq(r,theta) = \int J0.feq Jv dvpar dmu
  !        where Jv the Jacobian in velocity space 
  !        is equal to:   
  !          Jv = 2*pi*Bstar(r,theta,vpar) in 5D case
  !             = 1                        in 4D case
  !                        
  !    . n_eq(r)   = (1/2pi)\int  nGi_eq(r,theta) dtheta
  !    . P_eq(r)   = press_coeff*(1/2pi)
  !                  \int  M2_equil(r,theta) dtheta
  !       with press_coef = 2/3 if Nmu.ne.0        
  !                       = 2   otherwise
  !    . T_eq(r)   = P_eq(r)/n_eq(r)        
  !    . nbions_eq = \int J0.feq Js Jv dr dtheta d3v
  !                = \int n_eq dr                   
  !-------------------------------------------------------------
!baoter (ATTENTION OTER le calcul de M2equil)
!eaoter
  subroutine compute_equilibrium_part(geom,J0, &
    init_prof,init_magnet,init_curr,Sfmu_eq,SJ0fmu_eq, &
    M2eq_rtheta,nGieq_rtheta,neq_r,dneq_dr,Tieq_r,Snbions_eq)
    use globals, only : pglobal_id, Nbproc_mu, mu_id, &
      mpi_comm_intermu, min_err_rbar, max_err_rbar, &
      istart, iend, jstart, jend, ipow, a, alphaq, &
      Rarray1_NrNtheta, Rarray2_NrNtheta, &
      Rarray3_NrNtheta, Rarray4_NrNtheta
    use utils_module, only : deriv1
    use gyroaverage_class
    include "mpiperso.h"
    type(geometry)               , intent(in)    :: geom
    type(J0operator)             , intent(inout) :: J0
    type(init_profile)           , intent(in)    :: init_prof
    type(init_magnetic)          , intent(in)    :: init_magnet
    type(init_current)           , intent(in)    :: init_curr
    real(RKIND), dimension(:,:,:), pointer       :: Sfmu_eq
    real(RKIND), dimension(:,:,:), pointer       :: SJ0fmu_eq
    real(RKIND), dimension(0:,0:), intent(out)   :: M2eq_rtheta
    real(RKIND), dimension(0:,0:), intent(out)   :: nGieq_rtheta 
    real(RKIND), dimension(0:)   , intent(out)   :: neq_r 
    real(RKIND), dimension(0:)   , intent(out)   :: dneq_dr
    real(RKIND), dimension(0:)   , intent(out)   :: Tieq_r  
    real(RKIND)                  , intent(out)   :: Snbions_eq
      
    integer     :: ierr
    integer     :: ir, itheta, ivpar
    real(RKIND) :: fmu_eq_tmp , J0fmu_eq_tmp 
    real(RKIND) :: vparl2, B_ij, energy
    real(RKIND) :: intdtheta_Js_tmp
    real(RKIND) :: jacob_space_tmp, jacob_vel_tmp
    real(RKIND) :: coeff_intdr_tmp, coeff_intdtheta_tmp
    real(RKIND) :: coeff_intdvpar_tmp, coeff_intdmu_tmp
    real(RKIND) :: neq_tmp, Peq_tmp, press_coeff
      
    real(RKIND), dimension(:,:), pointer :: feq_rtheta
    real(RKIND), dimension(:,:), pointer :: J0feq_rtheta
    real(RKIND), dimension(:,:), pointer :: nGieq_rtheta_loc
    real(RKIND), dimension(:,:), pointer :: M2eq_rtheta_loc
      
    if (pglobal_id.eq.0) then
      write(6,'(2A45)') &
        '---------------------------------------------', &
        '---------------------------------------------'
      write(6,*) '--> EQUILIBRIUM COMPUTATION :'
      if (alphaq.ne.2._RKIND) then        
        write(6,*) '-->   - WARNING : the q profile is ', &
          'not parabolic, so rbar is '
        write(6,*) '-->               computed with an ', &
          'approximation at first order of psibar'
      else
        write(6,*) '-->   - rbar is the analytical ',&
          'expression of psibar'
      end if
    end if
      
    feq_rtheta       => Rarray1_NrNtheta
    J0feq_rtheta     => Rarray2_NrNtheta
    nGieq_rtheta_loc => Rarray3_NrNtheta
    M2eq_rtheta_loc  => Rarray4_NrNtheta
      
    !*** initialisation to 0 ***
    do itheta = 0,geom%Ntheta
      do ir = 0,geom%Nr
        nGieq_rtheta_loc(ir,itheta) = 0._RKIND
        M2eq_rtheta_loc(ir,itheta)  = 0._RKIND
      end do
    end do
      
    min_err_rbar     = a
    max_err_rbar     = 0._RKIND
    coeff_intdmu_tmp = geom%coeff_intdmu(mu_id)
      
    do ivpar = 0,geom%Nvpar
      vparl2             = geom%vparg(ivpar)**2
      coeff_intdvpar_tmp = geom%coeff_intdvpar(ivpar)
      !*** computation of feq ***
      do itheta = 0,geom%Ntheta
        do ir = 0,geom%Nr
          call compute_feq_value(geom,ir,itheta,ivpar,mu_id, &
            init_prof,init_magnet,fmu_eq_tmp)
          feq_rtheta(ir,itheta) = fmu_eq_tmp
        end do
      end do
      !*** save fmu_eq(istart:iend,jstart:jend,0:Nvpar) ***
      do itheta = jstart_fmu,jend_fmu
        do ir = istart_fmu,iend_fmu
          Sfmu_eq(ir,itheta,ivpar) = feq_rtheta(ir,itheta)
        end do
      end do
      !*** computation of J0.feq ***
      call compute_gyrofunction_2D(J0,geom,mu_id, &
        feq_rtheta,J0feq_rtheta)
      do itheta = jstart,jend
        do ir = istart,iend
          SJ0fmu_eq(ir,itheta,ivpar) = J0feq_rtheta(ir,itheta)
        end do
      end do
      !*** computation of nGieq(r,theta) ***
      do itheta = 0,geom%Ntheta
        do ir = 0,geom%Nr
          !--> J0.feq
          J0fmu_eq_tmp = J0feq_rtheta(ir,itheta)
          !-> jacobian_velocity = 2*pi*Bstar
          call compute_jacobian_velocity(geom,init_magnet, &
            init_curr,ir,itheta,ivpar,jacob_vel_tmp)
          !--> \int Jv d3v J0.feq
          nGieq_rtheta_loc(ir,itheta) = &
            nGieq_rtheta_loc(ir,itheta) + &
            J0fmu_eq_tmp*jacob_vel_tmp * &
            coeff_intdvpar_tmp*coeff_intdmu_tmp
          !--> \int[feq*(0.5*vpar^2+mu*B(r,theta))]dvpar
          B_ij   = init_magnet%B_norm(ir,itheta)
          energy = 0.5_RKIND*vparl2 + geom%mug(mu_id)*B_ij
          M2eq_rtheta_loc(ir,itheta) = &
            M2eq_rtheta_loc(ir,itheta) + &
            feq_rtheta(ir,itheta)*energy*jacob_vel_tmp * &
            coeff_intdvpar_tmp*coeff_intdmu_tmp
        end do
      end do
    end do
      
    !*** integration in mu direction for computing ***
    !***  \int Jv d3v J0(sqrt(2*mu)).feq           ***
    if (Nbproc_mu.ne.1) then
      call MPI_ALLREDUCE(nGieq_rtheta_loc,nGieq_rtheta, &
        (geom%Nr+1)*(geom%Ntheta+1),MPI_REAL8,MPI_SUM, &
        mpi_comm_intermu,ierr)
      call MPI_ALLREDUCE(M2eq_rtheta_loc,M2eq_rtheta, &
        (geom%Nr+1)*(geom%Ntheta+1), &
        MPI_REAL8,MPI_SUM,mpi_comm_intermu,ierr)
    else
      do itheta = 0,geom%Ntheta
        do ir = 0,geom%Nr
          nGieq_rtheta(ir,itheta) = nGieq_rtheta_loc(ir,itheta)
          M2eq_rtheta(ir,itheta)  = M2eq_rtheta_loc(ir,itheta)
        end do
      end do
    end if
      
    !***********************************************************
    !*** computation of the radial equilibrium profiles:     ***
    !***   . n_eq(r) = ((1/2pi)\int  nGi_eq(r,theta) dtheta  ***
    !***   . P_eq(r) = press_coeff*(1/2pi)                   ***
    !***               \int  M2_equil(r,theta) dtheta        ***
    !***       with press_coef = 2/3 if Nmu.ne.0             ***
    !***                       = 2   otherwise               ***
    !***   . T_eq(r) = P_eq(r)/n_eq(r)                       ***
    !*** computation of                                      ***
    !***   . nbions_eq = \int J0.feq Js Jv dr dtheta d3v     *** 
    !***               = \int n_eq dr                        ***
    !***********************************************************
    press_coeff = 0.5_RKIND+float(ipow)
    Snbions_eq  = 0._RKIND
    do ir = 0,geom%Nr
      neq_tmp         = 0._RKIND
      Peq_tmp         = 0._RKIND
      coeff_intdr_tmp = geom%coeff_intdr(ir)
      do itheta = 0,geom%Ntheta
        coeff_intdtheta_tmp = geom%coeff_intdtheta(itheta)
        jacob_space_tmp     = jacobian_space(ir,itheta)
        neq_tmp             = neq_tmp + &
          nGieq_rtheta(ir,itheta) * &
          jacob_space_tmp*coeff_intdtheta_tmp
        Peq_tmp             = Peq_tmp + &
          M2eq_rtheta(ir,itheta) * &
          jacob_space_tmp*coeff_intdtheta_tmp
      end do
      Peq_tmp          = Peq_tmp/press_coeff
      intdtheta_Js_tmp = intdtheta_Js(ir)
      neq_r(ir)        = neq_tmp/intdtheta_Js_tmp
      Peq_tmp          = Peq_tmp/intdtheta_Js(ir)
      Tieq_r(ir)       = Peq_tmp/neq_r(ir)
      Snbions_eq       = Snbions_eq + &
        neq_tmp*coeff_intdr_tmp*geom%Lphi
    end do
      
    !*** dn_eq/dr ***
    call deriv1(neq_r,dneq_dr,geom%Nr,geom%dr,0)
    
    if (pglobal_id.eq.0) then
      print*,' min err(r-rbar) = ', min_err_rbar
      print*,' max err(r-rbar) = ', max_err_rbar
      write(6,'(2A45)') &
        '---------------------------------------------', &
        '---------------------------------------------'
    end if
  end subroutine compute_equilibrium_part
end module fequil4d_module
