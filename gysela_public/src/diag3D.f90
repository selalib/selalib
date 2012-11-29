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
! file : diag3D.f90
! date : 18/04/2005
!  Computes and save the fluid moments of the distribution
!   function of particles:
!     -> J0M0 = \int J0.f Jv dvpar dmu  
!     -> J0M1 = \int J0.f vpar Jv dvpar dmu 
!     -> J0M2 = \int J0.f (0.5*vpar^2+mu B) Jv dvpar dmu 
!  
!  Computes and save the fluid moments of the distribution
!   function of guiding-centers:
!     -> M2_parallel = \int f*(0.5*vpar^2) Jv dvpar dmu 
!     -> M2_perp     = \int f*(mu*B) Jv dvpar dmu 
!
!  Rk: Jv is the jacobian in velocity, i.e : 
!   Jv = 2*pi*B(r,theta)
!
!  Save the electric potential Phi(r,theta,phi)
!--------------------------------------------------------------- 
module diag3D_module
  use prec_const
      
  implicit none
      
  !--> fmin and fmax need for CalviExport
  real(RKIND), public :: fmin_muid_diag
  real(RKIND), public :: fmax_muid_diag
  
  include "Bstar_inline.h"
  include "velocities_inline.h"
  !******************************
  contains
  !******************************
#include "Bstar_inline.f90"
#include "velocities_inline.f90"
      
  !--------------------------------------------------------
  ! Compute and save 3D fluid moments:
  !  - J0M0
  !  - J0M1
  !  - J0M2
  !  - M2_parallel
  !  - M2_perp
  ! Save 3D electric potential Phi
  !--------------------------------------------------------
  subroutine compute_save_3D(geom,init_prof,init_magnet, &
    init_curr,J0,f,sJ0ftmp,sPhi)
    use globals
    use OMPutils_module
    use clock_module
    use integration_module, only : compute_vpar_integral_colloc, &
      compute_omp_vpar_integral_CS
    use gyroaverage_class
    use poisson_class
    use fdistribu5d_class
    use init_magnetic_class
    use init_current_class
#ifndef NOHDF5
    use resu3D_saving_module
    use f5D_saving_module
#endif
    include "mpiperso.h"
    type(geometry)     , intent(in)    :: geom
    type(init_profile) , intent(in)    :: init_prof
    type(init_magnetic), intent(in)    :: init_magnet
    type(init_current) , intent(in)    :: init_curr
    type(J0operator)   , intent(inout) :: J0
    type(fdistribu5d)  , intent(inout) :: f
    type(fdistribu5d)  , intent(inout) :: sJ0ftmp
    real(RKIND), &
      dimension(:,:,:) , pointer       :: sPhi
    integer, parameter :: nbmoments = 8
    integer     :: ir, itheta, iphi, ivpar, imu, iphiloc
    real(RKIND) :: mu_value
    real(RKIND) :: ri, Bij
    real(RKIND) :: vpar, fval, J0fval, energy
    real(RKIND) :: coeff_intdmu
    real(RKIND) :: jacob_space_ij, jacob_vel_tmp
      
    real(RKIND),    dimension(:)    , pointer :: scoefvpar
    real(RKIND),    dimension(:)    , pointer :: fvperp2
    real(RKIND),    dimension(:)    , pointer :: fvpar2
    real(RKIND),    dimension(:)    , pointer :: J0f
    real(RKIND),    dimension(:)    , pointer :: J0fv
    real(RKIND),    dimension(:)    , pointer :: J0fE
      
    real(RKIND),    dimension(:,:)  , pointer :: intfvperp2_dvpar
    real(RKIND),    dimension(:,:)  , pointer :: intfvpar2_dvpar
    real(RKIND),    dimension(:,:)  , pointer :: intJ0f_dvpar
    real(RKIND),    dimension(:,:)  , pointer :: intJ0fv_dvpar
    real(RKIND),    dimension(:,:)  , pointer :: intJ0fE_dvpar
    real(RKIND),    dimension(:,:)  , pointer :: J0f_rtheta
    complex(CKIND), dimension(:,:)  , pointer :: Acomp
    real(RKIND),    dimension(:,:,:), pointer :: M2par_loc
    real(RKIND),    dimension(:,:,:), pointer :: M2perp_loc
    real(RKIND),    dimension(:,:,:), pointer :: J0M0_loc
    real(RKIND),    dimension(:,:,:), pointer :: J0M1_loc
    real(RKIND),    dimension(:,:,:), pointer :: J0M2_loc
    real(RKIND),    dimension(:,:,:), pointer :: target3D
      
    ! -> variable for parallelization
    integer :: no_operation
    integer :: ierr, mypid
    integer :: itag_phi, j
    integer :: i, nbreqr, local_id, dep, base
    integer :: global_id, startbuf
    integer :: l_istart, l_iend, l_jstart, l_jend
    integer :: nbelements, tid
      
    ! -> variables for time perfomance analysis
    integer(TIMEPREC) :: tdeb, tfin
    real(RKIND)       :: tcomm1, tcomm2, tsave1
    real(RKIND)       :: cexport, tivpar, timu
    
    
#ifndef NOHDF5
    call clck_time(tdeb)
    call clck_time(bclock_diag3D)
      
    call MPI_COMM_RANK(mpi_comm_moments,mypid,ierr)
    if ((FMoments3D_saving).and. &
      (mod(iter_glob,FMoments3D_nbstep).eq.0)) then
      if (pglobal_id.eq.outputproc) &
        write (uout_res,*) '---> saving of 3D moments'
      J0M0_loc   => Rarray1_NrNthetaDomphi
      J0M1_loc   => Rarray2_NrNthetaDomphi
      J0M2_loc   => Rarray3_NrNthetaDomphi
      M2perp_loc => Rarray4_NrNthetaDomphi
      M2par_loc  => Rarray5_NrNthetaDomphi
      target3D   => Rarray1_NrNthetaNphi
      
      tid        = 1
      Acomp      => Comp1_1Nrp1_1Ntheta(tid)%val
      J0f_rtheta => Romp1_1Nrp1_1Nthetap1(tid)%val
      
      !*** initialisation to 0 of the integrals ***
      call pptranspose_forward(f%values)
      do ivpar = lstart,lend
        do iphi = kstart,kend
          do itheta = 0,f%n2
            do ir = 0,f%n1
              J0f_rtheta(ir+1,itheta+1) = &
                f4D_transp(ir,itheta,iphi,ivpar)
            end do
          end do
          !--> computation of J0.f
          call omp_compute_gyrofunction_2D(J0,geom, &
            mu_id,Acomp,J0f_rtheta)
          do itheta = 0,f%n2
            do ir = 0,f%n1
              f4D_transp(ir,itheta,iphi,ivpar) = &
                J0f_rtheta(ir+1,itheta+1)
            end do
          end do
        end do
      end do
      !->  f4D_transp -> J0ftmp%values
      call pptranspose_backward(sJ0ftmp%values)
      
      mu_value  = f%mu_value
      
      !*** Parallel communication to prepare the      ***
      !***  working processors                        ***
      !***  (i.e processors having a phi value)       *** 
      !***  to receive the partial values of          ***
      !***  several moments of f                      ***
      !*** (RK: The processors are always prepared    ***
      !***   to the receiving before the sending)     ***
#ifdef _OPENMP
!$OMP PARALLEL private(tid,ir,itheta,iphi,ivpar, &
!$OMP scoefvpar,fvperp2,fvpar2,J0f,J0fv,J0fE, &
!$OMP ri,Bij,vpar,energy,fval,J0fval,jacob_space_ij, &
!$OMP jacob_vel_tmp) &
!$OMP default(shared)
!$OMP BARRIER
      tid = 1+omp_get_thread_num()
#else
      tid = 1
#endif
      
      scoefvpar => Romp_scoefvpar(tid)%val
      J0f       => Romp1_0Nvpar(tid)%val
      J0fv      => Romp2_0Nvpar(tid)%val
      J0fE      => Romp3_0Nvpar(tid)%val
      fvpar2    => Romp4_0Nvpar(tid)%val
      fvperp2   => Romp5_0Nvpar(tid)%val
      
!$OMP DO SCHEDULE(STATIC)
      do iphi = 0,f%n3-1
        do itheta = f%jstart,f%jend
          do ir = f%istart,f%iend
            Bij            = init_magnet%B_norm(ir,itheta)
            ri             = geom%rg(ir)
            jacob_space_ij = jacobian_space(ir,itheta)
            do ivpar = 0,f%n4
              vpar           = geom%vparg(ivpar)
              !-> jacobian_velocity = 2*pi*Bstar
              call compute_jacobian_velocity(geom,init_magnet, &
                init_curr,ir,itheta,ivpar,jacob_vel_tmp)
              energy         = 0.5_RKIND*vpar*vpar+mu_value*Bij
              fval           = f%values(ir,itheta,iphi,ivpar)
              J0fval         = sJ0ftmp%values(ir,itheta,iphi,ivpar)
              J0f(ivpar)     = J0fval*jacob_vel_tmp
              J0fv(ivpar)    = J0fval*vpar*jacob_vel_tmp
              J0fE(ivpar)    = J0fval*energy*jacob_vel_tmp
              fvpar2(ivpar)  = fval*(0.5_RKIND*vpar*vpar) * &
                jacob_vel_tmp
              fvperp2(ivpar) = fval*(mu_value*Bij)*jacob_vel_tmp
            end do
      
            !*** computation of the integrals in vpar direction ***
            !--> \int J0.f Jv dvpar 
            call compute_vpar_integral_colloc(J0f(0:f%n4),geom, &
              Rarray_PNrPNthetaNphi_nbM8(0,ir,itheta,iphi))
            !--> \int (J0.f)*vpar Jv dvpar
            call compute_vpar_integral_colloc(J0fv(0:f%n4),geom, &
              Rarray_PNrPNthetaNphi_nbM8(1,ir,itheta,iphi))
            !--> \int J0.f*(0.5*vpar^2+mu*B) Jv dvpar
            call compute_vpar_integral_colloc(J0fE(0:f%n4),geom, &
              Rarray_PNrPNthetaNphi_nbM8(2,ir,itheta,iphi))
            !--> \int f*(0.5*vpar^2) Jv dvpar
            call compute_vpar_integral_colloc(fvpar2(0:f%n4),geom, &
              Rarray_PNrPNthetaNphi_nbM8(3,ir,itheta,iphi))
            !--> \int f*(mu*B) Jv dvpar
            call compute_vpar_integral_colloc(fvperp2(0:f%n4), &
              geom,Rarray_PNrPNthetaNphi_nbM8(4,ir,itheta,iphi))
          enddo
        enddo
      end do
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL
      call clck_time(tfin)
      call clck_ldiff(tdeb,tfin,tivpar)
      
      call ppbarrier()
      call clck_time(tdeb)
      nbreqr = 0
      if (iphistart .ne. -1) then
        do imu = 0, Nmu
          iphiloc = 0
          do iphi = 0,f%n3-1
            if (moments_mapphi(iphi) .eq. pglobal_id) then
              do local_id = 0,Nbproc_loc-1
                global_id = imu*Nbproc_loc + local_id
                startbuf  = dom_theta * &
                  (local_id+Nbproc_loc*(imu+(Nmu+1)*iphiloc)) 
                nbelements = nbmoments * dom_r * dom_theta
                itag_phi   = 20 + iphi
                call MPI_IRECV( &
                  Rarray_NrNthetamuphi_nbM8(0,startbuf), &
                  nbelements,MPI_REAL8,global_id,itag_phi, &
                  MPI_COMM_WORLD,moments_reqr(nbreqr),ierr)
                nbreqr = nbreqr + 1
              enddo
              iphiloc = iphiloc + 1
            endif
          enddo
        enddo
      endif
      call clck_time(tfin)
      call clck_ldiff(tdeb,tfin,tcomm1)
      
      !*** Parallel communication to send the five  ***
      !***  integrals in vpar                       ***
      !***  [ \int f*(vpar^2) Jv dvpar,             ***
      !***    \int f*(mu*B) Jv dvpar,               ***
      !***    \int (J0.f) Jv dvpar,                 ***
      !***    \int (J0.f)*vpar Jv dvpar,            ***
      !***    \int (J0.f)*(vpar^2+mu*B) Jv dvpar,   ***
      !*** The parallelism is such that             ***
      !***   (r=block,theta=block,phi=*)            ***
      !*** RK : These 5 integrals are saved in an   ***
      !***      unique array                        ***
      call clck_time(tdeb)
      nbelements = nbmoments*dom_r*dom_theta
      do iphi = 0,f%n3-1
        itag_phi = 20 + iphi
        call MPI_SEND( &
          Rarray_PNrPNthetaNphi_nbM8(0,istart,jstart,iphi), &
          nbelements,MPI_REAL8,moments_mapphi(iphi),itag_phi, &
          MPI_COMM_WORLD,ierr)
      enddo
      if (nbreqr .ne. 0) then
        call MPI_WAITALL(nbreqr,moments_reqr(0), &
          moments_status(1,0),ierr)
      end if
      call clck_time(tfin)
      call clck_ldiff(tdeb,tfin,tcomm2)
      
      !*** Reconstruction of the complete (r,theta) domain    ***
      !***  RK:- Parallelization in phi on all                ***
      !***       the processors                               ***
      !***     - iphistart.eq.-1 for some processors          ***
      !***       if the total number of processors is         ***
      !***       greater than the number of phi values (Nphi) ***
      call clck_time(tdeb)
      timu = 0.
      
      if (iphistart .ne. -1) then
#ifdef _OPENMP
!$OMP PARALLEL private(tid,iphi,itheta,ir,ivpar, &
!$OMP iphiloc,imu,coeff_intdmu,local_id,base,dep, &
!$OMP l_iend,l_istart,startbuf,i,j,l_jstart,l_jend,Acomp, &
!$OMP intJ0f_dvpar,intJ0fv_dvpar,intJ0fE_dvpar, &
!$OMP intfvpar2_dvpar,intfvperp2_dvpar) default(shared) 
!$OMP BARRIER
        tid = 1+omp_get_thread_num()
#else
        tid = 1
#endif
        Acomp            => Comp1_1Nrp1_1Ntheta(tid)%val
        intJ0f_dvpar     => Romp1_1Nrp1_1Nthetap1(tid)%val
        intJ0fv_dvpar    => Romp2_1Nrp1_1Nthetap1(tid)%val
        intJ0fE_dvpar    => Romp3_1Nrp1_1Nthetap1(tid)%val
        intfvpar2_dvpar  => Romp4_1Nrp1_1Nthetap1(tid)%val
        intfvperp2_dvpar => Romp5_1Nrp1_1Nthetap1(tid)%val
      
        do iphi = 0,geom%Nphi-1
          if (moments_mapphi(iphi) .eq. pglobal_id) then
            iphiloc = iphi - iphistart
            do imu = 0,Nmu
              ! -> This if block replaces a OMP DO 
              !     SCHEDULE(static,1) directive
              !  (this for avoiding problem on BULL architecture)
              if (tid-1 .eq. mod(imu,Nbthread)) then
                coeff_intdmu = geom%coeff_intdmu(imu)
                do itheta = 0,f%n2
                  do ir = 0,f%n1
                    intJ0f_dvpar(1+ir,1+itheta)     = 0._RKIND
                    intJ0fv_dvpar(1+ir,1+itheta)    = 0._RKIND
                    intJ0fE_dvpar(1+ir,1+itheta)    = 0._RKIND
                    intfvpar2_dvpar(1+ir,1+itheta)  = 0._RKIND
                    intfvperp2_dvpar(1+ir,1+itheta) = 0._RKIND
                  end do
                end do
                do local_id = 0,Nbproc_loc-1
                  base     = (local_id/Nbproc_r)
                  dep      = mod(local_id,Nbproc_r)              
                  l_istart = dep *  dom_r
                  l_iend   = l_istart + dom_r - 1
                  l_jstart = base * dom_theta
                  l_jend   = l_jstart + dom_theta - 1
                  startbuf = dom_theta * (local_id + &
                    Nbproc_loc * (imu + (Nmu+1) * iphiloc)) 
                  do itheta = l_jstart, l_jend
                    do ir = l_istart, l_iend
                      intJ0f_dvpar(1+ir,1+itheta)     = &
                        Rarray_NrNthetamuphi_nbM8( &
                        nbmoments*(ir-l_istart)+0, &
                        startbuf+itheta-l_jstart) 
                      intJ0fv_dvpar(1+ir,1+itheta)    = &
                        Rarray_NrNthetamuphi_nbM8( &
                        nbmoments*(ir-l_istart)+1, &
                        startbuf+itheta-l_jstart) 
                      intJ0fE_dvpar(1+ir,1+itheta)    = &
                        Rarray_NrNthetamuphi_nbM8( &
                        nbmoments*(ir-l_istart)+2, &
                        startbuf+itheta-l_jstart)
                      intfvpar2_dvpar(1+ir,1+itheta)  = &
                        Rarray_NrNthetamuphi_nbM8( &
                        nbmoments*(ir-l_istart)+3, &
                        startbuf+itheta-l_jstart) 
                      intfvperp2_dvpar(1+ir,1+itheta) = &
                        Rarray_NrNthetamuphi_nbM8( &
                        nbmoments*(ir-l_istart)+4, &
                        startbuf+itheta-l_jstart) 
                    enddo
                  enddo
                enddo
                !*** periodic boundary conditions in theta ***
                do ir = 0,f%n1
                  intJ0f_dvpar(1+ir,1+f%n2)     = &
                    intJ0f_dvpar (1+ir,1+0)
                  intJ0fv_dvpar(1+ir,1+f%n2)    = &
                    intJ0fv_dvpar(1+ir,1+0)
                  intJ0fE_dvpar(1+ir,1+f%n2)    = &
                    intJ0fE_dvpar(1+ir,1+0)
                  intfvpar2_dvpar(1+ir,1+f%n2)  = &
                    intfvpar2_dvpar (1+ir,1+0)
                  intfvperp2_dvpar(1+ir,1+f%n2) = &
                    intfvperp2_dvpar (1+ir,1+0)
                end do
      
                !*** computation of the integrals in mu ***
                do itheta = 0,f%n2-1
                  do ir = 0,f%n1
                    !-> J0M0 = \int (\int J0.f Jv dvpar) dmu
                    Romp9_1Nrp1_1Nthetap1(imu)%val(1+ir,1+itheta) &
                      = intJ0f_dvpar(1+ir,1+itheta)*coeff_intdmu 
                    !-> J0M1 = \int (\int J0.f*vpar Jv dvpar) dmu
                    Romp10_1Nrp1_1Nthetap1(imu)%val(1+ir,1+itheta) &
                      = intJ0fv_dvpar(1+ir,1+itheta)*coeff_intdmu
                    !-> J0M2 = \int (\int 
                    !          J0.f*(0.5*vpar^2+mu*B) Jv dvpar) dmu 
                    Romp11_1Nrp1_1Nthetap1(imu)%val(1+ir,1+itheta) &
                      = intJ0fE_dvpar(1+ir,1+itheta)*coeff_intdmu
                    !-> M2vpar = \int \int 
                    !            f*(0.5*vpar^2) Jv dvpar dmu 
                    Romp12_1Nrp1_1Nthetap1(imu)%val(1+ir,1+itheta) &
                      = intfvpar2_dvpar(1+ir,1+itheta)*coeff_intdmu
                    !-> M2perp = \int \int f*(mu*B) Jv dvpar dmu 
                    Romp13_1Nrp1_1Nthetap1(imu)%val(1+ir,1+itheta) &
                      = intfvperp2_dvpar(1+ir,1+itheta)*coeff_intdmu
                  end do
                end do
              end if
            end do
!$OMP BARRIER
!$OMP MASTER 
            do itheta = 0,geom%Ntheta
              do ir = 0,geom%Nr
                J0M0_loc(ir,itheta,iphi)   = 0._RKIND
                J0M1_loc(ir,itheta,iphi)   = 0._RKIND
                J0M2_loc(ir,itheta,iphi)   = 0._RKIND
                M2par_loc(ir,itheta,iphi)  = 0._RKIND
                M2perp_loc(ir,itheta,iphi) = 0._RKIND
              end do
            end do
      
            !*** computation of the integrals in mu ***
            do j = 0,Nmu
              do itheta = 0,f%n2-1
                do ir = 0,f%n1
                  !-> J0M0 = \int (\int J0.f Jv dvpar) dmu
                  J0M0_loc(ir,itheta,iphi) = &
                    J0M0_loc(ir,itheta,iphi) + &
                    Romp9_1Nrp1_1Nthetap1(j)%val(1+ir,1+itheta)
                  !-> J0M1 = \int (\int J0.f*vpar Jv dvpar) dmu
                  J0M1_loc(ir,itheta,iphi) = &
                    J0M1_loc(ir,itheta,iphi) + &
                    Romp10_1Nrp1_1Nthetap1(j)%val(1+ir,1+itheta)
                  ! -> J0M2 = \int [\int J0.f *       
                  !    (0.5*vpar^2+mu*B) Jv dvpar] dmu 
                  J0M2_loc(ir,itheta,iphi) = &
                    J0M2_loc(ir,itheta,iphi) + &
                    Romp11_1Nrp1_1Nthetap1(j)%val(1+ir,1+itheta)
                  !-> M2_par = \int [\int f * 
                  !                  (0.5*vpar^2) Jv dvpar] dmu
                  M2par_loc(ir,itheta,iphi) = &
                    M2par_loc(ir,itheta,iphi) + &
                    Romp12_1Nrp1_1Nthetap1(j)%val(1+ir,1+itheta)
                  !-> M2_perp = \int [\int f * (mu*B) Jv dvpar] dmu
                  M2perp_loc(ir,itheta,iphi) = &
                    M2perp_loc(ir,itheta,iphi) + &
                    Romp13_1Nrp1_1Nthetap1(j)%val(1+ir,1+itheta)
                end do
              end do
            end do
      
            !*** periodic conditions in theta ***
            do ir = 0,f%n1
              J0M0_loc(ir,f%n2,iphi)   = J0M0_loc(ir,0,iphi)
              J0M1_loc(ir,f%n2,iphi)   = J0M1_loc(ir,0,iphi)
              J0M2_loc(ir,f%n2,iphi)   = J0M2_loc(ir,0,iphi)
              M2par_loc(ir,f%n2,iphi)  = M2par_loc(ir,0,iphi)
              M2perp_loc(ir,f%n2,iphi) = M2perp_loc(ir,0,iphi)
            end do
!$OMP END MASTER
          endif
          no_operation = 0
!$OMP BARRIER
        end do
!$OMP END PARALLEL
        call clck_time(tfin)
        call clck_ldiff(tdeb,tfin,timu)
        
        !*** distribution of the six moments  ***
        !***   J0M0(r=*,theta=*,phi=*), J0M1, ***
        !***   J0M2, M2perp, M2vpar           ***
        call clck_time(tdeb)
        nbelements = dom_mapphi*(Nr+1)*(Ntheta+1)
      
        if (diag_para) then
!bmodif 
!VG!      call MPI_GATHER(J0M0_loc(0,0,iphistart), &
!VG!        nbelements,MPI_REAL8,target3D(0,0,0),nbelements, &
!VG!        MPI_REAL8,1,mpi_comm_moments, ierr)
!VG!      call MPI_GATHER(J0M1_loc(0,0,iphistart), &
!VG!        nbelements,MPI_REAL8,target3D(0,0,0),nbelements, &
!VG!        MPI_REAL8,2,mpi_comm_moments, ierr)
!VG!      call MPI_GATHER(J0M2_loc(0,0,iphistart), &
!VG!        nbelements,MPI_REAL8,target3D(0,0,0),nbelements, &
!VG!        MPI_REAL8,3,mpi_comm_moments, ierr)
!VG!      call MPI_GATHER(M2perp_loc(0,0,iphistart), &
!VG!        nbelements,MPI_REAL8,target3D(0,0,0),nbelements, &
!VG!        MPI_REAL8,4,mpi_comm_moments, ierr)
!VG!      call MPI_GATHER(M2par_loc(0,0,iphistart), &
!VG!        nbelements,MPI_REAL8,target3D(0,0,0),nbelements, &
!VG!        MPI_REAL8,5,mpi_comm_moments, ierr)
      
          ! -> solution to avoid temporary the problem 
          !    of MPI_GATHER instruction
          call comm_gather3D(J0M0_loc,nbelements, &
            iphistart,dom_mapphi,1, &
            target3D,mpi_comm_moments) 
          call comm_gather3D(J0M1_loc,nbelements, &
            iphistart,dom_mapphi,2, &
            target3D,mpi_comm_moments) 
          call comm_gather3D(J0M2_loc,nbelements, &
            iphistart,dom_mapphi,3, &
            target3D,mpi_comm_moments) 
          call comm_gather3D(M2par_loc,nbelements, &
            iphistart,dom_mapphi,4, &
            target3D,mpi_comm_moments)
          call comm_gather3D(M2perp_loc,nbelements, &
            iphistart,dom_mapphi,5, &
            target3D,mpi_comm_moments) 
!emodif
      
          !*** periodic conditions in phi ***
          do itheta = 0,f%n2
            do ir = 0,f%n1
              target3D(ir,itheta,f%n3) = target3D(ir,itheta,0)
            end do
          end do
          select case (mypid)
          case(1) 
            call HDF5_resu3D_saving('J0M0',FMoments3D_nb_diag, &
              geom,target3D)
          case(2) 
            call HDF5_resu3D_saving('J0M1',FMoments3D_nb_diag, &
              geom,target3D)
          case(3) 
            call HDF5_resu3D_saving('J0M2',FMoments3D_nb_diag, &
              geom,target3D)
          case(4) 
            call HDF5_resu3D_saving('M2par',FMoments3D_nb_diag, &
              geom,target3D)
          case(5) 
            call HDF5_resu3D_saving('M2perp',FMoments3D_nb_diag, &
              geom,target3D)
          end select
        else
          !*** saving of J0M0 ***
          call comm_gather3D(J0M0_loc,nbelements, &
            iphistart,dom_mapphi,0,target3D,mpi_comm_moments) 
          !--> periodic conditions in phi
          target3D(0:f%n1,0:f%n2,f%n3) = target3D(0:f%n1,0:f%n2,0)
          if (mypid .eq. 0) &
            call HDF5_resu3D_saving('J0M0',FMoments3D_nb_diag, &
            geom,target3D)
      
          !*** saving of J0M1 ***
          call comm_gather3D(J0M1_loc,nbelements, &
            iphistart,dom_mapphi,0,target3D,mpi_comm_moments) 
          !--> periodic conditions in phi
          target3D(0:f%n1,0:f%n2,f%n3) = target3D(0:f%n1,0:f%n2,0)
          if (mypid .eq. 0) &
            call HDF5_resu3D_saving('J0M1',FMoments3D_nb_diag, &
            geom,target3D)
      
          !*** saving of J0M2 ***
          call comm_gather3D(J0M2_loc,nbelements, &
            iphistart,dom_mapphi,0,target3D,mpi_comm_moments) 
          !--> periodic conditions in phi
          target3D(0:f%n1,0:f%n2,f%n3) = target3D(0:f%n1,0:f%n2,0)
          if (mypid .eq. 0) &
            call HDF5_resu3D_saving('J0M2',FMoments3D_nb_diag, &
            geom,target3D)
      
          !*** saving of M2_parallel ***
          call comm_gather3D(M2par_loc,nbelements, &
            iphistart,dom_mapphi,0,target3D,mpi_comm_moments)
          !--> periodic conditions in phi
          target3D(0:f%n1,0:f%n2,f%n3) = target3D(0:f%n1,0:f%n2,0)
          if (mypid .eq. 0) &
            call HDF5_resu3D_saving('M2par',FMoments3D_nb_diag, &
            geom,target3D)
      
          !*** saving of M2_perp ***
          call comm_gather3D(M2perp_loc,nbelements, &
            iphistart,dom_mapphi,0,target3D,mpi_comm_moments) 
          !--> periodic conditions in phi
          target3D(0:f%n1,0:f%n2,f%n3) = target3D(0:f%n1,0:f%n2,0)
          if (mypid .eq. 0) &
            call HDF5_resu3D_saving('M2perp',FMoments3D_nb_diag, &
            geom,target3D)
        endif
      end if
      FMoments3D_nb_diag = FMoments3D_nb_diag + 1
    end if
      
    !*** Saving of Phi3D ***
    if ((Phi3D_saving) .and. &
      (mod(iter_glob,Phi3D_nbstep).eq.0)) then
      if (pglobal_id.eq.outputproc) &
        write (uout_res,*) '---> saving of Phi3D'
      if ((mypid .eq. 0).and.(iphistart .ne. -1)) then
        call HDF5_resu3D_saving('Phi', &
          Phi3D_nb_diag,geom,sPhi)
      end if
      Phi3D_nb_diag = Phi3D_nb_diag + 1
    end if
      
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tsave1)
    call clck_time(tdeb)
      
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,cexport)
    
    call clck_time(eclock_diag3D)
    call clck_diff(bclock_diag3D,eclock_diag3D,&
      global_time_diag3D)
      
#ifdef TIMER    
    write(6,'(I4,A,F10.4,A,F10.4,A,&
      F10.4,A,F10.4,A,F10.4,A,F10.4)') pglobal_id, &
      " Temps comp_save_3D ivpar ", tivpar," ivperp ", &
      timu," comm1 ",tcomm1, " comm2 ", tcomm2, &
      " savef ",tsave1," cexport ",cexport
#endif
#endif
  end subroutine compute_save_3D
end module diag3D_module
