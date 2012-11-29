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
      
!-------------------------------------------------
! file : diagnostics.f90
! date : 17/07/2000
!  used for all physic output computation
!------------------------------------------------
module diagnostics_module
  use prec_const
  use MPIutils_module
  use mem_alloc_module
  use geometry_class
      
  implicit none
      
  include "Bstar_inline.h"
  include "velocities_inline.h"
  !******************************
  contains
  !******************************
#include "Bstar_inline.f90"
#include "velocities_inline.f90"
  
  !-----------------------------------------------------
  ! Used to write screen infos for ne and ni conservation
  !    at each iteration
  !-----------------------------------------------------
  subroutine WriteScreen(iter_time,dt,diag_status)
    use globals, only : plocal_id, mu_id, restart, &
      iter_glob, iter_run, nbiter, &
      nbelectrons_diag, nbions_diag, diag_targ, &
      uout_res, outputproc
    real(RKIND), intent(in) :: iter_time
    real(RKIND), intent(in) :: dt
    logical    , intent(in) :: diag_status
    
    integer :: uout
    integer :: line_count
      
    if ((restart).and.(iter_run.eq.1)) then
      line_count = 0
    else
      line_count = modulo(iter_run,40)
    end if
      
    if (pglobal_id.eq.outputproc) then
      if (line_count==0) then
        write(uout_res,'(3A36)') &
          '------------------------------------', &
          '------------------------------------', &
          '------------------------------------'
        write(uout_res,'(3A8,2A12,3A20)') 'GStep','LStep','NbStep',&
          'Time','TimeStep','ni','ne','ni/ne'
        write(6,'(3A36)') &
          '------------------------------------', &
          '------------------------------------', &
          '------------------------------------'
        write(6,'(3A8,2A12,3A20)') 'GStep','LStep','NbStep',&
          'Time','TimeStep','ni','ne','ni/ne'
      end if
      if (.not.(diag_status)) then
        write(uout_res,'(3I8,2(1pe12.3),3(1pe20.12))') &
          iter_glob, iter_run, nbiter, iter_time, dt
        write(6,'(3I8,2(1pe12.3),3(1pe20.12))') &
          iter_glob, iter_run, nbiter, iter_time, dt
      else
        write(uout_res,'(3I8,2(1pe12.3),3(1pe20.12))') &
          iter_glob, iter_run, nbiter, iter_time, dt, &
          nbions_diag, nbelectrons_diag, &
          (nbions_diag/nbelectrons_diag)
        write(6,'(3I8,2(1pe12.3),3(1pe20.12))') &
          iter_glob, iter_run, nbiter, iter_time, dt, &
          nbions_diag, nbelectrons_diag, &
          (nbions_diag/nbelectrons_diag)
      end if
    end if
  end subroutine WriteScreen
      
  !-----------------------------------------------------
  ! Used to write results on the conservation laws
  !    in gysela_CL.out file at each iteration
  !-----------------------------------------------------
  subroutine WriteResults_CL(iter_time)
    use globals             , only : pglobal_id, plocal_id, &
      mu_id, iter_glob, iter_run, nbiter, restart, &
      nbelectrons_diag, nbions_diag, diag_targ
    use physics_module, only : L2norm_diag, &
      entropy_diag, Enkin_diag, Enpot_diag
    use read_write_module    , only : uout_CL
    real(RKIND), intent(in) :: iter_time
      
    character(LEN=20) :: file_name_CL
      
    if (diag_targ(1) .eq. pglobal_id) then
      write(file_name_CL,'(A)')"gysela_CL.out"
      open(uout_CL, file = file_name_CL, status = 'UNKNOWN', &
        position = 'APPEND', form = 'FORMATTED')
      !-> title writting
      if (iter_glob.eq.-1) then
        write(uout_CL,'(A8,A12,7A20)') &
          'GStep','Time','ni','ne','ni/ne', &
          'Enkin','Enpot','Entropy','L2norm'
      end if
      !-> writting of the different 0D physical quantities
      write(uout_CL,'(I8,1pe12.3,7(1pe20.12))') &
        iter_glob, iter_time, nbions_diag, nbelectrons_diag, &
        nbions_diag/nbelectrons_diag, Enkin_diag, Enpot_diag, &
        entropy_diag, L2norm_diag
      !-> file closing
      close(uout_CL)
    end if
  end subroutine WriteResults_CL
      
  !-----------------------------------------------------
  ! Used to write results 
  !    in gysela_rprof.out file at each iteration
  !-----------------------------------------------------
  subroutine WriteResults_rprof(iter_time,geom,Phi)
    use globals, only : plocal_id, mu_id, diag_targ, &
      iter_glob, iter_run, nbiter, diag_level, Rarray1_Nr
    use physics_module    , only : Qr_turb_diag 
    use poisson_class     , only : Phi00_diag
    use read_write_module , only : uout_rprof
    use utils_module      , only : deriv1
    use integration_module, only : compute_intrdr_colloc
    real(RKIND)                        , intent(in) :: iter_time
    type(geometry)                     , intent(in) :: geom
    real(RKIND)   , dimension(0:,0:,0:), intent(in) :: Phi
      
    character(LEN=20)                  :: file_name_rprof
    integer                            :: ir, irp, itheta, iphi
    real(RKIND)                        :: Phi_tmp
    real(RKIND)                        :: intPhi2
    real(RKIND)                        :: sqrt_intPhi2
    real(RKIND)                        :: max_Phi00
    real(RKIND)                        :: max_dPhi00dr
    real(RKIND)                        :: int_Qrdr
    real(RKIND), dimension(:), pointer :: dPhi00dr 
!R3 #include "r3_info.h" !R3
    
    if (diag_level.ge.3) then
!R3 call r3_info_begin (r3_info_index_0, 'CPU_WriteResRprof') !R3
      if (diag_targ(3) .eq. pglobal_id) then
        !*** computation of                                   ***
        !***   sqrt(\int(Phi(rpeak,theta,phi)^2 dtheta dphi)) ***
        irp     = int(geom%Nr/2)
        intPhi2 = 0._RKIND
        do iphi = 0,geom%Nphi
          do itheta = 0,geom%Ntheta
            Phi_tmp = Phi(irp,itheta,iphi)
            intPhi2 = intPhi2 + Phi_tmp*Phi_tmp * &
              geom%coeff_intdtheta(itheta)*geom%coeff_intdphi(iphi)
          end do
        end do
        sqrt_intPhi2 = sqrt(intPhi2)
        !*** computation of the maximal value of ***
        !***  the (0,0) mode of Phi              ***
        max_Phi00    = maxval(Phi00_diag)
        !*** computation of maximal velocity due to zonal flow ***
        dPhi00dr => Rarray1_Nr
        call deriv1(Phi00_diag(0:),dPhi00dr(0:),geom%Nr,geom%dr,0)
        max_dPhi00dr = maxval(dPhi00dr)
        !*** computation of \int Q(r) r dr ***
        call compute_intrdr_colloc(Qr_turb_diag(0:),geom,int_Qrdr)
      
        !-> opening of the file 'gysela_rprof.out'
        write(file_name_rprof,'(A)')"gysela_rprof.out"
        open(uout_rprof, file = file_name_rprof, &
          status = 'UNKNOWN', position = 'APPEND', &
          form = 'FORMATTED')
        !-> title writting
        if (iter_glob.eq.-1) then
          write(uout_rprof,'(A8,A12,4A20)') &
            'GStep','Time','sqrt(intPhi2)', &
            'max(Phi00)', 'max(dPhi00dr)','\intQ(r)rdr'
        end if
        !-> writting of the different 0D physical quantities
        write(uout_rprof,'(I8,1pe12.3,4(1pe20.12))') &
          iter_glob, iter_time, sqrt_intPhi2, &
          max_Phi00, max_dPhi00dr, int_Qrdr
        !-> file closing
        close(uout_rprof)
      end if
!R3 call r3_info_end (r3_info_index_0) !R3
    end if
  end subroutine WriteResults_rprof
      
  !---------------------------------------------- 
  ! array saving for all the diagnostics in time 
  !---------------------------------------------- 
  subroutine diagnostics(iter_time,idiag_num,dt_diag, &
    geom,coord_sys,init_prof,init_magnet,init_curr, &
    fn,fnp1,sJ0ftmp,Sfmu_eq,sJ0fmu_eq,poiss,J0)
    use globals, only : uout_res, outputproc
    use coord_system_class
    use init_profile_class
    use init_magnetic_class
    use init_current_class
    use gyroaverage_class
    use fdistribu5d_class
    use poisson_class
    use output_saving_module
    use diag3D_module, only : compute_save_3D
    use Pcross_section_module
    use physics_module, only : compute_GC_moments, &
      compute_part_moments
#ifndef NOHDF5
    use resu3D_saving_module
    use f5D_saving_module
#endif
    use clock_module
!R3 #include "r3_info.h" !R3
    real(RKIND)        , intent(in)    :: iter_time
    integer            , intent(in)    :: idiag_num
    real(RKIND)        , intent(in)    :: dt_diag
    type(geometry)     , intent(in)    :: geom
    type(coord_system) , intent(in)    :: coord_sys
    type(init_profile) , intent(in)    :: init_prof
    type(init_magnetic), intent(in)    :: init_magnet
    type(init_current) , intent(in)    :: init_curr
    type(fdistribu5d)  , intent(inout) :: fn
    type(fdistribu5d)  , intent(inout) :: fnp1
    type(fdistribu5d)  , intent(inout) :: sJ0ftmp
    real(RKIND), &
      dimension(:,:,:) , pointer       :: Sfmu_eq
    real(RKIND), &
      dimension(:,:,:) , pointer       :: sJ0fmu_eq
    type(poisson)      , intent(inout) :: poiss
    type(J0operator)   , intent(inout) :: J0
    integer(TIMEPREC) :: tdeb,tfin
    real(RKIND)       :: tdiag(1:12)
    real(RKIND)       :: thdf5_diag(1:12)
!R3 call r3_info_begin (r3_info_index_0, 'CPU_diagnostics') !R3
      
    !*****************************************************
    !*** Parallelisation of the output saving          ***
    !***  - HDF5_CL_saving                             *** 
    !***     --> performed by the proc. diag_targ(1)   ***
    !***  - HDF5_rprof_GC_saving                       ***
    !***     --> performed by the proc. diag_targ(2)   ***
    !***  - HDF5_rprof_part_saving                     ***
    !***     --> performed by the proc. diag_targ(3)   ***
    !***  - HDF5_Phi2D_saving                          ***
    !***     --> performed by the proc. diag_targ(4)   ***
    !***  - HDF5_f2D_saving                            ***
    !***     --> performed by the proc. diag_targ(5)   ***
    !*****************************************************
      
    time_diag = iter_time
      
    call ppbarrier_timer()
    call clck_time(tdeb)   
      
    call clck_time(bclock_diagnostics)
    call clck_time(bclock_physics)
    !*** computation of the moments associated  ***
    !***  to the guiding-center (GC)            *** 
    call compute_GC_moments(geom,coord_sys, &
      init_prof,init_magnet,init_curr, &
      fn,fnp1,Sfmu_eq,dJ0Phidr,dJ0Phidtheta, &
      dJ0Phidphi,dt_diag,iter_time)
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tdiag(1))    
    tdeb = tfin
      
    !*** computation of the moments associated ***
    !***   to the particles                    *** 
    call compute_part_moments(geom,coord_sys, &
      init_prof,init_magnet,init_curr,J0,fnp1,SJ0ftmp, &
      SJ0fmu_eq,poiss%Phi)
      
    call ppbarrier()
    call clck_time(eclock_physics)
    call clck_diff(bclock_physics,eclock_physics, &
      global_time_physics)
      
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tdiag(2))    
    tdeb = tfin
      
    !*** computation of the cross-sections of ***
    !***  the distribution function and the   ***
    !***  potential                           ***
    call clck_time(bclock_crosssection)
    call compute_ftrapped(fnp1,geom)
    call compute_fpassing(fnp1,geom)
      
    call ppbarrier()
    call clck_time(eclock_crosssection)
    call clck_diff(bclock_crosssection,eclock_crosssection, &
      global_time_crosssection)
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tdiag(3))    
    tdeb = tfin
      
    !*** diagnostics saving ***
    call clck_time(bclock_HDF5_diag)
      
    thdf5_diag(1:12) = 0.0
    call clck_time(tdeb)
    !-> Writing of "conservation_laws_<num_diag>.h5" files 
    if (diag_targ(1) .eq. pglobal_id) then
      call HDF5_CL_saving(geom,idiag_num)
    endif
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,thdf5_diag(1))
    tdeb = tfin
      
    !-> Writing of "rprof_GC_<num_diag>.h5" files
    if (diag_targ(2) .eq. pglobal_id) then
      call HDF5_rprof_GC_saving(geom,idiag_num)
    end if
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,thdf5_diag(2))
    tdeb = tfin
      
    !-> Writing of "rprof_part_<num_diag>.h5" files
    if (diag_targ(3) .eq. pglobal_id) then
      call HDF5_rprof_part_saving(geom,idiag_num)
    endif
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,thdf5_diag(3))
    tdeb = tfin
    
    !-> Writing of "Phi2D_<num_diag>.h5" files
    if (diag_targ(4) .eq. pglobal_id) then
      call compute_Phi2DCS(geom,poiss)
      call HDF5_Phi2D_saving(geom,idiag_num)  
    endif
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,thdf5_diag(4))
    tdeb = tfin
    
    !-> Writing of "f2D_<num_diag>.h5" files
    !   ATTENTION: it should be the processor 'cible'
    !              that performs this writing
    if (diag_targ(5) .eq. pglobal_id) then
      call HDF5_f2D_saving(geom,idiag_num)
    endif
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,thdf5_diag(5))
    tdeb = tfin
      
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,thdf5_diag(6))
      
    call ppbarrier()
    call clck_time(eclock_HDF5_diag)
    call clck_diff(bclock_HDF5_diag,eclock_HDF5_diag, &
      global_time_1D2D_diag)
    call clck_time(tfin)
    call clck_ldiff(tdeb,tfin,tdiag(4))    
    
    !*** saving of the 3D results in (r,theta,phi) ***
    !***  in HDF5 format                           ***
#ifndef NOHDF5
    if (idiag_num.ne.-1) then
      !-> Saving of 3D fluid moments and Phi3D(r,theta,phi)
      if ( ((Phi3D_saving) .and. &
        (mod(iter_glob,Phi3D_nbstep).eq.0)) .or. &
        ((FMoments3D_saving) .and. &
        (mod(iter_glob,FMoments3D_nbstep).eq.0)) ) then         
        if (pglobal_id.eq.outputproc) &
          write (uout_res,*) '---> saving of 3D'
        call compute_save_3D(geom,init_prof,init_magnet, &
          init_curr,J0,fnp1,sJ0ftmp,poiss%Phi)
      end if
      
      !*** save f(r,theta,phi,vpar,mu) in HDF5 format ***
      if ((f5D_saving).and.(mod(iter_glob,f5D_nbstep).eq.0)) then
        if (pglobal_id.eq.outputproc) &
          write (uout_res,*) '---> saving of f5D'
        call HDF5_f5D_saving(f5D_nb_diag,geom,fnp1%values)
        f5D_nb_diag = f5D_nb_diag + 1
      end if
    end if
    call ppbarrier()
#endif
      
    call clck_ldiff(bclock_HDF5,eclock_HDF5,tdiag(5))    
    
    call clck_time(eclock_diagnostics)
    call clck_diff(bclock_diagnostics,eclock_diagnostics, &
      global_time_diagnostics)
      
#ifdef TIMER
    write(6,'(I4,A,7F13.6)') pglobal_id, &
      " Temps HDF5 diagnostics ", thdf5_diag(1:6)
    call clck_ldiff(bclock_diagnostics,&
      eclock_diagnostics,tdiag(6))    
    write(6,'(I4,A,5F11.4,A,F11.4)') pglobal_id, &
      " Temps diagnostics ", tdiag(1:5)," tot ",tdiag(6)
#endif
      
!R3 call r3_info_end (r3_info_index_0) !R3
  end subroutine diagnostics
end module diagnostics_module
