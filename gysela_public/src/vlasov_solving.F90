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
! file : vlasov_solving.f90
! date :  16/09/2005
!  Global algorithm for solving the vlasov system
!  which is divided into :
!   - time-splitting algorithm and 
!   - leap-frog or predictor algorithm
!-------------------------------------------------
module vlasov_solving_module
  use globals
  use prec_const
  use geometry_class
  use init_profile_class
  use init_magnetic_class
  use init_current_class
  use fdistribu5d_class
  use poisson_class
!CP!  use diag3D_module
  use diagnostics_module
  use clock_module
  implicit none
      
  !******************************
  contains
  !******************************
  
  !-----------------------------------------------------
  !  Mesh and profile initialisations
  !   - grid definition in (r,theta,phi)
  !   - initialisation of the density and temperature
  !     profiles
  !   - grid definition in vparallel direction
  !  (step only used for the initialisation)
  !-----------------------------------------------------
  subroutine init_mesh_profile(geom,init_prof)
    use globals, only : memory_test, &
      rhomin, rhomax, Ltheta, Lphi, nb_vth0
    type(geometry)    , intent(inout) :: geom
    type(init_profile), intent(inout) :: init_prof
      
    real(RKIND) :: rmin, rmax, thetamin, phimin, vparmin
    real(RKIND) :: Lr, Lvpar
      
    !*** (r,theta,phi) mesh definition ***
    rmin     = a*rhomin
    rmax     = a*rhomax
    Lr       = abs(rmax-rmin)
    thetamin = 0._RKIND
    phimin   = 0._RKIND
    call new_space_coord(geom,rmin,thetamin,phimin, &
      Lr,Ltheta,Lphi,Nr,Ntheta,Nphi)
      
    ! *** Definition of the major radius ***
    R0 = a*aspect_ratio
    
    !*** velocity mesh definition ***
    Lvpar   = 2_RKIND*nb_vth0
    vparmin = -nb_vth0
    call new_velocity_coord(geom,vparmin,mumin,Lvpar,Lmu,Nvpar,Nmu)
      
    !*** profile definition *** 
    call init_all_profile(init_prof,geom)
  end subroutine init_mesh_profile
      
  !--------------------------------------------------------
  !  - initialisation of the distribution function 
  ! at the initial time (with the equilibrium Maxwellian)
  !-------------------------------------------------------- 
  subroutine initialisation_fstart(geom,init_prof,time_init, &
    Sfmu_eq,fstart)
    use globals, only : uout_res, outputproc
    use Pcross_section_module
    type(geometry)    , intent(in)    :: geom
    type(init_profile), intent(in)    :: init_prof
    real(RKIND)       , intent(inout) :: time_init
    real(RKIND), &
      dimension(:,:,:), pointer       :: Sfmu_eq
    type(fdistribu5d) , intent(inout) :: fstart
      
    !*** distribution function initialisation ***
    nbiter_prev = 0
    nb_restart  = 0
    time_init   = 0._RKIND
    iter_glob   = -1
      
    if (pglobal_id.eq.outputproc) then
      write(uout_res,*) '---> Distribution function ', &
        'initialisation : mode NO RESTART'
      write(uout_res,*) ' '
    end if
    call f_initialisation(fstart,geom,init_prof,Sfmu_eq)
  end subroutine initialisation_fstart
      
  !-----------------------------------------------------
  !  - initialisation of the distribution function 
  ! at time 'start+ delta t' (for leap-frog case)
  !----------------------------------------------------- 
  subroutine initialisation_f0(geom,init_prof,fm1,poiss,J0,f0)
    use globals, only : istart, iend, jstart, jend, mu_id
    use advec1D_BSL_module
    use advec2D_BSL_module
    use gyroaverage_class
    type(geometry)    , intent(in)    :: geom
    type(init_profile), intent(in)    :: init_prof
    type(fdistribu5d) , intent(inout) :: fm1
    type(poisson)     , intent(in)    :: poiss
    type(J0operator)  , intent(inout) :: J0
    type(fdistribu5d) , intent(inout) :: f0
      
    integer :: ir, itheta, iphi, ivpar
      
    !*** f0%values = fm1%values ***
    do ivpar = 0,geom%Nvpar
      do iphi = 0,geom%Nphi
        do itheta = jstart,jend
          do ir = istart,iend
            f0%values(ir,itheta,iphi,ivpar) = &
              fm1%values(ir,itheta,iphi,ivpar)    
          end do
        end do
      end do
    end do
  end subroutine initialisation_f0
      
  !-----------------------------------------------------
  !  - initialisation of the distribution function at 
  !      the time t-dt and t
  !  These functions are read into two set 
  !  of restart files as:
  !     . gysela.rst.n<0 or 1>.fnm1.p<num_proc>.h5
  !     . gysela.rst.n<0 or 1>.fn.p<num_proc>.h5
  !  (used in leap-frog case)
  !----------------------------------------------------- 
  subroutine init_fnm1fn_rst(geom,init_prof, &
    time_init,fnm1,fn)
    use globals, only : nb_restart, uout_res, outputproc
    use read_write_module
    
    type(geometry)    , intent(in)    :: geom
    type(init_profile), intent(in)    :: init_prof
    real(RKIND)       , intent(inout) :: time_init
    type(fdistribu5d) , intent(out)   :: fnm1
    type(fdistribu5d) , intent(out)   :: fn
      
    integer           :: ios
    character(LEN=50) :: restart_file
      
    !*** reading of the distribution function fnm1 ***
    write(restart_file(1:),'(A,I1,A,I5.5,A)') &
      "gysela.rst.n",num_restart,".fnm1.p",pglobal_id,".h5"
    if (pglobal_id.eq.outputproc) then
      write(uout_res,'(A,A,A,A)') '---> Distribution function ', &
        'initialisation : mode RESTART with ', &
        restart_file(1:len_trim(restart_file)),'/'
      write(uout_res,*) ' '
    end if
    call read_f5d_rst(restart_file(1:len_trim(restart_file)), &
      time_init,geom,init_prof,"fnm1",fnm1)
      
    !*** reading of the distribution function fn ***
    write(restart_file(1:),'(A,I1,A,I5.5,A)') &
      "gysela.rst.n",num_restart,".fn.p",pglobal_id,".h5"
    if (pglobal_id.eq.outputproc) then
      write(uout_res,'(A,A,A,A)') '---> Distribution function ', &
        'initialisation : mode RESTART with ', &
        restart_file(1:len_trim(restart_file)),'/'
      write(uout_res,*) ' '
    end if  
    call read_f5d_rst(restart_file(1:len_trim(restart_file)), &
      time_init,geom,init_prof,"fn",fn)
      
    nb_restart = nb_restart + 1
    iter_glob  = nbiter_prev
  end subroutine init_fnm1fn_rst
      
  !-----------------------------------------------------
  !  - initialisation of the distribution function at 
  !      the time t
  !  These functions are read into one set 
  !  of restart files as:
  !     . gysela.rst.n<0 or 1>.fn.p<num_proc>.h5
  !  (used in predictor-corrector case)
  !----------------------------------------------------- 
  subroutine init_fn_rst(geom,init_prof, &
    time_init,fn)
    use globals, only : uout_res, outputproc
    use read_write_module
    use globals, only : nb_restart
      
    type(geometry)    , intent(in)    :: geom
    type(init_profile), intent(in)    :: init_prof
    real(RKIND)       , intent(inout) :: time_init
    type(fdistribu5d) , intent(out)   :: fn
      
    integer           :: ios
    character(LEN=50) :: restart_file
      
    !*** reading of the distribution function fn ***
    write(restart_file(1:),'(A,I1,A,I5.5,A)') &
      "gysela.rst.n",num_restart,".fn.p",pglobal_id,".h5"
    if (pglobal_id.eq.outputproc) then
      write(uout_res,'(A,A,A,A)') '---> Distribution function ', &
        'initialisation : mode RESTART with ', &
        restart_file(1:len_trim(restart_file)),'/'
      write(uout_res,*) ' '
    end if  
    call read_f5d_rst(restart_file(1:len_trim(restart_file)), &
      time_init,geom,init_prof,"fn",fn)
    
    nb_restart = nb_restart + 1
    iter_glob  = nbiter_prev
  end subroutine init_fn_rst
      
  !-------------------------------------------------------------
  ! TIME-SPLITTING between t and t+dt for the
  !  SEMI-LAGRANGIAN method IN ALL DIRECTIONS
  ! Computation at the time t, with the shift sequence :
  !     vpar/2,phi/2,(r,theta),phi/2,vpar/2 
  ! to keep a shift symmetry
  !------------------------------------------------------------
  subroutine timesplitting_vlasov_BSL(geom, &
    init_prof,init_magnet,init_curr, &
    dt,fn,Sfmu_eq,activatefilter, &
    maxalprth,maxalpphi,maxalpvpar)
    use globals, only : transp_velocity, advec2Ddt, &
      Rarray1_NrNthetaNphi, Rarray1_NrNthetaNphi
    use MPIutils_module
    use advec1D_BSL_module
    use advec2D_BSL_module
    use filter_module
    type(geometry)       , intent(in)    :: geom
    type(init_profile)   , intent(in)    :: init_prof
    type(init_magnetic)  , intent(in)    :: init_magnet
    type(init_current)   , intent(in)    :: init_curr
    real(RKIND)          , intent(in)    :: dt
    type(fdistribu5d)    , intent(inout) :: fn
    real(RKIND), &
      dimension(:,:,:)   , pointer       :: Sfmu_eq
    logical              , intent(in)    :: activatefilter
    real(RKIND)          , intent(out)   :: maxalprth
    real(RKIND)          , intent(out)   :: maxalpphi
    real(RKIND)          , intent(out)   :: maxalpvpar
    
    integer     :: ivpar, nbsubit
    real(RKIND) :: subdt
      
    if (advec2Ddt .ne. 0._RKIND) then
       nbsubit = MAX(1, nint(dt / advec2Ddt))
       subdt   = dt / nbsubit
    else
       nbsubit = 1
       subdt   = dt
    end if
      
    !*** step 1 : splitting in v parallel ***
    !***  direction (on 0.5*dt)           ***
    call advec1D_BSL_vpar(geom,0.5_RKIND*dt, &
      init_prof,init_magnet,init_curr,fn,maxalpvpar)
    !*** step 2 : splitting in phi direction (on 0.5*dt) ***
    call advec1D_BSL_phi(geom,0.5_RKIND*dt, &
      init_prof,init_magnet,init_curr,fn,maxalpphi)
    !*** step 3 : splitting in (r,theta) (on dt) ***
    if (transpose4D) then
      call transpose_advec2D_xy(geom,subdt, &
        nbsubit,init_prof,init_magnet,init_curr, &
        fn,Sfmu_eq,activatefilter,maxalprth)
    else
      call local_advec2D_xy(geom,dt, &
        init_prof,init_magnet,init_curr,fn,maxalprth)
      if (activatefilter) call hffilter(fn,Sfmu_eq)
    end if
    !*** step 4 : splitting in phi direction (on 0.5*dt) ***
    call advec1D_BSL_phi(geom,0.5_RKIND*dt, &
      init_prof,init_magnet,init_curr,fn,maxalpphi)
    !*** step 5 : splitting in v parallel ***
    !***  direction (on 0.5*dt)           ***
    call advec1D_BSL_vpar(geom,0.5_RKIND*dt, &
      init_prof,init_magnet,init_curr,fn,maxalpvpar)
      
    !-> maximum of the normalised displacments in phi and vpar
    !-> (multiplied by 2, because advection performed on 1/2 dt)
    maxalpvpar = 2._RKIND*maxalpvpar
    maxalpphi  = 2._RKIND*maxalpphi
  end subroutine timesplitting_vlasov_BSL
      
  !****************************************************
  !   LEAP-FROG ALGORITHM FOR VLASOV EQUATION
  !****************************************************
  !---------------------------------------------------------
  ! step 0 : (step only used for the initialisation)
  !  - f initialisation,
  !  - computation of the corresponding Phi,
  !  (used in leap-frog case for the Vlasov system, 
  !    i.e RHS=0)  
  !--------------------------------------------------------- 
  subroutine step0_leapfrog_vlasov(geom,coord_sys, &
    init_prof,init_magnet,init_curr, &
    dt,fn,fnp1,SJ0ftmp,Sfmu_eq,SJ0fmu_eq,poiss,J0)
    use globals          , only : mu_id, &
      Rarray1_NrNthetaNphi, uout_res, outputproc
    use coord_system_class
    use poisson_class
!CP!    use QNsolver_module
    use gyroaverage_class
    use efield_module, only : compute_dPhidx
!R3 #include "r3_info.h" !R3
    type(geometry)     , intent(in)    :: geom
    type(coord_system) , intent(in)    :: coord_sys
    type(init_profile) , intent(in)    :: init_prof
    type(init_magnetic), intent(in)    :: init_magnet
    type(init_current) , intent(in)    :: init_curr
    real(RKIND)        , intent(in)    :: dt
    type(fdistribu5d)  , intent(inout) :: fn
    type(fdistribu5d)  , intent(inout) :: fnp1
    type(fdistribu5d)  , intent(inout) :: SJ0ftmp
    real(RKIND), &
      dimension(:,:,:) , pointer       :: Sfmu_eq
    real(RKIND), &
      dimension(:,:,:) , pointer       :: SJ0fmu_eq
    type(poisson)      , intent(inout) :: poiss
    type(J0operator)   , intent(inout) :: J0
      
    real(RKIND), dimension(:,:,:), pointer :: gyroPhi
      
!R3 call r3_info_begin (r3_info_index_0, 'CPU_step0') !R3
      
    gyroPhi => Rarray1_NrNthetaNphi
    if (pglobal_id.eq.outputproc) then
      write(uout_res,*) ' '
      write(uout_res,*) ' '
      write(uout_res,*) '**************************', &
        '*****************************************'
      write(uout_res,*) '******         ', &
        '  VLASOV SOLVING WITH LEAP-FROG ALGORITHM     ******'
      if (geom%Nmu.eq.0) then
        if (geom%Lmu.eq.0._RKIND) then
          write(uout_res,*) &
            '******                   NO GYROAVERAGE', &
            '                      ******'
        else
          write(uout_res,*) '******           ', &
            'GYROAVERAGE WITH ONE VALUE OF MU            ******'
        end if
      else
        write(uout_res,*) '******                 ', &
          ' WITH GYROAVERAGE                     ******'
      end if
    end if
      
    if (.not.restart) then
      !*** initialisation of f(-1) ***
      call initialisation_fstart(geom,init_prof,init_time, &
        Sfmu_eq,fn)
      !*** computation of Phi(-1) by solving   ***
      !***   the quasi-neutrality equation     ***
      call solve_poisson(poiss,geom,init_prof, &
        init_magnet,init_curr,J0,fn,Sfmu_eq)
      
      !*** computation of the velocity components at time t=-1 ***
      call compute_gyrofunction_3D(J0,geom,mu_id,poiss%Phi,gyroPhi)
      call compute_dPhidx(geom,gyroPhi, &
        dJ0Phidr,dJ0Phidtheta,dJ0Phidphi)
      
      !*** initialisation of f(0) ***
      call initialisation_f0(geom,init_prof,fn,poiss,J0,fnp1)
      
      !*** saving at the initial time t=-1 ***
      call diagnostics(-dt,iter_glob,dt,geom,coord_sys, &
        init_prof,init_magnet,init_curr,fn,fnp1,SJ0ftmp, &
        Sfmu_eq,SJ0fmu_eq,poiss,J0)
      call WriteScreen(-dt,dt,.true.)
      call WriteResults_CL(-dt)
      call WriteResults_rprof(-dt,geom,poiss%Phi)
      iter_glob = iter_glob + 1
    else
      call init_fnm1fn_rst(geom,init_prof,init_time,fn,fnp1)
    end if
      
    !*** computation of Phi(0) by solving    ***
    !***  the quasi-neutrality equation      ***
    call solve_poisson(poiss,geom,init_prof, &
      init_magnet,init_curr,J0,fnp1,Sfmu_eq)
!R3 call r3_info_end (r3_info_index_0) !R3
  end subroutine step0_leapfrog_vlasov
      
  !------------------------------------------------------------
  ! Solve the Vlasov equation for the distribution function f
  !  between t-dt and t+dt.
  ! Computation at the time t, with the shift sequence :
  !   vpar,phi,2*(r,theta),phi,vpar to keep a shift symmetry
  ! Use a LEAP-FROG algorithm
  !---------------------------------------------------------  
  subroutine global_leapfrog_vlasov(geom,init_prof,init_magnet, &
    init_curr,dt,fnm1,fn,fnp1,Sfmu_eq,J0,poiss)
    use globals, only : mu_id, Rarray1_NrNthetaNphi, &
      uout_res, outputproc
    use poisson_class
!CP!    use QNsolver_module
    use gyroaverage_class
    use efield_module         , only : compute_dPhidx
    type(geometry)     , intent(in)    :: geom
    type(init_profile) , intent(in)    :: init_prof
    type(init_magnetic), intent(in)    :: init_magnet
    type(init_current) , intent(in)    :: init_curr
    real(RKIND)        , intent(inout) :: dt
    type(fdistribu5d)  , pointer       :: fnm1
    type(fdistribu5d)  , pointer       :: fn
    type(fdistribu5d)  , pointer       :: fnp1
    real(RKIND), &
      dimension(:,:,:) , pointer       :: Sfmu_eq
    type(J0operator)   , intent(inout) :: J0
    type(poisson)      , intent(inout) :: poiss      
      
    real(RKIND), &
      dimension(:,:,:), pointer :: gyroPhi
    integer                     :: ir, itheta, iphi, ivpar
    integer(TIMEPREC)           :: tdeb, tfin
    type(fdistribu5d) , pointer :: swap
    !-> for computation of the maximum normalised displacements
    !->  (they are computed but not used in the leap-frog
    !->   algorithm because no adaptative time step is possible
    !->   in leap-frog global scheme)
    real(RKIND) :: maxalprth, maxalpphi, maxalpvpar
      
!R3 #include "r3_info.h" !R3
      
    !*** leap-frog iteration ***
    ! => fn%values = fnm1%values (if average)
    ! (prepared for computing 
    !   fn%values = 0.5*(fnm1%values+fnp1%values) later)
!R3 call r3_info_begin(r3_info_index_0, 'CPU_copy_fdistrib') !R3
    call clck_time(tdeb)
    gyroPhi => Rarray1_NrNthetaNphi
      
    ! => fnm1%values = fn%values (for the next step)
    swap => fnm1
    fnm1 => fn
    fn   => swap
    ! => fn%values = fnp1%values
    swap => fn
    fn   => fnp1
    fnp1 => swap
      
    if (mod(iter_glob,LF_average).eq.0) then
      do ivpar = 0,geom%Nvpar
        do iphi = 0,geom%Nphi
          do itheta = fn%jstart_buf,fn%jend_buf
            do ir = fn%istart_buf,fn%iend_buf
              fn%values(ir,itheta,iphi,ivpar) = &
                fnm1%values(ir,itheta,iphi,ivpar)
            end do
          end do
        end do
      end do
    end if
!R3 call r3_info_end(r3_info_index_0) !R3
      
    !*** fnp1%values = fnm1%values ***
    call clck_time(tdeb)
    swap => fnp1
    fnp1 => fnm1
    fnm1 => swap
    call clck_time(tfin)
    call clck_diff(tdeb,tfin,global_time_fncopy)
      
    !*** time-splitting iteration                            ***
    !***  . f(tn+1) computed with f(tn-1) on 2dt using E(tn) ***
    call compute_gyrofunction_3D(J0,geom,mu_id,poiss%Phi,gyroPhi)
    call compute_dPhidx(geom,gyroPhi, &
      dJ0Phidr,dJ0Phidtheta,dJ0Phidphi)
    call timesplitting_vlasov_BSL(geom, &
      init_prof,init_magnet,init_curr, &
      2._RKIND*dt,fnp1,Sfmu_eq,.true., &
      maxalprth,maxalpphi,maxalpvpar)
      
    !*** computation of Phi(tn+1) by solving    ***
    !***  the quasi-neutrality equation         ***
    call solve_poisson(poiss,geom,init_prof, &
      init_magnet,init_curr,J0,fnp1,Sfmu_eq)
      
    !*** leap-frog iteration ***
    ! => fn%values =  0.5*(fnm1%values+fnp1%values) (if average)
    call clck_time(tdeb)
    if (mod(iter_glob,LF_average).eq.0) then
      if (pglobal_id.eq.outputproc) then
        write(uout_res,'(A17,I5,A30)') &
          '---> iteration = ',iter_glob, &
          ' Average between fn and fnp1'
      end if
!R3 call r3_info_begin(r3_info_index_0, 'CPU_copy_fdistrib') !R3
      do ivpar = 0,geom%Nvpar
        do iphi = 0,geom%Nphi-1
          do itheta = fn%jstart_buf,fn%jend_buf
            do ir = fn%istart_buf,fn%iend_buf
              fn%values(ir,itheta,iphi,ivpar) = &
                0.5_RKIND*(fn%values(ir,itheta,iphi,ivpar) + &    
                fnp1%values(ir,itheta,iphi,ivpar))
            end do
          end do
        end do
        !=> periodicity in phi direction
        do itheta = fn%jstart_buf,fn%jend_buf
          do ir = fn%istart_buf,fn%iend_buf
            fn%values(ir,itheta,geom%Nphi,ivpar) = &
              fn%values(ir,itheta,0,ivpar)
          end do
        end do
      end do
!R3 call r3_info_end(r3_info_index_0) !R3
    end if
      
    call clck_time(tfin)
    call clck_diff(tdeb,tfin,global_time_fncopy)
  end subroutine global_leapfrog_vlasov
      
  !****************************************************
  !   PREDICTOR-CORRECTOR ALGORITHM FOR THE BSL METHOD
  !****************************************************
  !---------------------------------------------------------
  ! step 0 : (step only used for the initialisation)
  !  - f initialisation,
  !  - computation of the corresponding Phi,
  !  (used in predictor-corrector case)  
  !--------------------------------------------------------- 
  subroutine step0_predcorr_vlasov(geom,coord_sys, &
    init_prof,init_magnet,init_curr, &
    dt,fnp1,SJ0ftmp,Sfmu_eq,SJ0fmu_eq,poiss,J0)
    use globals          , only : mu_id, &
      Rarray1_NrNthetaNphi, uout_res, outputproc
    use coord_system_class
    use poisson_class
!CP!    use QNsolver_module
    use gyroaverage_class
    use efield_module         , only : compute_dPhidx
!R3 #include "r3_info.h" !R3
    type(geometry)     , intent(in)    :: geom
    type(coord_system) , intent(in)    :: coord_sys
    type(init_profile) , intent(in)    :: init_prof
    type(init_magnetic), intent(in)    :: init_magnet
    type(init_current) , intent(in)    :: init_curr
    real(RKIND)        , intent(in)    :: dt
    type(fdistribu5d)  , intent(inout) :: fnp1
    type(fdistribu5d)  , intent(inout) :: SJ0ftmp
    real(RKIND), &
      dimension(:,:,:) , pointer       :: Sfmu_eq
    real(RKIND), &
      dimension(:,:,:) , pointer       :: SJ0fmu_eq
    type(poisson)      , intent(inout) :: poiss
    type(J0operator)   , intent(inout) :: J0
      
    real(RKIND), dimension(:,:,:), pointer :: gyroPhi
      
!R3 call r3_info_begin (r3_info_index_0, 'CPU_step0') !R3
 
    gyroPhi => Rarray1_NrNthetaNphi
    if (pglobal_id.eq.outputproc) then
      write(uout_res,*) ' '
      write(uout_res,*) ' '
      write(uout_res,*) '**************************', &
        '*****************************************'
      write(uout_res,*) '******    ', &
        '  SOLVING WITH PREDICTOR-CORRECTOR ALGORITHM       ******'
      if (geom%Nmu.eq.0) then
        if (geom%Lmu.eq.0._RKIND) then
          write(uout_res,*) &
            '******                   NO GYROAVERAGE', &
            '                      ******'
        else
          write(uout_res,*) '******           ', &
            'GYROAVERAGE WITH ONE VALUE OF MU            ******'
        end if
      else
        write(uout_res,*) '******                 ', &
          ' WITH GYROAVERAGE                     ******'
      end if
      write(uout_res,*) '******                  ', &
        '                                     ******'
      write(uout_res,*) &
        '******       WITH PRESCRIBED TEMPERATURE', &
        ' AT THE EDGE         ******'
      write(uout_res,*) '********************************', &
        '***********************************'
      write(uout_res,*) ' '
    end if
      
    if (.not.restart) then
      !*** initialisation of f(-1) ***
      call initialisation_fstart(geom,init_prof,init_time, &
        Sfmu_eq,fnp1)
      !*** computation of Phi(-1) by solving   ***
      !***   the quasi-neutrality equation     ***
      call solve_poisson(poiss,geom,init_prof, &
        init_magnet,init_curr,J0,fnp1,Sfmu_eq)
      
      !*** computation of the velocity components at time t=-1 ***
      call compute_gyrofunction_3D(J0,geom,mu_id,poiss%Phi,gyroPhi)
      call compute_dPhidx(geom,gyroPhi, &
        dJ0Phidr,dJ0Phidtheta,dJ0Phidphi)
      
      !*** saving at the initial time t=-1 ***
      call diagnostics(-dt,iter_glob,dt,geom,coord_sys, &
        init_prof,init_magnet,init_curr,fnp1,fnp1,SJ0ftmp, &
        Sfmu_eq,SJ0fmu_eq,poiss,J0)
      call WriteScreen(-dt,dt,.true.)
      call WriteResults_CL(-dt)
      call WriteResults_rprof(-dt,geom,poiss%Phi)
      
      !*** initialisation of f(0) ***
      iter_glob = iter_glob + 1
    else
      call init_fn_rst(geom,init_prof,init_time,fnp1)
      
      !***  the quasi-neutrality equation      ***
      call solve_poisson(poiss,geom,init_prof, &
        init_magnet,init_curr,J0,fnp1,Sfmu_eq)
    end if
!R3 call r3_info_end (r3_info_index_0) !R3
  end subroutine step0_predcorr_vlasov
      
  !------------------------------------------------------------
  ! Solve the Vlasov equation for the distribution function f
  !  between t and t+dt by using a predictor-corrector
  !  scheme
  !---------------------------------------------------------  
  subroutine global_predcorr_vlasov(geom,init_prof, &
    init_magnet,init_curr,dt,iter_time,time_save, &
    fn,fnp1,Sfmu_eq,J0,poiss)
    use globals
    use gyroaverage_class
    use poisson_class
!CP!    use QNsolver_module
    use efield_module, only : compute_dPhidx
    include "mpiperso.h"
    type(geometry)     , intent(in)    :: geom
    type(init_profile) , intent(in)    :: init_prof
    type(init_magnetic), intent(in)    :: init_magnet
    type(init_current) , intent(in)    :: init_curr
    real(RKIND)        , intent(in)    :: iter_time
    real(RKIND)        , intent(in)    :: time_save
    real(RKIND)        , intent(out)   :: dt
    type(fdistribu5d)  , pointer       :: fn
    type(fdistribu5d)  , pointer       :: fnp1
    real(RKIND), &
      dimension(:,:,:) , pointer       :: Sfmu_eq
    type(J0operator)   , intent(inout) :: J0
    type(poisson)      , intent(inout) :: poiss
      
    integer            :: ir, itheta, iphi, ivpar, ipred, ierr
    integer(TIMEPREC)  :: tdeb, tfin
    real(RKIND)        :: predratio, effectivedt, coeffadapt
    integer            :: itadapt, nbadapt
    integer, parameter :: MAXPRED = 2
    logical            :: activatefilter
    real(RKIND), &
      dimension(:,:,:), pointer :: gyroPhi
    type(fdistribu5d) , pointer :: swap
    !-> for time step adaptation
    real(RKIND) :: maxalprth, maxalpphi, maxalpvpar, maxalp
    real(RKIND) :: dt_muid_reduce, dt_muid_diag
    integer     :: dt_direction_reduce
      
    gyroPhi => Rarray1_NrNthetaNphi
      
    if ((reduceddt .ne. 0._RKIND) .and. &
      (iter_time .ge. reducedbegin) .and. &
      (iter_time .le. reducedend)) then
      nbadapt    = MAX(1, nint(dt / reduceddt))
      coeffadapt = 1._RKIND / nbadapt
    else
      nbadapt    = 1
      coeffadapt = 1._RKIND
    end if
    
    do itadapt = 1,nbadapt
      swap => fn
      fn   => fnp1
      fnp1 => swap
      do ipred = MAXPRED,1,-1
        activatefilter = (itadapt .eq. nbadapt) .and. &
          (ipred .eq. 1)
        predratio   = (2._RKIND ** (ipred-1))
        predratio   = coeffadapt / predratio
        effectivedt = predratio * dt
        !*** time-splitting iteration to compute ***
        !***  f(tn+1/(2^ipred))                  ***
        call compute_gyrofunction_3D(J0,geom,mu_id, &
          poiss%Phi,gyroPhi)
        call compute_dPhidx(geom,gyroPhi, &
          dJ0Phidr,dJ0Phidtheta,dJ0Phidphi)
        do ivpar = 0,geom%Nvpar
          do iphi = 0,geom%Nphi
            do itheta = fn%jstart_buf,fn%jend_buf
              do ir = fn%istart_buf,fn%iend_buf
                fnp1%values(ir,itheta,iphi,ivpar) = &
                  fn%values(ir,itheta,iphi,ivpar)
              end do
            end do
          end do
        end do
        call timesplitting_vlasov_BSL(geom,init_prof,init_magnet, &
          init_curr,effectivedt,fnp1,Sfmu_eq,activatefilter, &
          maxalprth,maxalpphi,maxalpvpar)
      
        !*** computation of Phi(tn+1) by solving ***
        !***  the quasi-neutrality equation      ***
        call solve_poisson(poiss,geom,init_prof, &
          init_magnet,init_curr,J0,fnp1,Sfmu_eq)
      end do
    end do
      
    !*** computation of the adaptative time step ***
    !-> local maximum displacement 
    maxalp = max(maxalpphi/8._RKIND,maxalprth/0.5_RKIND)
    maxalp = max(maxalp,maxalpvpar/8._RKIND)
    !-> give the direction of the biggest normalized displacement
    if (maxalp == maxalprth/0.5_RKIND)  &
      dt_direction_reduce = 1
    if (maxalp == maxalpphi/8._RKIND)  &
      dt_direction_reduce = 3
    if (maxalp == maxalpvpar/8._RKIND)  &
      dt_direction_reduce = 4
    !-> local adaptative time step
    dt = min(dt_save,dt/maxalp)
    if (iter_time+dt>time_save.and.time_save.ne.iter_time) &
      dt = time_save - iter_time
    !-> global adaptative time step
    dt_muid_reduce = dt
    call MPI_ALLREDUCE(dt_muid_reduce,dt_muid_diag,1, &
      MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ierr)
    dt = dt_muid_diag
    !-> print the direction responsible of 
    !->  the time step changing 
    !->  (most restrictive CFL condition)
    if (dt.eq.dt_muid_reduce) then
      write(6,*) "==> dt_direction_reduce = ",dt_direction_reduce
    end if
  end subroutine global_predcorr_vlasov
end module vlasov_solving_module
