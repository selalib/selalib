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
      
!---------------------------------------------------------------
! file : gysela_pub.f90
! date : 24/11/2010
! main program for public version
!                      GYROKINETIC CODE
!  - f(r,theta,phi,vpar,mu) : 5D ionic distribution function
!  - Phi(r,theta,phi)       : 3D electric potential
!  Solving of a 5D gyroaverage vlasov equation coupled to a 3D 
!   quasi-neutrality equation
!---------------------------------------------------------------
program gysela
  use globals  
  use clock_module
  use OMPutils_module
  use MPIutils_module
  use geometry_class
  use gys_alloc_module
  use init_profile_class
  use coord_system_class
  use init_magnetic_class
  use init_current_class
  use read_write_module
  use fdistribu5d_class
  use interpolation_module
  use integration_module
  use gyroaverage_class
!baoter
  use poisson_class
!CP!  use QNsolver_module
!eaoter
  use advec2D_BSL_module
  use vlasov_solving_module
  use physics_module
  use Pcross_section_module
  use output_saving_module
  use efield_module
  use diagnostics_module
#ifndef NOHDF5
  use f5D_saving_module
#endif
      
  implicit none
      
  include "mpiperso.h"
  !*** variable definition ***
  !-> mesh
  type(geometry)              :: geom
  !-> initial profiles
  type(init_profile)          :: init_prof
  !-> definition of the coordinate system
  type(coord_system)          :: coord_sys
  !-> initial plasma current
  type(init_current)          :: init_curr
  !-> initial magnetic configuration   
  type(init_magnetic)         :: init_magnet
  !-> distribution function f(tn-1)
  type(fdistribu5d), pointer :: pfnm1
  !-> distribution function f(tn)
  type(fdistribu5d), pointer :: pfn
  !-> distribution function f(tn+1)
  type(fdistribu5d), pointer :: pfnp1
  !-> distribution function 1
  type(fdistribu5d), target  :: fna
  !-> distribution function 2
  type(fdistribu5d), target  :: fnb
  !-> distribution function 3
  type(fdistribu5d), target  :: fnc
  !-> J0.f
  type(fdistribu5d), pointer :: J0ftmp
  !-> gyroaveraged distribution function
  type(fdistribu5d)          :: gyrofn
  type(poisson)              :: poiss
  !-> equilibrium distribution functions 
  real(RKIND), dimension(:,:,:), pointer :: fmu_equil
  real(RKIND), dimension(:,:,:), pointer :: J0fmu_equil
  
  real(RKIND) :: iter_time, time_save
      
  ! gyroaverage operator
  type(J0operator) :: J0
      
  ! used for output
  integer     :: nbsave, Nbt
  real(RKIND) :: max_memory
  integer     :: ierr_ce, ios, ierr
      
  ! used for interrupt the run
  logical :: fexist, run
  integer :: status
      
  !*** initialisation for the parallelisation ***
  call ppinit()
      
!R3 file_r3_open = .false. !R3
!R3 call r3_info_init () !R3
      
  max_allocate = 0
  nb_allocate  = 0
    
  !*** restart files exist ***
  call init_restart()
      
  !*** open output file results
  if (pglobal_id.eq.outputproc) &
    open(unit=uout_res,file=file_name_res,status='UNKNOWN', &
    position='APPEND',form='FORMATTED')
      
  !*** input file reading ***
  call read_input()
      
  !*** creation of the file for the memory size saving ***
  call new_allocate()
      
  !*** initialisation of the communicators ***
  call ppinit_comm()
      
  !*** allocation of the OpenMP arrays ***
  call init_OMP()
      
  !*** input parameter saving ***
  call Write_parameters()
  
  call clck_time(bclock_init)
      
  if (pglobal_id.eq.outputproc) &
    close(uout_res)
      
  !*************************************************
  ! ARRAY ALLOCATION AND INITIALIZATION
  !*************************************************
  iter_run        = 0
      
  if (pglobal_id.eq.outputproc) then
    !*** open output file results
    open(unit=uout_res,file=file_name_res,status='OLD',&
      position='APPEND', form='FORMATTED')
  end if
      
  !*** power for the integral computation ***
  ipow = 1
  if ((Nmu.eq.0).and.(Lmu.eq.0)) ipow = 0
      
  !*** mesh and profile initializations ***
  call init_mesh_profile(geom,init_prof)
  
  !*** initialization of the coordinate system ***
  call init_coordinate_system(coord_sys,geom)
      
  !*** initialization of magnetic configuration ***
  call init_magnetic_config(init_magnet,geom, &
    coord_sys,init_prof,init_curr)
      
  !*** initialization of the plasma current ***
  call init_plasma_current(init_curr,geom,coord_sys, &
    init_prof,init_magnet%zeta)
  
  call init_precompute_Bstar(init_magnet,init_curr, &
    geom,coord_sys,init_prof)
      
  !*** allocation of the distribution functions ***
  call new_f5D(fna,geom)
  call new_f5D(fnb,geom)
  if (leapfrog) &
    call new_f5D(fnc,geom)
  
  !*** temporary array allocations ***
  call allocate_temporary_arrays(geom) 
  
  !*** J0(Phi) allocation ***
  call new_J0operator(J0,geom)
      
  !*** Phi and Efield allocation ***
  call new_poisson(poiss,geom,init_prof)
      
  !*** allocation for physic value computations ***
  call new_physics(geom)
  call new_efield(geom)
  call new_integration(geom)
  call new_Pcross_section()
      
  nbsave = int((nbiter)/diag_nbstep)-1 + &
    min(1,mod(nbiter,diag_nbstep))
  if (.not.restart) nbsave = nbsave+1
      
  !*** allocation and computation of the equilibrium part ***
  call allocate_equilibrium_part(geom,fmu_equil,J0fmu_equil, &
    M2_equil,nGieq_rtheta,neq_r,dneq_dr,Tieq_r)
  if (.not.memory_test) then
    call compute_equilibrium_part(geom,J0,init_prof, &
      init_magnet,init_curr,fmu_equil,J0fmu_equil, &
      M2_equil,nGieq_rtheta,neq_r,dneq_dr,Tieq_r,nbions_eq)
  end if
      
  !*** allocation for the HDF5 file saving ***
#ifndef NOHDF5
  if (f5D_saving) call new_f5D_saving()
#endif
      
  !*** initialization of the time variables ***
  call initialization_time
      
  !*** Saving of the initial state ***
  if (.not.memory_test) then
    if ((.not.restart).and.(pglobal_id.eq.0)) then
      call HDF5_coordsys_saving(geom,coord_sys)
    end if
    if (pglobal_id.eq.0) then
#ifndef NOHDF5
        call HDF5_init_saving(geom,init_prof,init_magnet, &
          init_curr)
#else
        call ASCII_init1D_saving(geom,init_prof,init_magnet)
#endif
    end if
  end if
      
  !*** write max memory allocation ***
  max_memory = max_allocate
  call print_memory_size(max_memory)
      
  if (pglobal_id.eq.outputproc) &
    close(uout_res)
      
  if (memory_test) then
    call ppexit()
    status = -1
    if (pglobal_id.eq.0) print *,' STOP ! memory_test=.true. '
    call exit(status)
  end if
      
  !*************************************************
  ! FIRST STEP
  !*************************************************
      
!bmodif (just done temporary for solving BULL bug)
#ifdef _OPENMP
  call omp_set_dynamic(.false.)
  call omp_set_num_threads(Nbthread)
!$OMP PARALLEL default(shared)
!$OMP BARRIER
!$OMP MASTER
  Nbt = omp_get_num_threads()
!$OMP END MASTER
!$OMP BARRIER
!$OMP END PARALLEL    
#else 
  Nbthread = 1
#endif
!emodif
 
  if (pglobal_id.eq.outputproc) then
    !*** open output file results
    open(unit=uout_res,file=file_name_res,status='OLD',&
      position='APPEND', form='FORMATTED')
  end if
      
  if (leapfrog) then
    pfnm1  => fna
    pfn    => fnb
    pfnp1  => fnc
    J0ftmp => pfnm1
    ! -> initialization of fm1 and computation of f0
    call step0_leapfrog_vlasov(geom,coord_sys, &
      init_prof,init_magnet,init_curr, &
      deltat,pfn,pfnp1,J0ftmp,fmu_equil,J0fmu_equil,poiss,J0)
  else
    pfn    => fna
    pfnp1  => fnb
    J0ftmp => pfn
    ! -> initialization of fm1 and computation of f0
    call step0_predcorr_vlasov(geom,coord_sys, &
      init_prof,init_magnet,init_curr, &
      deltat,pfnp1,J0ftmp,fmu_equil,J0fmu_equil,poiss,J0)
  end if
      
  iter_time = init_time
  dt_save   = deltat*diag_nbstep
  time_save = init_time+dt_save
      
  if (.not. restart) then
    !*** saving at the time t=0  ***
    if (leapfrog) then
      J0ftmp => pfnm1
    else
      J0ftmp => pfn
    end if
    call diagnostics(iter_time,0,deltat,geom,coord_sys, &
      init_prof,init_magnet,init_curr,pfn,pfnp1,J0ftmp, &
      fmu_equil,J0fmu_equil,poiss,J0)
    call WriteScreen(iter_time,deltat,.true.)
    call WriteResults_CL(iter_time)
    call WriteResults_rprof(iter_time,geom,poiss%Phi)
  end if
      
  if (pglobal_id.eq.outputproc) &
    close(uout_res)
      
  !*************************************************
  ! GLOBAL SOLVING BEGINNING
  !*************************************************
!R3 call r3_info_print (-2, -2, 'INITIALIZATION') !R3
      
  call clck_time(eclock_init)
  iter_run  = 1
  run       = .true.
  status    = 0
      
  do while (run)
    if (pglobal_id.eq.outputproc) then
      !*** open file results
      open(unit=uout_res,file=file_name_res,status='OLD',&
        position='APPEND', form='FORMATTED')
    end if
      
    iter_glob = iter_glob + 1
    iter_time = iter_time + deltat
    if (leapfrog) then
      call global_leapfrog_vlasov(geom,init_prof,init_magnet, &
        init_curr,deltat,pfnm1,pfn,pfnp1,fmu_equil,J0,poiss)
    else
      call global_predcorr_vlasov(geom, &
        init_prof,init_magnet,init_curr, &
        deltat,iter_time,time_save,pfn,pfnp1,fmu_equil,J0,poiss)
    end if
    if ( (iter_time.ge.time_save) &
      .or.(iter_run.eq.nbiter) ) then
      nb_diag   = nb_diag + 1
      time_save = iter_time + dt_save
      !*** saving at the time t=tn+1 ***
      if (leapfrog) then
        J0ftmp => pfnm1
      else
        J0ftmp => pfn
      end if
      call diagnostics(iter_time,nb_diag,deltat,geom,coord_sys, &
        init_prof,init_magnet,init_curr,pfn,pfnp1,J0ftmp, &
        fmu_equil,J0fmu_equil,poiss,J0)
      call WriteScreen(iter_time,deltat,.true.)
      call WriteResults_CL(iter_time)
      call WriteResults_rprof(iter_time,geom,poiss%Phi)
      inquire(file='gysela.stop', exist=fexist)
      if (fexist) then
        if (pglobal_id.eq.0) then
          write(6,*) ' '
          write(6,'(A)') &
            ' ===> The signal file "gysela.stop" exist '
          write(6,'(A,I6,A,E12.5)') &
            ' ===> So! stop the run at iter_glob = ', &
            iter_glob,' iter_time = ',iter_time
        end if
        run    = .false.
        call MPI_BCAST(run,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        status = 1
      end if
    else
      call WriteScreen(iter_time,deltat,.false.)
    endif
!R3    if ((iter_run .eq. 1).and.(.not.restart)) then     !R3
!R3      call r3_info_print (-3, -2, 'ITERATION    1')    !R3
!R3    else                                               !R3
!R3      call r3_info_print (iter_glob, -2, 'ITERATION')  !R3
!R3    end if                                             !R3
      
    iter_run = iter_run+1
    if (iter_run.gt.nbiter) run = .false.
    if (pglobal_id.eq.outputproc) &
      close(uout_res)
  enddo
      
  !*************************************************
  ! GLOBAL SOLVING END
  !************************************************* 
  if (pglobal_id.eq.outputproc) then
    open(unit=uout_res,file=file_name_res,status='OLD',&
      position='APPEND', form='FORMATTED')
  end if
      
  nbiter_prev = iter_glob  
  if (pglobal_id.eq.0) num_restart = mod(num_restart+1,2)
  call MPI_BCAST(num_restart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
  if (leapfrog) then
    call write_f5d_rst(iter_time,"fnm1",pfn)
    call write_f5d_rst(iter_time,"fn",pfnp1)
  else
    call write_f5d_rst(iter_time,"fn",pfnp1)
  end if
      
  !*** Global barrier in order to guaranty that ***
  !*** restart files are well written           *** 
  call ppbarrier()
      
  !*** writting of the files 'file_list.out' and       ***
  !*** "num_restart.out' (used for the restart script) ***
  call end_restart()
      
  !*** write infos on the CPU time ***
  call write_time()
!R3 call r3_info_summary () !R3
      
  !*** write the number of particles out of domain ***
  if (plocal_id.eq.0) then
    call print_particles_outofdomain()
  end if
      
  !*** array deallocations ***
  !-> deallocate arrays used for OpenMP parallelization
  call del_OMP()
  call del_geometry(geom)
  call del_init_profile(init_prof)
  call del_coord_system(coord_sys)
  call del_magnet_config(init_magnet)
  call del_init_current(init_curr)
  call deallocate_equilibrium_part(fmu_equil,J0fmu_equil, &
    M2_equil,nGieq_rtheta,neq_r,dneq_dr,Tieq_r)
  !-> deallocate distribution functions
  if (leapfrog) then
    call del_f5D(fna)
    call del_f5D(fnb)
    call del_f5D(fnc)
  else
    call del_f5D(fna)
    call del_f5D(fnb)
  end if
      
  !-> deallocate poisson + gyroaverage and collision operators
  call del_poisson(poiss)
  call del_J0operator(J0)
  !-> deallocate physical quantities 
  call delete_integration()
  call delete_physics()
  call delete_efield()
  call delete_Pcross_section()
  !-> deallocate 5D arrays
#ifndef NOHDF5
  if (f5D_saving) call delete_f5D_saving()
#endif
  !-> deallocate temporary arrays
  call deallocate_temporary_arrays()
  !-> deallocate parallel buffers 
  call ppdeallocate()
      
  if (pglobal_id.eq.outputproc) then
    write(uout_res,*) 'memory size at the end of the run = ', &
    nb_allocate, 'Bytes'
    close(uout_res)
  end if
      
  !*** end of the parallelisation ***
  call ppexit()
  if (.not.run) call exit(status)
end program gysela
