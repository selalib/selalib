!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
!> @author
!> Michel Mehrenberger (mehrenbe@math.unistra.fr)
!> Edwin Chacon Golcher
!> @brief 
!> Simulation class to solve slab drift kinetic equation in polar coordinates
!> (3d space (x1=r,x2=theta,x3=z) 1d velocity (x4=v))
!> translation of slv2d/src/vp4d_dk.F90 program in simulation class
!> intended to be close to sll_simulation_4d_DK_hybrid_module
!> but specific use of polar coordinates
!> should merge with sll_simulation_4d_DK_hybrid_module in future
!> @details
!> Example of use in test program: see unit_test_4d_dk_polar.F90 file
!> 
!> \code
!>
!>  use sll_simulation_4d_drift_kinetic_polar
!>  type(sll_simulation_4d_vp_polar)    :: simulation
!>  call simulation%init_from_file(trim(filename))
!>  call simulation%run()
!>  call delete(simulation)
!> \endcode


module sll_simulation_4d_drift_kinetic_polar_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
  use sll_collective
  use sll_remapper
  use sll_constants
  use sll_cubic_spline_interpolator_1d
  use sll_test_4d_initializer
  use sll_poisson_2d_periodic_cartesian_par
  use sll_cubic_spline_interpolator_1d
  use sll_simulation_base
  use sll_fdistribu4D_DK
  use sll_logical_meshes
  use polar_operators
  use polar_advection
  use sll_reduction_module

  implicit none

!! choice of QNS solver
!! should be else where
  sll_int32, parameter :: SLL_NO_QUASI_NEUTRAL = 0 
  sll_int32, parameter :: SLL_QUASI_NEUTRAL_WITH_ZONAL_FLOW = 1 
  sll_int32, parameter :: SLL_QUASI_NEUTRAL_WITHOUT_ZONAL_FLOW = 2 

!! choice of time scheme solver
!! should be else where
  sll_int32, parameter :: SLL_TIME_LOOP_EULER = 0 
  sll_int32, parameter :: SLL_TIME_LOOP_PREDICTOR_CORRECTOR = 1 

!! choice of characteristics scheme
!! should be else where
  sll_int32, parameter :: SLL_CARAC_EULER = 0 
  sll_int32, parameter :: SLL_CARAC_VERLET = 1 



  
  type, extends(sll_simulation_base_class) :: &
    sll_simulation_4d_drift_kinetic_polar

     ! Parallel environment parameters
     sll_int32  :: world_size
     sll_int32  :: my_rank
     sll_int32  :: power2 ! 2^power2 = number of processes available
     ! Processor mesh sizes
     sll_int32  :: nproc_x1
     sll_int32  :: nproc_x2
     sll_int32  :: nproc_x3
     sll_int32  :: nproc_x4 
     ! Mesh parameters
     type(sll_logical_mesh_1d), pointer :: m_x1
     type(sll_logical_mesh_1d), pointer :: m_x2
     type(sll_logical_mesh_1d), pointer :: m_x3
     type(sll_logical_mesh_1d), pointer :: m_x4
     !sll_real64 :: r_min
     !sll_real64 :: r_max
     !sll_real64 :: phi_min
     !sll_real64 :: phi_max
     !sll_real64 :: vpar_min
     !sll_real64 :: vpar_max
     ! Physics/numerical parameters
     sll_real64 :: dt
     sll_int32  :: num_iterations
     sll_int32  :: time_case
     sll_int32  :: carac_case
     !sll_int32  :: spline_degree_eta1, spline_degree_eta2
     !sll_int32  :: spline_degree_eta3, spline_degree_eta4
     !--> Equilibrium
     sll_real64 :: tau0      !-> tau0 = Ti(rpeak)/Te(rpeak)
     sll_real64 :: rho_peak    
     sll_real64 :: kappan   
     sll_real64 :: deltarn  
     sll_real64 :: kappaTi  
     sll_real64 :: deltarTi 
     sll_real64 :: kappaTe  
     sll_real64 :: deltarTe
     sll_int32  :: QN_case     
     sll_real64 :: n0_at_rpeak     
     !--> Pertubation
     sll_int32  :: perturb_choice
     sll_int32  :: mmode
     sll_int32  :: nmode
     sll_real64 :: eps_perturb   

     !--> 4D logical mesh (r,theta,phi,vpar)
     !type(sll_logical_mesh_4d), pointer :: logical_mesh4d


     !--> Density and temperature profiles
     sll_real64, dimension(:)  , pointer :: n0_r
     sll_real64, dimension(:)  , pointer :: Ti_r
     sll_real64, dimension(:)  , pointer :: Te_r

     !--> Equilibrium distribution function
     sll_real64, dimension(:,:), pointer :: feq_x1x4


     !--> 4D distribution function 
     !----> sequential in (x1,x2,x4) and parallel in (x3)
     type(layout_4D), pointer :: layout4d_x1x2x4
     sll_real64, dimension(:,:,:,:), pointer :: f4d_x1x2x4 
     !----> parallel in (x3) and sequential in (x1,x2,x4) 
     type(layout_4D), pointer :: layout4d_x3
     sll_real64, dimension(:,:,:,:), pointer :: f4d_x3
     !----> definition of remap
     type(remap_plan_4D_real64), pointer ::remap_plan_x1x2x4_x3
     type(remap_plan_4D_real64), pointer ::remap_plan_x3_x1x2x4
     

     !--> 3D charge density and 3D electric potential
     !----> sequential in (x1,x2)
     type(layout_3D), pointer :: layout3d_x1x2
     sll_real64, dimension(:,:,:), pointer :: rho3d_x1x2 
     sll_real64, dimension(:,:,:), pointer :: phi3d_x1x2 
     sll_real64, dimension(:,:,:), pointer :: dx1_phi3d_x1x2 
     sll_real64, dimension(:,:,:), pointer :: dx2_phi3d_x1x2 
     !----> sequential in x3
     type(layout_3D), pointer :: layout3d_x3
     sll_real64, dimension(:,:,:), pointer :: rho3d_x3
     sll_real64, dimension(:,:,:), pointer :: phi3d_x3
     sll_real64, dimension(:,:,:), pointer :: dx3_phi3d_x3
     !----> definition of remap
     type(remap_plan_3D_real64), pointer ::remap_plan_x1x2_x3
     type(remap_plan_3D_real64), pointer ::remap_plan_x3_x1x2

     !--> cubic splines interpolation
    type(sll_cubic_spline_2d), pointer :: interp_x1x2
    type(sll_cubic_spline_1d), pointer :: interp_x3
    type(sll_cubic_spline_1d), pointer :: interp_x4



     !--> temporary structures that are used in CG_polar
     !type(sll_SL_polar), pointer :: plan_sl_polar


   contains
     procedure, pass(sim) :: run => run_dk4d_polar
     procedure, pass(sim) :: init_from_file => init_dk4d_polar
  end type sll_simulation_4d_drift_kinetic_polar

  interface delete
     module procedure delete_dk4d_polar
  end interface delete

contains

  subroutine init_dk4d_polar( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
    sll_int32            :: IO_stat
    sll_int32, parameter :: input_file = 99

    !--> Mesh
    sll_int32  :: num_cells_x1
    sll_int32  :: num_cells_x2
    sll_int32  :: num_cells_x3
    sll_int32  :: num_cells_x4
    sll_real64 :: r_min
    sll_real64 :: r_max
    sll_real64 :: z_min
    sll_real64 :: z_max
    sll_real64 :: v_min
    sll_real64 :: v_max
    !--> Equilibrium
    sll_real64 :: tau0
    sll_real64 :: rho_peak    
    sll_real64 :: kappan   
    sll_real64 :: deltarn  
    sll_real64 :: kappaTi  
    sll_real64 :: deltarTi 
    sll_real64 :: kappaTe  
    sll_real64 :: deltarTe
    sll_int32  :: QN_case
    !--> Pertubation
    sll_int32  :: perturb_choice
    sll_int32  :: mmode
    sll_int32  :: nmode
    sll_real64 :: eps_perturb   
    !--> Algorithm
    sll_real64 :: dt
    sll_int32  :: number_iterations
    sll_int32  :: carac_case
    sll_int32  :: time_case    
    !sll_int32  :: spline_degree
    
    !--> temporary variables for using cg_polar structures
    sll_int32  :: bc_cg(2)
    sll_int32  :: grad_cg
    sll_int32  :: carac_cg
    

    namelist /mesh/ num_cells_x1, num_cells_x2, &
      num_cells_x3, num_cells_x4, &
      r_min, r_max, z_min, z_max, &
      v_min, v_max
    namelist /equilibrium/ tau0, rho_peak, kappan, deltarn, &
      kappaTi, deltarTi, kappaTe, deltarTe, QN_case
    namelist /perturbation/ perturb_choice, mmode, nmode, eps_perturb
    namelist /sim_params/ dt, number_iterations, carac_case, time_case
      !, spline_degree

    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, '#init_dk4d_polar() failed to open file ', filename
       STOP
    end if
    read(input_file,mesh)
    read(input_file,equilibrium)
    read(input_file,perturbation)
    read(input_file,sim_params)
    close(input_file)

    !--> Mesh
    sim%m_x1 => new_logical_mesh_1d(num_cells_x1,eta_min=r_min,eta_max=r_max)
    sim%m_x2 => new_logical_mesh_1d(num_cells_x2,&
      eta_min=0._f64,eta_max=2._f64*sll_pi)
    sim%m_x3 => new_logical_mesh_1d(num_cells_x3,eta_min=z_min,eta_max=z_max)
    sim%m_x4 => new_logical_mesh_1d(num_cells_x4,eta_min=v_min,eta_max=v_max)
    
    !--> Equilibrium
    sim%tau0     = tau0
    sim%rho_peak = rho_peak 
    sim%kappan   = kappan
    sim%deltarn  = deltarn
    sim%kappaTi  = kappaTi
    sim%deltarTi = deltarTi
    sim%kappaTe  = kappaTe
    sim%deltarTe = deltarTe
    
    select case (QN_case)
      case (0)
        sim%QN_case = SLL_NO_QUASI_NEUTRAL
      case (1)
        sim%QN_case = SLL_QUASI_NEUTRAL_WITH_ZONAL_FLOW 
      case (2)
        sim%QN_case = SLL_QUASI_NEUTRAL_WITHOUT_ZONAL_FLOW 
      case default
        print *,'#bad choice for QN_case', QN_case
        print *,'#in init_dk4d_polar'
        stop
    end select

    select case (time_case)
      case (0)
        sim%time_case = SLL_TIME_LOOP_EULER
      case (1)
        sim%time_case = SLL_TIME_LOOP_PREDICTOR_CORRECTOR
      case default
        print *,'#bad choice for time_case', time_case
        print *,'#in init_dk4d_polar'
         stop
    end select
    select case (carac_case)
      case (0)
        sim%time_case = SLL_CARAC_EULER
      case (1)
        sim%time_case = SLL_CARAC_VERLET 
      case default
        print *,'#bad choice for carac_case', carac_case
        print *,'#in init_dk4d_polar'
         stop
    end select



    
    !--> Pertubation
    sim%perturb_choice = perturb_choice
    sim%mmode          = mmode
    sim%nmode          = nmode
    sim%eps_perturb    = eps_perturb
    !--> Algorithm
    sim%dt                 = dt
    sim%num_iterations     = number_iterations
    !sim%spline_degree_eta1 = spline_degree
    !sim%spline_degree_eta2 = spline_degree
    !sim%spline_degree_eta3 = spline_degree
    !sim%spline_degree_eta4 = spline_degree


    if(sll_get_collective_rank(sll_world_collective)==0)then
      print *,'##Mesh'
      print *,'#num_cells_x1=',num_cells_x1
      print *,'#num_cells_x2=',num_cells_x2
      print *,'#num_cells_x3=',num_cells_x3
      print *,'#num_cells_x4=',num_cells_x4
      print *,'#r_min=',r_min
      print *,'#r_max=',r_max
      print *,'#z_min=',z_min
      print *,'#z_max=',z_max
      print *,'#v_min=',v_min
      print *,'#v_max=',v_max
      print *,'##equilibrium'
      print *,'#tau0=',tau0
      print *,'#rho_peak=',rho_peak
      print *,'#kappan=',kappan
      print *,'#deltarn=',deltarn
      print *,'#kappaTi=',kappaTi
      print *,'#deltarTi=',deltarTi
      print *,'#kappaTe=',kappaTe
      print *,'#deltarTe=',deltarTe
      print *,'#QN_case=',QN_case
      print *,'##perturbation'
      print *,'#perturb_choice=',perturb_choice
      print *,'#mmode=',mmode
      print *,'#nmode=',nmode
      print *,'#eps_perturb=',eps_perturb
      print *,'#dt=',dt
      print *,'#number_iterations=',number_iterations
      print *,'#time_case=',time_case
      print *,'#carac_case=',carac_case
    endif
    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)

    
    call initialize_profiles_analytic(sim)    
    call allocate_fdistribu4d_DK(sim)
    call allocate_QN_DK( sim )
    
    sim%interp_x1x2 => new_spline_2d( sim%m_x1%num_cells+1,&
     sim%m_x2%num_cells+1,&
     sim%m_x1%eta_min,&
     sim%m_x1%eta_max,&
     sim%m_x2%eta_min,&
     sim%m_x2%eta_max,&
     SLL_HERMITE,&
     SLL_PERIODIC )

    sim%interp_x3 => new_spline_1d( sim%m_x3%num_cells+1,&
     sim%m_x3%eta_min,&
     sim%m_x3%eta_max,&
     SLL_PERIODIC )

    sim%interp_x4 => new_spline_1d( sim%m_x4%num_cells+1,&
     sim%m_x4%eta_min,&
     sim%m_x4%eta_max,&
     SLL_PERIODIC )


    
    

    grad_cg = 3  !3: splines 1: finite diff order 2 
    
    select case (sim%carac_case)
    
      case (SLL_CARAC_EULER)
        carac_cg = 1
      case (SLL_CARAC_VERLET)
        carac_cg = 4
      case default
        print *,'#bad value of sim%carac_case', sim%carac_case
        print *,'#in init_dk4d_polar'  
    end select
    
    bc_cg = (/SLL_DIRICHLET,SLL_DIRICHLET/)

    
!    sim%plan_sl_polar => &
!      new_SL(sim%m_x1%eta_min,&
!        sim%m_x1%eta_max,&
!        sim%m_x1%delta_eta,&
!        sim%m_x1%delta_eta,&
!        sim%dt,&
!        sim%m_x1%num_cells,&
!        sim%m_x2%num_cells,&
!        grad_cg,&
!        carac_cg,&
!        bc_cg)
!  
    !plan_poisson => new_plan_poisson_polar(geomx%dx,geomx%x0,geomx%nx-1,geomx%ny,bc,&
    !  &dlog_density,inv_Te)

    
    

  end subroutine init_dk4d_polar

  subroutine run_dk4d_polar(sim)
    class(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim
    !--> For initial profile HDF5 saving
    integer                      :: file_err
    sll_int32                    :: file_id
    character(len=12), parameter :: filename_prof = "init_prof.h5"
    sll_int32 :: loc4d_sz_x1
    sll_int32 :: loc4d_sz_x2
    sll_int32 :: loc4d_sz_x3
    sll_int32 :: loc4d_sz_x4

    !*** Saving of the radial profiles in HDF5 file ***
    if (sll_get_collective_rank(sll_world_collective)==0) then
      call sll_hdf5_file_create(filename_prof,file_id,file_err)
      call sll_hdf5_write_array_1d(file_id,sim%n0_r,'n0_r',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%Ti_r,'Ti_r',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%Te_r,'Te_r',file_err)
      call sll_hdf5_file_close(file_id,file_err)
    end if


    !*** Initialization of the distribution function ***
    !***  i.e f4d(t=t0)                              ***
    call initialize_fdistribu4d_DK(sim,sim%layout4d_x1x2x4,sim%f4d_x1x2x4)
!    call compute_reduction_4d_to_3d(&
!      sim%m_x4, &
!      sim%f4d_x1x2x4, &
!      sim%rho3d_x1x2, &
!      sim%layout4d_x1x2x4)
    call compute_local_sizes_4d( sim%layout4d_x1x2x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )


    call compute_reduction_4d_to_3d_direction4(&
      sim%f4d_x1x2x4, &
      sim%rho3d_x1x2, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4, &
      sim%m_x4%delta_eta)

 
 
    call apply_remap_3D( sim%remap_plan_x1x2_x3, sim%rho3d_x1x2, sim%rho3d_x3 )

    call solve_quasi_neutral(sim)

  end subroutine run_dk4d_polar


  subroutine delete_dk4d_polar( sim )
    class(sll_simulation_4d_drift_kinetic_polar) :: sim
    sll_int32 :: ierr
  end subroutine delete_dk4d_polar
  
  
  subroutine initialize_profiles_analytic(sim)
    class(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim
    sll_int32 :: i,ierr,nc_x1
    sll_real64 :: x1,delta_x1,rpeak,tmp,x1_min,x1_max
    
    nc_x1 = sim%m_x1%num_cells
    delta_x1 = sim%m_x1%delta_eta
    x1_min = sim%m_x1%eta_min
    x1_max = sim%m_x1%eta_max
    
    SLL_ALLOCATE(sim%n0_r(nc_x1+1),ierr)
    SLL_ALLOCATE(sim%Ti_r(nc_x1+1),ierr)
    SLL_ALLOCATE(sim%Te_r(nc_x1+1),ierr)
    
    rpeak = x1_min+sim%rho_peak*(x1_max-x1_min)
    do i=1,nc_x1+1
      x1 = x1_min+real(i-1,f64)*delta_x1
      sim%n0_r(i) = exp(-sim%kappan*sim%deltarn*tanh((x1-rpeak)/sim%deltarn))
      sim%Ti_r(i)=exp(-sim%kappaTi*sim%deltarTi*tanh((x1-rpeak)/sim%deltarTi))    
      sim%Te_r(i)=exp(-sim%kappaTe*sim%deltarTe*tanh((x1-rpeak)/sim%deltarTe))    
    enddo
    
    !we then change the normalization for n0_r
    tmp = 0.5_f64*(sim%n0_r(1)*x1_min+sim%n0_r(nc_x1+1)*x1_max)
    do i = 2,nc_x1
      x1 = x1_min+real(i-1,f64)*delta_x1
      tmp = tmp + sim%n0_r(i)*x1
    enddo
    tmp = tmp/real(nc_x1,f64)
    sim%n0_r = sim%n0_r/tmp
    sim%n0_at_rpeak = 1._f64/tmp      
  
  end subroutine initialize_profiles_analytic
  



  !----------------------------------------------------
  ! Allocation of the distribution function for
  !   drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine allocate_fdistribu4d_DK( sim )
    class(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim

    sll_int32 :: ierr, itemp
    sll_int32 :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4

    ! layout for sequential operations in x3 and x4. 
    ! Make an even split for x1 and x2, or as close as 
    ! even if the power of 2 is odd. This should 
    ! be packaged in some sort of routine and set up 
    ! at initialization time.
    sim%nproc_x1 = 1
    sim%nproc_x2 = 1
    sim%nproc_x3 = sim%world_size
    sim%nproc_x4 = 1
   
    
    
    sim%power2 = int(log(real(sim%world_size))/log(2.0))

    !--> Initialization of parallel layout of f4d in (x3,x4) directions
    !-->  (x1,x2) : sequential
    !-->  (x3,x4) : parallelized layout
    sim%layout4d_x1x2x4  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells+1, & 
      sim%m_x3%num_cells+1, &
      sim%m_x4%num_cells+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%layout4d_x1x2x4 )
    
    ! Allocate the array needed to store the local chunk 
    ! of the distribution function data. First compute the 
    ! local sizes. Since the remap operations
    ! are out-of-place, we will allocate two different arrays, 
    ! one for each layout.
    call compute_local_sizes_4d( sim%layout4d_x1x2x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
    SLL_ALLOCATE(sim%f4d_x1x2x4(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)

    !--> Initialization of parallel layout of f4d in (x1,x2,x4) directions
    !-->  (x1,x2,x4) : parallelized layout
    !-->  (x3) : sequential
    
    sim%nproc_x1 = 2**(sim%power2/3)
    sim%nproc_x2 = 2**(sim%power2/3)
    sim%nproc_x3 = 1
    sim%nproc_x4 = 2**(sim%power2-2*(sim%power2/3))
     

    sim%layout4d_x3  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells+1, & 
      sim%m_x3%num_cells+1, &
      sim%m_x4%num_cells+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%layout4d_x3 )
        
    call compute_local_sizes_4d( sim%layout4d_x3, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )    
    SLL_ALLOCATE(sim%f4d_x3(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)
    
    
    sim%remap_plan_x1x2x4_x3 => NEW_REMAP_PLAN(sim%layout4d_x1x2x4, &
      sim%layout4d_x3, sim%f4d_x1x2x4)
    sim%remap_plan_x3_x1x2x4 => NEW_REMAP_PLAN(sim%layout4d_x3, &
      sim%layout4d_x1x2x4, sim%f4d_x3)

    
    
  end subroutine allocate_fdistribu4d_DK


  !----------------------------------------------------
  ! Initialization of the distribution function for
  !   drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine initialize_fdistribu4d_DK(sim,layout,f4d)
    class(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim
    type(layout_4D), pointer :: layout
    sll_real64, dimension(:,:,:,:), pointer :: f4d
    sll_int32  :: ierr
    sll_int32  :: i1, i2, i3, i4
    sll_int32  :: iloc1, iloc2, iloc3, iloc4
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4
    sll_int32, dimension(1:4) :: glob_ind
    sll_real64, dimension(:), pointer :: x1_node,x2_node,x3_node,x4_node
    sll_real64 :: rpeak,k_x2,k_x3
    sll_real64 :: tmp_mode,tmp
    sll_real64 :: x1_min,x1_max
          

    !--> Initialization of the equilibrium distribution function
    SLL_ALLOCATE(sim%feq_x1x4(sim%m_x1%num_cells+1,sim%m_x4%num_cells+1),ierr)
    
    
    x1_min = sim%m_x1%eta_min
    x1_max = sim%m_x1%eta_max
        
    call initialize_eta1_node_1d(sim%m_x1,x1_node)
    call initialize_eta1_node_1d(sim%m_x2,x2_node)
    call initialize_eta1_node_1d(sim%m_x3,x3_node)
    call initialize_eta1_node_1d(sim%m_x4,x4_node)
    
    call init_fequilibrium( sim%m_x1%num_cells+1, &
      sim%m_x4%num_cells+1, &
      x1_node, &
      x4_node, &
      sim%n0_r, &
      sim%Ti_r, &
      sim%feq_x1x4 )

    !--> Initialization of the distribution function f4d_x3x4
    call compute_local_sizes_4d( layout, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
   
    k_x2  = 2._f64*sll_pi/(sim%m_x2%eta_max - sim%m_x2%eta_min)
    k_x3  = 2._f64*sll_pi/(sim%m_x3%eta_max - sim%m_x3%eta_min)
      
    rpeak = x1_min+sim%rho_peak*(x1_max-x1_min) 
    
    do iloc4 = 1,loc4d_sz_x4
      do iloc3 = 1,loc4d_sz_x3
        do iloc2 = 1,loc4d_sz_x2
          do iloc1 = 1,loc4d_sz_x1
            glob_ind(:) = local_to_global_4D(layout, &
              (/iloc1,iloc2,iloc3,iloc4/))
            i1 = glob_ind(1)
            i2 = glob_ind(2)
            i3 = glob_ind(3)
            i4 = glob_ind(4)
            tmp_mode = cos(real(sim%nmode,f64)*k_x3*x3_node(i3)&
               +real(sim%mmode,f64)*k_x2*x2_node(i2))
            tmp = exp(-(x1_node(i1)-rpeak)**2/(4._f64*sim%deltarn/sim%deltarTi))   
            f4d(iloc1,iloc2,iloc3,iloc4) = &
              (1._f64+tmp_mode*sim%eps_perturb*tmp)*sim%feq_x1x4(i1,i4)            
          end do
        end do
      end do
    end do
    SLL_DEALLOCATE(x1_node,ierr)
    SLL_DEALLOCATE(x2_node,ierr)
    SLL_DEALLOCATE(x3_node,ierr)
    SLL_DEALLOCATE(x4_node,ierr)
  end subroutine initialize_fdistribu4d_DK


  
  function compute_equil_analytic(sim,x1,x4)
    class(sll_simulation_4d_drift_kinetic_polar), intent(in) :: sim
    sll_real64,intent(in)::x1,x4
    sll_real64::compute_equil_analytic
    sll_real64:: tmp(2),rpeak,x1_min,x1_max
    x1_min = sim%m_x1%eta_min
    x1_max = sim%m_x1%eta_max
    
    rpeak = x1_min+sim%rho_peak*(x1_max-x1_min)
    tmp(1) = sim%n0_at_rpeak*exp(-sim%kappan*sim%deltarn*tanh((x1-rpeak)/sim%deltarn))
    tmp(2) = exp(-sim%kappaTi*sim%deltarTi*tanh((x1-rpeak)/sim%deltarTi))  
    compute_equil_analytic = tmp(1)/sqrt(2._f64*sll_pi*tmp(2))*exp(-0.5_f64*x4**2/tmp(2))
  
  end function compute_equil_analytic



  !----------------------------------------------------
  ! Allocation for QN solver
  !----------------------------------------------------
  subroutine allocate_QN_DK( sim )
    class(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim

    type(sll_logical_mesh_2d), pointer :: logical_mesh2d
    sll_int32 :: ierr, itemp
    !sll_int32 :: i1, i2, i3, i4
    !sll_int32 :: iloc1, iloc2, iloc3, iloc4
    sll_int32 :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3
    sll_int32 :: nproc3d_x3


    ! layout for sequential operations in x3 
    sim%power2 = int(log(real(sim%world_size))/log(2.0))
    !--> special case N = 1, so power2 = 0
    if(sim%power2 == 0) then
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 1
       sim%nproc_x4 = 1
    end if
    
    if(is_even(sim%power2)) then
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 2**(sim%power2/2)
       sim%nproc_x4 = 2**(sim%power2/2)
    else 
       sim%nproc_x1 = 1
       sim%nproc_x2 = 1
       sim%nproc_x3 = 2**((sim%power2-1)/2)
       sim%nproc_x4 = 2**((sim%power2+1)/2)
    end if

    !--> Initialization of rho3d_x1x2 and phi3d_x1x2
    !-->  (x1,x2) : sequential
    !-->  x3 : parallelized layout    
    sim%layout3d_x1x2  => new_layout_3D( sll_world_collective )
    nproc3d_x3 = sim%nproc_x3*sim%nproc_x4
    call initialize_layout_with_distributed_3D_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells+1, & 
      sim%m_x3%num_cells+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      nproc3d_x3, &
      sim%layout3d_x1x2 )
    call compute_local_sizes_3d( sim%layout3d_x1x2, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    SLL_ALLOCATE(sim%rho3d_x1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_x1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)

    !--> Initialization of rho3d_x3 and phi3d_x3
    !-->  (x1,x2) : parallelized layout
    !-->  x3 : sequential
    ! switch x3 and x1:
    itemp        = sim%nproc_x1
    sim%nproc_x1 = sim%nproc_x3
    sim%nproc_x3 = itemp
    ! switch x4 and x2
    itemp        = sim%nproc_x2
    sim%nproc_x2 = sim%nproc_x4 
    sim%nproc_x4 = itemp
        
    sim%layout3d_x3  => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells+1, & 
      sim%m_x3%num_cells+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%layout3d_x3 )
    call compute_local_sizes_3d( sim%layout3d_x3, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    SLL_ALLOCATE(sim%rho3d_x3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_x3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    
    
    sim%remap_plan_x1x2_x3 => NEW_REMAP_PLAN(sim%layout3d_x1x2, &
      sim%layout3d_x3, sim%rho3d_x1x2)
    sim%remap_plan_x3_x1x2 => NEW_REMAP_PLAN(sim%layout3d_x3, &
      sim%layout3d_x1x2, sim%rho3d_x3)

    
    
    
    
!
!    !----->
!    sim%rho2d => new_scalar_field_2d_discrete_alt( &
!      "rho2d_x1x2", &
!      sim%interp_rho2d, &     
!      sim%transf_xy, &
!      sim%bc_left_eta1, &
!      sim%bc_right_eta1, &
!      sim%bc_left_eta2, &
!      sim%bc_right_eta2)
!
!    sim%phi2d => new_scalar_field_2d_discrete_alt( &
!      "phi2d_x1x2", &
!      sim%interp_phi2d, &     
!      sim%transf_xy, &
!      sim%bc_left_eta1, &
!      sim%bc_right_eta1, &
!      sim%bc_left_eta2, &
!      sim%bc_right_eta2)
  end subroutine allocate_QN_DK

 
 
  
  subroutine solve_quasi_neutral(sim)
    class(sll_simulation_4d_drift_kinetic_polar), intent(inout) :: sim
    sll_int32 :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3
    sll_int32 :: iloc1, iloc2, iloc3
    sll_int32 :: i1, i2, i3
    sll_real64 :: tmp
    sll_int32 :: glob_ind(3)    
    call compute_local_sizes_3d( sim%layout3d_x3, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3 )



    
    select case (sim%QN_case)
      case (SLL_NO_QUASI_NEUTRAL)
      ! no quasi neutral solver as in CRPP-CONF-2001-069
        
        if((loc3d_sz_x3).ne.(sim%m_x3%num_cells+1))then
          print *,'#Problem of parallelization dimension in solve_quasi_neutral'
          print *,'#sll_simulation_4d_drift_kinetic_polar type simulation'
          stop
        endif
        
        do iloc2 = 1, loc3d_sz_x2
          do iloc1 = 1, loc3d_sz_x1          
            tmp = sum(sim%rho3d_x3(iloc1,iloc2,1:sim%m_x3%num_cells))&
              /real(sim%m_x3%num_cells,f64)
            SLL_ASSERT(loc3d_sz_x3==sim%m_x3%num_cells+1)
            do i3 = 1,sim%m_x3%num_cells+1
              glob_ind(:) = local_to_global_3D(sim%layout3d_x3, &
                (/iloc1,iloc2,i3/))                        
              sim%phi3d_x3(iloc1,iloc2,i3) = (sim%rho3d_x3(iloc1,iloc2,i3)-tmp)&
                *sim%Te_r(glob_ind(1))/sim%n0_r(glob_ind(1))
            enddo    
          enddo
        enddo  
        call apply_remap_3D( sim%remap_plan_x3_x1x2, sim%phi3d_x3, sim%phi3d_x1x2 )  
      case (SLL_QUASI_NEUTRAL_WITHOUT_ZONAL_FLOW)
        print *,'#SLL_QUASI_NEUTRAL_WITHOUT_ZONAL_FLOW'
        print *,'not implemented yet '
        stop
      case (SLL_QUASI_NEUTRAL_WITH_ZONAL_FLOW)
        print *,'#SLL_QUASI_NEUTRAL_WITH_ZONAL_FLOW'
        print *,'not implemented yet '
        stop      
      case default
        print *,'#bad value for sim%QN_case'
        stop  
    end select        
  
  end subroutine solve_quasi_neutral
  
  


!
!  subroutine compute_characteristics2D_verlet( A1,&
!    A2,&
!    dt, &
!    input1,&
!    intput2,&
!    output1,&
!    output2,&
!    Npts1,&
!    Npts2,&
!    interp1,&
!    interp2,&
!    maxiter,&
!    tol_input_x1,&
!    tol_input_x2)
!    
!    sll_real64, dimension(:,:), intent(in) :: A1
!    sll_real64, dimension(:,:), intent(in) :: A2
!    sll_real64, intent(in) :: dt
!    sll_real64, dimension(:), intent(in) ::  input1
!    sll_real64, dimension(:), intent(in) ::  input2
!    sll_int32, intent(in) :: Npts1    
!    sll_int32, intent(in) :: Npts2    
!    sll_real64, dimension(:,:), intent(out) :: output1
!    sll_real64, dimension(:,:), intent(out) :: output2
!    class(sll_interpolator_2d_base), pointer :: interp1
!    class(sll_interpolator_2d_base), pointer :: interp2
!    sll_int32,intent(in), optional :: maxiter_input
!    sll_real64,intent(in), optional :: tol_input_x1     
!    sll_real64,intent(in), optional :: tol_input_x2     
!    sll_int32 :: i
!    sll_int32 :: j
!    sll_int32 :: maxiter
!    sll_real64 :: tol_x1
!    sll_real64 :: tol_x2
!    sll_real64 :: x1
!    sll_real64 :: x2
!    sll_int32 :: iter
!    sll_real64 :: x1_old
!    
!    maxiter = 1000
!    tol_input = 1.e-12_f64
!    
!    
!    if(present(maxiter_input))then
!      maxiter = maxiter_input
!    endif
!    if(present(tol_input_x1))then
!      tol_x1 = tol_input_x1
!    endif
!    if(present(tol_input_x2))then
!      tol_x2 = tol_input_x2
!    endif
!    
!    
!    SLL_ASSERT(size(A1,1)>=Npts1)
!    SLL_ASSERT(size(A1,2)>=Npts2)
!    SLL_ASSERT(size(A2,1)>=Npts1)
!    SLL_ASSERT(size(A2,2)>=Npts2)
!    SLL_ASSERT(size(input1)>=Npts1)
!    SLL_ASSERT(size(input2)>=Npts2)
!    SLL_ASSERT(size(output1,1)>=Npts1)
!    SLL_ASSERT(size(output1,2)>=Npts2)
!    SLL_ASSERT(size(output2,1)>=Npts1)
!    SLL_ASSERT(size(output2,2)>=Npts2)
!    
!    interp1%compute_interpolants(A1)
!    interp2%compute_interpolants(A2)
!    
!    
!    do j=1,Npts2
!      do i=1,Npts1
!        !initialization for x1 interpolation
!        x1 = input1(i)-0.5_f64*dt*A1(i,j)
!        x1_old = 0._f64
!        do while (iter<maxiter .and. abs(x1_old-x1)>tol_x1)
!          x1_old = x1
!          x1 = input1(i)-0.5_f64*dt*interpolate_value(interp1, x1, input2(j))
!        end do
!        if (iter==maxiter .and. abs(x1_old-x1)>tol_x1) then
!          print*,'#not enough iterations for compute_characteristics2D_verlet',iter,abs(x1_old-x1)
!          stop
!        end if
!        !initialization for x2 interpolation
!        x2 = input2(j)-dt*A2(i,j)
!        x2_old = 0._f64
!        do while (iter<maxiter .and. abs(x2_old-x2)>tol_x2)
!          x2_old = x2
!          x2 = input2(j)-0.5_f64*dt*(interpolate_value(interp2, x1, x2)+interpolate_value(interp2, x1, input2(j)))
!        end do
!        if (iter==maxiter .and. abs(x2_old-x2)>tol_x2) then
!          print*,'#not enough iterations for compute_characteristics2D_verlet',iter,abs(x2_old-x2)
!          stop
!        end if
!        !initialization for x1 interpolation
!        x1 = x1-0.5_f64*dt*interpolate_value(interp1, x1, x2)
!        output1(i,j) = x1 
!        output2(i,j) = x2  
!      enddo
!    enddo        
!       
!  end subroutine compute_characteristics2D_verlet
!



  

end module sll_simulation_4d_drift_kinetic_polar_module



