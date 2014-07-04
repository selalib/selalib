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


module sll_simulation_4d_drift_kinetic_field_aligned_polar_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
  use sll_collective
  use sll_remapper
  use sll_constants
  use sll_test_4d_initializer
  use sll_poisson_2d_periodic_cartesian_par
  use sll_cubic_spline_interpolator_1d
  use sll_simulation_base
  use sll_fdistribu4D_DK
  use sll_logical_meshes
  use polar_operators
  use polar_advection
  use sll_reduction_module
  use sll_module_advection_2d_BSL
  use sll_module_characteristics_2d_explicit_euler
  use sll_module_characteristics_2d_verlet
  use sll_cubic_spline_interpolator_2d
  use sll_module_advection_1d_periodic
  use sll_module_poisson_2d_polar_solver
  use sll_qn_solver_3d_polar_parallel_x1_wrapper_module
  use sll_fcisl_module


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
  sll_int32, parameter :: SLL_TIME_LOOP_PREDICTOR2_CORRECTOR = 2 




  
  type, extends(sll_simulation_base_class) :: &
    sll_simulation_4d_drift_kinetic_field_aligned_polar

     ! Parallel environment parameters
     sll_int32  :: world_size
     sll_int32  :: my_rank
     sll_int32  :: power2 ! 2^power2 = number of processes available
     ! Processor mesh sizes
     sll_int32  :: nproc_x1
     sll_int32  :: nproc_x2
     sll_int32  :: nproc_x3
     sll_int32  :: nproc_x4 
     type(sll_collective_t), pointer :: new_collective_per_locx3
     type(sll_collective_t), pointer :: new_collective_per_locx4

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
     sll_int32  :: freq_diag_time
     sll_int32  :: freq_diag
     sll_int32  :: time_case
     sll_int32  :: charac_case
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
     sll_real64, dimension(:)  , pointer :: dlog_density_r

     !--> Magnetic field
     !sll_real64 :: q0
     !sll_real64 :: Dr_q0
     sll_real64, dimension(:), pointer :: iota_r
     sll_real64, dimension(:), pointer :: Diota_r
     sll_real64 :: B_norm_exponent
     sll_real64, dimension(:), pointer :: B_norm_r
     sll_real64, dimension(:), pointer :: Bstar_par_v_r
     sll_real64, dimension(:), pointer :: c_r
     sll_real64, dimension(:), pointer :: sigma_r
     sll_real64, dimension(:), pointer :: tau_r
     sll_real64, dimension(:), pointer :: iota_for_sigma
     sll_real64, dimension(:), pointer :: iota_for_tau
     


     !--> Equilibrium distribution function
     sll_real64, dimension(:,:), pointer :: feq_x1x4




     !----> parallel in x1
     type(layout_4D), pointer :: layout4d_parx1
     sll_real64, dimension(:,:,:,:), pointer :: f4d_parx1

     !----> parallel in (x3,x4)
     type(layout_4D), pointer :: layout4d_parx3x4
     sll_real64, dimension(:,:,:,:), pointer :: f4d_parx3x4

     !----> definition of remap
     type(remap_plan_4D_real64), pointer ::remap_plan_parx1_to_parx3x4
     type(remap_plan_4D_real64), pointer ::remap_plan_parx3x4_to_parx1

     !----> parallel in x1
     type(layout_3D), pointer :: layout3d_parx1
     sll_real64, dimension(:,:,:), pointer :: rho3d_parx1 
     sll_real64, dimension(:,:,:), pointer :: phi3d_parx1
     sll_real64, dimension(:,:,:), pointer :: Daligned_phi3d_parx1
     sll_real64, dimension(:,:,:), pointer :: A3_parx1 
     !----> parallel in x3
     type(layout_3D), pointer :: layout3d_parx3
     sll_real64, dimension(:,:,:), pointer :: rho3d_parx3
     sll_real64, dimension(:,:,:), pointer :: phi3d_parx3
     sll_real64, dimension(:,:,:), pointer :: A1_parx3
     sll_real64, dimension(:,:,:), pointer :: A2_parx3
     type(remap_plan_3D_real64), pointer ::remap_plan_parx1_to_parx3

     
     !----> for Poisson 
     type(layout_2D), pointer :: layout2d_parx1
     type(layout_2D), pointer :: layout2d_parx2

     
     
     
     

     !--> cubic splines interpolation
    !type(sll_cubic_spline_2d), pointer :: interp_x1x2
    !type(sll_cubic_spline_1d), pointer :: interp_x3
    !type(sll_cubic_spline_1d), pointer :: interp_x4

    sll_real64, dimension(:), pointer :: x1_node
    sll_real64, dimension(:), pointer :: x2_node
    sll_real64, dimension(:), pointer :: x3_node
    sll_real64, dimension(:), pointer :: x4_node





    class(sll_advection_2d_base), pointer :: adv_x1x2
    !class(sll_interpolator_2d_base), pointer :: interp_x1x2
    class(sll_characteristics_2d_base), pointer :: charac_x1x2
    class(sll_advection_1d_base), pointer :: adv_x3
    class(sll_advection_1d_base), pointer :: adv_x4


    class(sll_poisson_2d_base), pointer   :: poisson2d
    class(sll_poisson_2d_base), pointer   :: poisson2d_mean
    class(sll_poisson_3d_base), pointer :: poisson3d

    !for computing advection field from phi
    class(sll_interpolator_2d_base), pointer   :: phi_interp_x1x2
    class(sll_interpolator_1d_base), pointer   :: phi_interp_x3


     !--> temporary structures that are used in CG_polar
     !type(sll_SL_polar), pointer :: plan_sl_polar


   contains
     procedure, pass(sim) :: run => run_dk4d_field_aligned_polar
     procedure, pass(sim) :: init_from_file => init_dk4d_field_aligned_polar
  end type sll_simulation_4d_drift_kinetic_field_aligned_polar
  
   
  interface delete
     module procedure delete_dk4d_field_aligned_polar
  end interface delete

contains






!we should not give directly the file here
!but a long list of parameters that would be initialized with
!a read_from_file routine

  subroutine init_dk4d_field_aligned_polar( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename
    sll_int32            :: IO_stat
    sll_int32, parameter :: input_file = 99
    class(sll_characteristics_2d_base), pointer :: charac2d
    class(sll_interpolator_2d_base), pointer   :: A1_interp2d
    class(sll_interpolator_2d_base), pointer   :: A2_interp2d
    class(sll_interpolator_1d_base), pointer   :: A1_interp1d_x1
    class(sll_interpolator_1d_base), pointer   :: A2_interp1d_x1
    class(sll_interpolator_2d_base), pointer   :: f_interp2d
    sll_real64 :: charac2d_tol
    sll_int32 :: charac2d_maxiter




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
    sll_real64 :: iota0
    sll_real64 :: Dr_iota0
    character(len=256) :: iota_file
    character(len=256) :: Diota_file
    sll_int32 :: size_iota_file
    sll_int32 :: size_Diota_file
    logical :: is_iota_file
    logical :: is_Diota_file

    sll_real64 :: B_norm_exponent
    !sll_int32  :: QN_case
    !--> Pertubation
    sll_int32  :: perturb_choice
    sll_int32  :: mmode
    sll_int32  :: nmode
    sll_real64 :: eps_perturb   
    !--> Algorithm
    sll_real64 :: dt
    sll_int32  :: number_iterations
    sll_int32  :: freq_diag_time
    sll_int32  :: freq_diag
    !sll_int32  :: charac_case
    !sll_int32  :: time_case    
    character(len=256)      :: advect2d_case 
    character(len=256)      :: charac2d_case
    !character(len=256)      :: f_interp2d_case 
    !character(len=256)      :: phi_interpx1x2
    !character(len=256)      :: phi_interpx3
    !character(len=256)      :: A_interp_case 
    !character(len=256)      :: initial_function_case 
    character(len=256)      :: time_loop_case 
    character(len=256)      :: poisson2d_case 
    character(len=256)      :: QN_case 
    character(len=256)      :: advector_x3 
    character(len=256)      :: advector_x4
    character(len=256)      :: interp_x1x2
    character(len=256)      :: phi_interp_x1x2
    character(len=256)      :: phi_interp_x3
    character(len=256)      :: poisson2d_BC_rmin
    character(len=256)      :: poisson2d_BC_rmax
    
    sll_int32               :: order_x3
    sll_int32               :: order_x4
    sll_int32               :: poisson2d_BC(2)
    sll_real64, dimension(:,:), allocatable :: tmp_r
    sll_int32 :: i
    sll_int32 :: ierr
    !sll_int32  :: spline_degree
    
    
    
    
    !--> temporary variables for using cg_polar structures
    !sll_int32  :: bc_cg(2)
    !sll_int32  :: grad_cg
    !sll_int32  :: carac_cg
    

    namelist /mesh/ &
      num_cells_x1, &
      num_cells_x2, &
      num_cells_x3, &
      num_cells_x4, &
      r_min, &
      r_max, &
      z_min, &
      z_max, &
      v_min, &
      v_max
    namelist /equilibrium/ & 
      tau0, &
      rho_peak, &
      kappan, &
      deltarn, &
      kappaTi, &
      deltarTi, &
      kappaTe, &
      deltarTe, &
      QN_case, &
      poisson2d_BC_rmin, &
      poisson2d_BC_rmax, &
      iota0, &
      Dr_iota0, &
      iota_file, &
      Diota_file, &
      size_iota_file, &
      size_Diota_file, &
      is_iota_file, &
      is_Diota_file, &
      B_norm_exponent      
    namelist /perturbation/ &
      perturb_choice, &
      mmode, &
      nmode, &
      eps_perturb
    namelist /sim_params/ &
      dt, & 
      number_iterations, &
      freq_diag_time, &
      freq_diag, &
      charac2d_case, &
      time_loop_case, &
      advect2d_case, &
      charac2d_tol, &
      charac2d_maxiter, &
      interp_x1x2, &
      phi_interp_x1x2, &
      phi_interp_x3, &
      advector_x3, &
      advector_x4, &
      order_x3, &
      order_x4, &
      poisson2d_case
      
      !, spline_degree
    
    !default parameters
    
    iota0 = 0._f64
    Dr_iota0 = 0._f64
    iota_file = "no_q_file.dat"
    Diota_file = "no_q_file.dat"
    is_iota_file = .false.
    is_Diota_file = .false.
    size_iota_file = 0
    size_Diota_file = 0
    
    
    B_norm_exponent = -0.5_f64
    
     
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
    !sim%q0 = q0
    !sim%Dr_q0 = Dr_q0
    sim%B_norm_exponent = B_norm_exponent

    
    SLL_ALLOCATE(tmp_r(num_cells_x1+1,2),ierr)



 !   call compute_spaghetti_and_shift_from_guess( &
!      Nc_x1, &
!      Nc_x2, &
!      iota, &
!      spaghetti_size_guess, &
!      shift, &
!      spaghetti_size)
!
!    print *,'#shift=',shift
!    print *,'#spaghetti_size=',spaghetti_size
!    call compute_iota_from_shift(Nc_x1,shift,iota_modif)
!    print *,'#iota_modif=',iota_modif

    
    
    select case (poisson2d_BC_rmin)
      case ("SLL_DIRICHLET")
        poisson2d_BC(1) = SLL_DIRICHLET
      case ("SLL_NEUMANN")
        poisson2d_BC(1) = SLL_NEUMANN
      case ("SLL_NEUMANN_MODE_0")
        poisson2d_BC(1) = SLL_NEUMANN_MODE_0      
      case default
        print *,'#bad choice for poisson2d_BC_rmin'
        print *,'#in init_dk4d_polar'
        stop
    end select   


    select case (poisson2d_BC_rmax)
      case ("SLL_DIRICHLET")
        poisson2d_BC(2) = SLL_DIRICHLET
      case ("SLL_NEUMANN")
        poisson2d_BC(2) = SLL_NEUMANN
      case ("SLL_NEUMANN_MODE_0")
        poisson2d_BC(2) = SLL_NEUMANN_MODE_0      
      case default
        print *,'#bad choice for poisson2d_BC_rmax'
        print *,'#in init_dk4d_polar'
        stop
    end select   

    
    
    select case (QN_case)
      case ("SLL_NO_QUASI_NEUTRAL")
        sim%QN_case = SLL_NO_QUASI_NEUTRAL
      case ("SLL_QUASI_NEUTRAL_WITH_ZONAL_FLOW")
        sim%QN_case = SLL_QUASI_NEUTRAL_WITH_ZONAL_FLOW 
      case ("SLL_QUASI_NEUTRAL_WITHOUT_ZONAL_FLOW")
        sim%QN_case = SLL_QUASI_NEUTRAL_WITHOUT_ZONAL_FLOW 
      case default
        print *,'#bad choice for QN_case', QN_case
        print *,'#in init_dk4d_polar'
        stop
    end select

    select case (time_loop_case)
      case ("SLL_TIME_LOOP_EULER")
        sim%time_case = SLL_TIME_LOOP_EULER
      case ("SLL_TIME_LOOP_PREDICTOR_CORRECTOR")
        sim%time_case = SLL_TIME_LOOP_PREDICTOR_CORRECTOR
      case ("SLL_TIME_LOOP_PREDICTOR2_CORRECTOR")
        sim%time_case = SLL_TIME_LOOP_PREDICTOR2_CORRECTOR
      case default
        print *,'#bad choice for time_loop_case', time_loop_case
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
    sim%freq_diag_time     = freq_diag_time
    sim%freq_diag     = freq_diag
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
      print *,'#time_loop_case=',time_loop_case
      print *,'#charac2d_case=',charac2d_case
    endif
    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)

    call initialize_eta1_node_1d(sim%m_x1,sim%x1_node)
    call initialize_eta1_node_1d(sim%m_x2,sim%x2_node)
    call initialize_eta1_node_1d(sim%m_x3,sim%x3_node)
    call initialize_eta1_node_1d(sim%m_x4,sim%x4_node)

    SLL_ALLOCATE(sim%iota_r(num_cells_x1+1),ierr)
    SLL_ALLOCATE(sim%Diota_r(num_cells_x1+1),ierr)
   
    
    call initialize_iota_profile( &
      iota0, &    
      Dr_iota0, &
      iota_file, &
      Diota_file, &
      size_iota_file, &
      size_Diota_file, &
      is_iota_file, &
      is_Diota_file, &
      sim%m_x1%num_cells+1, &
      sim%x1_node, &
      sim%iota_r, &
      sim%Diota_r )


    
    call initialize_profiles_analytic(sim)    
    !call allocate_fdistribu4d_DK(sim)
    !call allocate_QN_DK( sim )

    call allocate_fdistribu4d_and_QN_DK_parx1(sim)

    
    
    
    
    select case (poisson2d_case)
      case ("POLAR_FFT")     
        
        do i=1,num_cells_x1+1
          tmp_r(i,1) = 1._f64/sim%Te_r(i)
        enddo  
        
        sim%poisson2d_mean =>new_poisson_2d_polar_solver( &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sim%m_x1%num_cells, &
          sim%m_x2%num_cells, &
          poisson2d_BC)

        sim%poisson2d =>new_poisson_2d_polar_solver( &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sim%m_x1%num_cells, &
          sim%m_x2%num_cells, &
          poisson2d_BC, &
          dlog_density=sim%dlog_density_r, &
          inv_Te=tmp_r(1:num_cells_x1+1,1), &
          poisson_case=SLL_POISSON_DRIFT_KINETIC)

        sim%poisson3d => new_qn_solver_3d_polar_parallel_x1_wrapper( &
          sim%layout2d_parx2, &
          sim%layout2d_parx1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sim%m_x1%num_cells, &
          sim%m_x2%num_cells, &
          sim%m_x3%num_cells, &
          poisson2d_BC(1), &
          poisson2d_BC(2), &
          dlog_density=sim%dlog_density_r, &
          inv_Te=tmp_r(1:num_cells_x1+1,1))




!        sim%poisson3d =>new_poisson_2d_polar_solver( &
!          sim%m_x1%eta_min, &
!          sim%m_x1%eta_max, &
!          sim%m_x1%num_cells, &
!          sim%m_x2%num_cells, &
!          poisson2d_BC, &
!          dlog_density=sim%dlog_density_r, &
!          inv_Te=tmp_r(1:num_cells_x1+1,1), &
!          poisson_case=SLL_POISSON_DRIFT_KINETIC)


          
      case default
        print *,'#bad poisson2d_case',poisson2d_case
        print *,'#not implemented'
        print *,'#in init_dk4d_polar'
        stop
    end select
    

!    select case (sim%QN_case)
!      case (SLL_NO_QUASI_NEUTRAL)
!      case (SLL_QUASI_NEUTRAL_WITHOUT_ZONAL_FLOW)
!      case (SLL_QUASI_NEUTRAL_WITH_ZONAL_FLOW)
!      case default
!        print *,'#bad value for sim%QN_case'
!        stop  
!    end select        



    select case (interp_x1x2)
      case ("SLL_CUBIC_SPLINES")
        f_interp2d => new_cubic_spline_2d_interpolator( &
          sim%m_x1%num_cells+1, &
          sim%m_x2%num_cells+1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sim%m_x2%eta_min, &
          sim%m_x2%eta_max, &
          SLL_HERMITE, &
          SLL_PERIODIC)
        A1_interp1d_x1 => new_cubic_spline_1d_interpolator( &
          sim%m_x1%num_cells+1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          SLL_HERMITE)
        A2_interp1d_x1 => new_cubic_spline_1d_interpolator( &
          sim%m_x1%num_cells+1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          SLL_HERMITE)
        A1_interp2d => new_cubic_spline_2d_interpolator( &
          sim%m_x1%num_cells+1, &
          sim%m_x2%num_cells+1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sim%m_x2%eta_min, &
          sim%m_x2%eta_max, &
          SLL_HERMITE, &
          SLL_PERIODIC)
        A2_interp2d => new_cubic_spline_2d_interpolator( &
          sim%m_x1%num_cells+1, &
          sim%m_x2%num_cells+1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sim%m_x2%eta_min, &
          sim%m_x2%eta_max, &
          SLL_HERMITE, &
          SLL_PERIODIC)
      case default
        print *,'#bad interp_x1x2',interp_x1x2
        print *,'#not implemented'
        print *,'#in init_dk4d_polar'
        stop
    end select


    select case (phi_interp_x1x2)
      case ("SLL_CUBIC_SPLINES")
         sim%phi_interp_x1x2 => new_cubic_spline_2d_interpolator( &
          sim%m_x1%num_cells+1, &
          sim%m_x2%num_cells+1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sim%m_x2%eta_min, &
          sim%m_x2%eta_max, &
          SLL_HERMITE, &
          SLL_PERIODIC)
      case default
        print *,'#bad phi_interp_x1x2',phi_interp_x1x2
        print *,'#not implemented'
        print *,'#in init_dk4d_polar'
        stop
    end select




    select case (phi_interp_x3)
      case ("SLL_CUBIC_SPLINES")
        sim%phi_interp_x3 => new_cubic_spline_1d_interpolator( &
          sim%m_x3%num_cells+1, &
          sim%m_x3%eta_min, &
          sim%m_x3%eta_max, &
          SLL_PERIODIC)
      case default
        print *,'#bad phi_interp_x3',phi_interp_x3
        print *,'#not implemented'
        print *,'#in init_dk4d_polar'
        stop
    end select




    
    select case (charac2d_case)   
      case("SLL_CHARAC_EULER")
        charac2d => new_explicit_euler_2d_charac(&
          sim%m_x1%num_cells+1, &
          sim%m_x2%num_cells+1, &
          bc_type_1=SLL_SET_TO_LIMIT, &
          bc_type_2=SLL_PERIODIC, &
          eta1_min = sim%m_x1%eta_min, &
          eta1_max = sim%m_x1%eta_max, &
          eta2_min = sim%m_x2%eta_min, &
          eta2_max = sim%m_x2%eta_max)
      case("SLL_CHARAC_VERLET")
        charac2d => new_verlet_2d_charac(&
          sim%m_x1%num_cells+1, &
          sim%m_x2%num_cells+1, &
          A1_interp2d, &
          A2_interp2d, &
          A1_interp1d_x1, &
          A2_interp1d_x1, &
          bc_type_1=SLL_SET_TO_LIMIT, &
          bc_type_2=SLL_PERIODIC, &
          eta1_min = sim%m_x1%eta_min, &
          eta1_max = sim%m_x1%eta_max, &
          eta2_min = sim%m_x2%eta_min, &
          eta2_max = sim%m_x2%eta_max, &
          x1_maxiter = charac2d_maxiter, &
          x2_maxiter = charac2d_maxiter, &
          x1_tol = charac2d_tol, &
          x2_tol = charac2d_tol)

      case default
        print *,'#bad choice for charac_case', charac2d_case
        print *,'#in init_dk4d_polar'
        print *,'#should be: SLL_CHARAC_EULER'
        print *,'#or: SLL_CHARAC_VERLET'
        stop
    end select


    select case(advect2d_case)
      case ("SLL_BSL")
      sim%adv_x1x2 => new_BSL_2d_advector(&
        f_interp2d, &
        charac2d, &
        sim%m_x1%num_cells+1, &
        sim%m_x2%num_cells+1, &
        eta1_min = sim%m_x1%eta_min, &
        eta1_max = sim%m_x1%eta_max, &
        eta2_min = sim%m_x2%eta_min, &
        eta2_max = sim%m_x2%eta_max)
      case default
        print *,'#bad advect_case',advect2d_case
        print *,'#not implemented'
        print *,'#in initialize_guiding_center_2d_polar'
        stop
    end select

      
    select case (advector_x3)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%adv_x3 => new_periodic_1d_advector( &
          sim%m_x3%num_cells, &
          sim%m_x3%eta_min, &
          sim%m_x3%eta_max, &
          SPLINE, & 
          order_x3) 
      case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation
        sim%adv_x3 => new_periodic_1d_advector( &
          sim%m_x3%num_cells, &
          sim%m_x3%eta_min, &
          sim%m_x3%eta_max, &
          LAGRANGE, & 
          order_x3) 
       case default
         print*,'#advector in x3', advector_x3, ' not implemented'
         stop 
    end select

    select case (advector_x4)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%adv_x4 => new_periodic_1d_advector( &
          sim%m_x4%num_cells, &
          sim%m_x4%eta_min, &
          sim%m_x4%eta_max, &
          SPLINE, & 
          order_x4) 
      case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation
        sim%adv_x4 => new_periodic_1d_advector( &
          sim%m_x4%num_cells, &
          sim%m_x4%eta_min, &
          sim%m_x4%eta_max, &
          LAGRANGE, & 
          order_x4) 
       case default
         print*,'#advector in x4', advector_x4, ' not implemented'
         stop 
    end select
        
  end subroutine init_dk4d_field_aligned_polar

  subroutine run_dk4d_field_aligned_polar(sim)
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar), intent(inout) :: sim
    !--> For initial profile HDF5 saving
    integer                      :: file_err
    sll_int32                    :: file_id
    character(len=12), parameter :: filename_prof = "init_prof.h5"
    sll_real64,dimension(:,:,:,:), allocatable :: f4d_store
    sll_int32 :: loc4d_sz_x1
    sll_int32 :: loc4d_sz_x2
    sll_int32 :: loc4d_sz_x3
    sll_int32 :: loc4d_sz_x4
    sll_int32 :: loc4d(4)
    sll_int32 :: iter
    sll_int32 :: nc_x1
    sll_int32 :: nc_x2
    sll_int32 :: nc_x3
    sll_int32 :: nc_x4
    !sll_int32 :: i2
    !sll_int32 :: i3
    !sll_int32 :: i4
    sll_int32 :: ierr
    sll_real64 :: dt
    sll_int32 :: th_diag_id
    sll_int32 :: i_plot 
    
    dt = sim%dt    
    
    nc_x1 = sim%m_x1%num_cells
    nc_x2 = sim%m_x2%num_cells
    nc_x3 = sim%m_x3%num_cells
    nc_x4 = sim%m_x4%num_cells
    


    !*** Saving of the radial profiles in HDF5 file ***
    if (sll_get_collective_rank(sll_world_collective)==0) then
      call sll_hdf5_file_create(filename_prof,file_id,file_err)
      call sll_hdf5_write_array_1d(file_id,sim%n0_r,'n0_r',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%Ti_r,'Ti_r',file_err)
      call sll_hdf5_write_array_1d(file_id,sim%Te_r,'Te_r',file_err)
      call sll_hdf5_file_close(file_id,file_err)
      
      call sll_gnuplot_write(sim%n0_r,'n0_r_init',ierr)
      call sll_gnuplot_write(sim%Ti_r,'Ti_r_init',ierr)
      call sll_gnuplot_write(sim%Te_r,'Te_r_init',ierr)
      call sll_gnuplot_write(sim%iota_r,'iota_r_init',ierr)
      call sll_gnuplot_write(sim%c_r,'c_r_init',ierr)
      call sll_gnuplot_write(sim%B_norm_r,'B_norm_r_init',ierr)
      call sll_gnuplot_write(sim%Bstar_par_v_r,'Bstar_par_v_r_init',ierr)
    end if


    call compute_local_sizes_4d( sim%layout4d_parx1, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
    SLL_ALLOCATE(f4d_store(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)

    if(sll_get_collective_rank(sll_world_collective)==0) then
      call sll_ascii_file_create('thdiag.dat', th_diag_id, ierr)
    endif




    call initialize_fdistribu4d_DK(sim,sim%layout4d_parx1,sim%f4d_parx1)




    i_plot = 0
        
    do iter=1,sim%num_iterations    


      call compute_local_sizes_4d( sim%layout4d_parx1, &
        loc4d_sz_x1, &
        loc4d_sz_x2, &
        loc4d_sz_x3, &
        loc4d_sz_x4 )
    
    !print *,'#locsize it',iter,sim%my_rank,loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4

   
      
      call compute_rho_dk(sim)    


      if(modulo(iter,sim%freq_diag)==0)then
        if(sll_get_collective_rank(sll_world_collective)==0) then
          print*,'#iteration=',iter      
        endif
      endif
      call solve_quasi_neutral_parx1( sim )
      call compute_field_dk_parx1( sim )

      if(modulo(iter,sim%freq_diag_time)==0)then
        call time_history_diagnostic_dk_polar( &
          sim, &
          th_diag_id, &    
          iter-1)
      endif            


      select case (sim%time_case)
        case (SLL_TIME_LOOP_EULER)
          call advection_x3( sim, dt )
          call advection_x4( sim, dt )
          call advection_x1x2( sim, dt )

        case (SLL_TIME_LOOP_PREDICTOR_CORRECTOR)
          !prediction
          f4d_store = sim%f4d_parx1
          call advection_x3( sim, 0.5_f64*dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x1x2( sim, 0.5_f64*dt )
          call compute_rho_dk(sim)  
          call solve_quasi_neutral_parx1( sim )
          call compute_field_dk_parx1( sim )          
          !correction
          sim%f4d_parx1 = f4d_store
          call advection_x3( sim, 0.5_f64*dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x1x2( sim, dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x3( sim, 0.5_f64*dt )
        case (SLL_TIME_LOOP_PREDICTOR2_CORRECTOR)
          !prediction order "2"
          f4d_store = sim%f4d_parx1
          call advection_x3( sim, 0.25_f64*dt )
          call advection_x4( sim, 0.25_f64*dt )
          call advection_x1x2( sim, 0.5_f64*dt )
          call advection_x4( sim, 0.25_f64*dt )
          call advection_x3( sim, 0.25_f64*dt )
          call compute_rho_dk(sim)  
          call solve_quasi_neutral_parx1( sim )
          call compute_field_dk_parx1( sim )          
          !correction
          sim%f4d_parx1 = f4d_store
          call advection_x3( sim, 0.5_f64*dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x1x2( sim, dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x3( sim, 0.5_f64*dt )
        case default
          call sll_halt_collective()
          print *,'#sim%time_case=',sim%time_case
          print *, '#not implemented'
          print *,'#in run_dk4d_polar'
          stop
      end select          

    

      if(modulo(iter,sim%freq_diag)==0) then

        i_plot = i_plot+1
        loc4d(1:4)  = global_to_local_4D( &
          sim%layout4d_parx1, &
          (/nc_x1/2+1,1,1,nc_x4/2+1/))
        if(loc4d(1) > 0) then
#ifndef NOHDF5
          call plot_f_cartesian( &
            i_plot, &
            sim%f4d_parx1(loc4d(1),:,:,loc4d(4)), &
            sim%m_x2,sim%m_x3)
#endif
        endif    


        call apply_remap_4D( &
          sim%remap_plan_parx1_to_parx3x4, &
          sim%f4d_parx1, &
          sim%f4d_parx3x4 )

        loc4d(1:4)  = global_to_local_4D( &
          sim%layout4d_parx3x4, &
          (/1,1,1,nc_x4/2+1/))
        if(loc4d(3) > 0) then
!          call sll_gnuplot_corect_2d( &
!            sim%m_x1%eta_min, &
!            sim%m_x1%eta_max, &
!            nc_x1+1, &
!            sim%m_x2%eta_min, &
!            sim%m_x2%eta_max, &
!            nc_x2+1, &
!            sim%f4d_parx3x4(:,:,loc4d(3),loc4d(4)), &
!            'fdist', &
!            i_plot, &
!            ierr)
#ifndef NOHDF5
          call plot_f_polar( &
            i_plot, &
            sim%f4d_parx3x4(:,:,loc4d(3),loc4d(4)), &
            sim%m_x1,sim%m_x2)
#endif
        endif

        call apply_remap_4D( &
          sim%remap_plan_parx3x4_to_parx1, &
          sim%f4d_parx3x4, &
          sim%f4d_parx1 )
          
                    
      endif



        

   enddo
   
   

  end subroutine run_dk4d_field_aligned_polar

!> initialize q_profile
!> if is_q_file=.false. and is_Dq_file=.false., takes q(r) = q0_r+r*Dq0_r
!> if is_q_file=.false. and is_Dq_file=.true., error
!> if is_q_file=.true. and is_Dq_file=.false., takes q(r) from q_file Dq(r) by interp of q 
!> if is_q_file=.true. and is_Dq_file=.false., takes q(r) from q_file Dq(r) from Dq_file 
!> when q (or Dq) are taken from q_file, we use cubic spline interpolation for getting
!> values on grid points
  
  subroutine initialize_iota_profile( &
    iota0_r, &    
    Diota0_r, &
    iota_file, &
    Diota_file, &
    size_iota_file, &
    size_Diota_file, &
    is_iota_file, &
    is_Diota_file, &
    num_points_r, &
    r_array, &
    iota, &
    Diota )
    sll_real64, intent(in) :: iota0_r
    sll_real64, intent(in) :: Diota0_r
    character(len=256), intent(in) :: iota_file
    character(len=256), intent(in) :: Diota_file
    sll_int32, intent(in) :: size_iota_file
    sll_int32, intent(in) :: size_Diota_file
    logical, intent(in) :: is_iota_file
    logical, intent(in) :: is_Diota_file
    sll_int32, intent(in) :: num_points_r
    sll_real64, dimension(:), intent(in) :: r_array
    sll_real64, dimension(:), intent(out) :: iota
    sll_real64, dimension(:), intent(out) :: Diota
    !local variables
    sll_int32 :: i
    
    ! some checking
    if(size(iota)<num_points_r)then
      call sll_halt_collective()
      print *,'#bad size for iota in initialize_iota_profile'
      stop
    endif
    if(size(Diota)<num_points_r)then
      call sll_halt_collective()
      print *,'#bad size for Diota in initialize_iota_profile'
      stop
    endif
    
    if ((is_iota_file .eqv. .false.) .and. (is_Diota_file .eqv. .false.)) then
      if(size(r_array)<num_points_r)then
        call sll_halt_collective()
        print *,'#bad size for r_array in initialize_iota_profile'
        stop
      endif
      do i=1,num_points_r
        iota(i) = iota0_r+Diota0_r*r_array(i)
        Diota(i) = Diota0_r        
      enddo
    endif


    if ((is_iota_file .eqv. .false.) .and. (is_Diota_file .eqv. .true.)) then
      call sll_halt_collective()
      print *,'#bad value for is_iota_file and is_Diota_file in initialize_iota_profile'
      stop      
    endif

    if ((is_iota_file .eqv. .true.) .and. (is_Diota_file .eqv. .false.))then
      call sll_halt_collective()
      print *,'#not implemented for the moment'
      print *,'#in initialize_iota_profile'
      stop
    endif

    if ((is_iota_file .eqv. .true.) .and. (is_Diota_file .eqv. .true.))then
      call sll_halt_collective()
      print *,'#not implemented for the moment'
      print *,'#in initialize_iota_profile'
      stop
    endif
    
    
    
  end subroutine initialize_iota_profile
  
  subroutine initialize_c_from_iota_profile( &
    iota, &
    r_array, &
    num_points_r, &
    L, &
    c_r)
    sll_real64, dimension(:), intent(in) :: iota
    sll_real64, dimension(:), intent(in) :: r_array
    sll_int32, intent(in) :: num_points_r
    sll_real64, intent(in) :: L
    sll_real64, dimension(:), intent(out) :: c_r
    !local variables
    sll_int32 :: i
    sll_real64 :: big_R

    ! some checking
    if(size(iota)<num_points_r)then
      call sll_halt_collective()
      print *,'#bad size for iota in initialize_c_from_iota_profile'
      stop
    endif

    if(size(r_array)<num_points_r)then
      call sll_halt_collective()
      print *,'#bad size for r_array in initialize_c_from_iota_profile'
      stop
    endif

    if(size(c_r)<num_points_r)then
      call sll_halt_collective()
      print *,'#bad size for c_r in initialize_c_from_iota_profile'
      stop
    endif
    
    big_R = L/(2._f64*sll_pi)
    
    do i=1,num_points_r
      c_r(i) = r_array(i)*iota(i)/big_R
    enddo
    
  end subroutine initialize_c_from_iota_profile

  subroutine initialize_iota_modif( &
    Nc_x1, &
    Nc_x2, &
    iota, &
    num_points_r, &
    spaghetti_size_guess, &
    spaghetti_size, &
    shift_r)
    sll_int32, intent(in) :: Nc_x1
    sll_int32, intent(in) :: Nc_x2
    sll_real64, dimension(:), intent(in) :: iota
    sll_int32, intent(in) :: num_points_r
    sll_int32, intent(in) :: spaghetti_size_guess
    sll_int32, intent(out) :: spaghetti_size
    sll_int32, dimension(:), intent(out) :: shift_r
    !local variables
    sll_int32 :: i
    sll_int32 :: spaghetti_size0

    if(size(iota)<num_points_r)then
      call sll_halt_collective()
      print *,'#bad size for iota in initialize_iota_modif'
      stop
    endif

    if(size(shift_r)<num_points_r)then
      call sll_halt_collective()
      print *,'#bad size for shift_R in initialize_iota_modif'
      stop
    endif

    do i=1,num_points_r
      call compute_spaghetti_and_shift_from_guess( &
        Nc_x1, &
        Nc_x2, &
        iota(i), &
        spaghetti_size_guess, &
        shift_r(i), &
        spaghetti_size)
      if(i==1)then
        spaghetti_size0 = spaghetti_size
      endif
      if(spaghetti_size .ne. spaghetti_size0)then
        print *,'#bad spaghetti size in initialize_iota_modif'
        print *,'#we want to have same spaghetti_size for all the r'
      endif  
    enddo
    

  end subroutine initialize_iota_modif



  subroutine time_history_diagnostic_dk_polar( &
    sim, &
    file_id, &    
    step)
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar) :: sim
    sll_int32, intent(in) :: file_id
    sll_int32, intent(in) :: step
    sll_real64 :: dt
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_int32 :: Nc_x3
    sll_int32 :: Nc_x4
    sll_real64 :: nrj
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: delta3
    sll_real64 :: delta4
    sll_int32 :: ierr
    dt = sim%dt

    Nc_x1 = sim%m_x1%num_cells
    Nc_x2 = sim%m_x2%num_cells
    Nc_x3 = sim%m_x3%num_cells
    Nc_x4 = sim%m_x4%num_cells

    delta1 = sim%m_x1%delta_eta
    delta2 = sim%m_x2%delta_eta
    delta3 = sim%m_x3%delta_eta
    delta4 = sim%m_x4%delta_eta



    nrj = 0._f64
    call compute_reduction_2d_to_0d(&
      sim%phi3d_parx3(:,:,1)**2, &
      nrj, &
      Nc_x1+1, &
      Nc_x2+1, &
      delta1, &    
      delta2)
    !nrj = sum(phi(this%geomx%nx/2,1:this%geomx%ny,1:this%geomv%nx)**2)*this%geomx%dy*this%geomv%dx


    if(sll_get_collective_rank(sll_world_collective)==0) then
      write(file_id,'(f12.5,2g20.12)') &
        real(step,f64)*dt, &
        nrj

      if(step==0)then    
        call sll_gnuplot_write(sim%phi3d_parx3(:,1,1),'phi_0',ierr)
        !call sll_gnuplot_write(sim%rho3d_seqx1x2(:,1,1)/sim%n0_r(:)-1._f64,'rho_0',ierr)
        call sll_gnuplot_write(sim%Ti_r(:),'Ti_r',ierr)
        call sll_gnuplot_write(sim%Te_r(:),'Te_r',ierr)
        call sll_gnuplot_write(sim%n0_r(:),'n0_r',ierr)
      endif


    endif
    
    
  end subroutine time_history_diagnostic_dk_polar



  subroutine compute_field_from_phi_polar(phi,mesh1,mesh2,A1,A2,interp2d)
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: A1
    sll_real64, dimension(:,:), intent(out) :: A2
    type(sll_logical_mesh_1d), pointer :: mesh1
    type(sll_logical_mesh_1d), pointer :: mesh2
    class(sll_interpolator_2d_base), pointer   :: interp2d
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: x1_min
    sll_real64 :: x2_min
    sll_real64 :: delta_x1
    sll_real64 :: delta_x2
    sll_real64 :: x1
    sll_real64 :: x2
    sll_int32 :: i1
    sll_int32 :: i2
    
    Nc_x1 = mesh1%num_cells
    Nc_x2 = mesh2%num_cells
    x1_min = mesh1%eta_min
    x2_min = mesh2%eta_min
    delta_x1 = mesh1%delta_eta
    delta_x2 = mesh2%delta_eta

    call interp2d%compute_interpolants(phi)

    do i2=1,Nc_x2+1
      x2=x2_min+real(i2-1,f64)*delta_x2
      do i1=1,Nc_x1+1
        x1=x1_min+real(i1-1,f64)*delta_x1
        A1(i1,i2)=interp2d%interpolate_derivative_eta2(x1,x2)/x1
        A2(i1,i2)=-interp2d%interpolate_derivative_eta1(x1,x2)/x1
      end do
    end do
    
    
    
  end subroutine compute_field_from_phi_polar



  subroutine compute_field_from_phi_cartesian_1d(phi,mesh,A,interp)
    sll_real64, dimension(:), intent(in) :: phi
    sll_real64, dimension(:), intent(out) :: A
    type(sll_logical_mesh_1d), pointer :: mesh
    class(sll_interpolator_1d_base), pointer   :: interp
    sll_int32 :: Nc_x1
    sll_real64 :: x1_min
    sll_real64 :: delta_x1
    sll_real64 :: x1
    sll_int32 :: i1
    
    Nc_x1 = mesh%num_cells
    x1_min = mesh%eta_min
    delta_x1 = mesh%delta_eta

    call interp%compute_interpolants(phi)

    do i1=1,Nc_x1+1
      x1=x1_min+real(i1-1,f64)*delta_x1
      A(i1)=interp%interpolate_derivative_eta1(x1)
    end do
  end subroutine compute_field_from_phi_cartesian_1d

  
  subroutine compute_field_dk_parx1( sim )
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar) :: sim
    sll_int32 :: loc4d_sz_x1
    sll_int32 :: loc4d_sz_x2
    sll_int32 :: loc4d_sz_x3
    sll_int32 :: loc4d_sz_x4
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: i3
    sll_int32 :: nc_x1
    sll_int32 :: nc_x2
    sll_int32 :: nc_x3
    
    
    nc_x1 = sim%m_x1%num_cells
    nc_x2 = sim%m_x2%num_cells
    nc_x3 = sim%m_x3%num_cells
    
    
    call compute_local_sizes_4d( &
      sim%layout4d_parx3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
    
    !print *,'#phi3d_parx3=',maxval(sim%phi3d_parx3),minval(sim%phi3d_parx3), &
    !  sll_get_collective_rank(sll_world_collective)
    
    do i3 = 1,loc4d_sz_x3
      call compute_field_from_phi_polar( &
        sim%phi3d_parx3(1:nc_x1+1,1:nc_x2+1,i3), &
        sim%m_x1, &
        sim%m_x2, &
        sim%A1_parx3(1:nc_x1+1,1:nc_x2+1,i3), &
        sim%A2_parx3(1:nc_x1+1,1:nc_x2+1,i3), &
        sim%phi_interp_x1x2)
    enddo

    call compute_local_sizes_4d( &
      sim%layout4d_parx1, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )

    do i2=1, loc4d_sz_x2
      do i1=1, loc4d_sz_x1
        call compute_field_from_phi_cartesian_1d( &
          sim%phi3d_parx1(i1,i2,1:nc_x3+1), &
          sim%m_x3, &
          sim%A3_parx1(i1,i2,1:nc_x3+1), &
          sim%phi_interp_x3)
        sim%A3_parx1(i1,i2,1:nc_x3+1)=-sim%A3_parx1(i1,i2,1:nc_x3+1)   
      enddo
    enddo
 
!    call compute_oblic_derivative( &
!      tau, &
!      phi, &
!      sim%m_x1, &
!      sim%m_x2, &    
!      sim%Dtau_degree, &
!      D_phi, &
!      sim%adv )   
!    print *,'#max A1=',maxval(sim%A1_parx3),minval(sim%A1_parx3),sll_get_Collective_rank(sll_world_collective)
!    print *,'#max A2=',maxval(sim%A2_parx3),minval(sim%A2_parx3),sll_get_Collective_rank(sll_world_collective)
!    print *,'#max A3=',maxval(sim%A3_parx1),minval(sim%A3_parx1),sll_get_Collective_rank(sll_world_collective)
!    print *,'#max phi3d_parx1=',maxval(sim%phi3d_parx1),minval(sim%phi3d_parx1),sll_get_Collective_rank(sll_world_collective)
!    print *,'#max phi3d_parx3=',maxval(sim%phi3d_parx3),minval(sim%phi3d_parx3),sll_get_Collective_rank(sll_world_collective)

  end subroutine compute_field_dk_parx1




  subroutine compute_rho_dk( sim )
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar) :: sim
    sll_int32 :: loc_sz_x1
    sll_int32 :: loc_sz_x2
    sll_int32 :: loc_sz_x3
    sll_int32 :: loc_sz_x4
    
    
    
    
     call compute_local_sizes_4d( sim%layout4d_parx1, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4 )

    !print *,sim%my_rank,loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4



    call compute_reduction_4d_to_3d_direction4(&
      sim%f4d_parx1, &
      sim%rho3d_parx1, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4, &
      sim%m_x4%delta_eta)

 

    
    

  end subroutine compute_rho_dk





  subroutine advection_x3( sim, dt )
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar) :: sim
    sll_real64,dimension(:), allocatable ::  f1d
    sll_real64, intent(in) :: dt
    sll_int32 :: nc_x1
    sll_int32 :: nc_x2
    sll_int32 :: nc_x3
    sll_int32 :: nc_x4
    sll_int32 :: i1
    sll_int32 :: i2
    !sll_int32 :: i3
    sll_int32 :: i4
    sll_int32 :: ierr
    sll_int32 :: global_indices(4)
    sll_real64 :: alpha
    sll_int32 :: loc_sz_x1
    sll_int32 :: loc_sz_x2
    sll_int32 :: loc_sz_x3
    sll_int32 :: loc_sz_x4
    

    nc_x1 = sim%m_x1%num_cells
    nc_x2 = sim%m_x2%num_cells
    nc_x3 = sim%m_x3%num_cells
    nc_x4 = sim%m_x4%num_cells

   
    
    SLL_ALLOCATE(f1d(nc_x3+1),ierr)  
      
      
    call compute_local_sizes_4d( sim%layout4d_parx1, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4 )

     
    do i2=1,loc_sz_x2
      do i1=1,loc_sz_x1
        do i4=1,loc_sz_x4
          global_indices(1:4) = local_to_global_4D( &
            sim%layout4d_parx1, &
            (/i1, i2, 1, i4/) )
          alpha = sim%m_x4%eta_min+real(global_indices(4)-1,f64)*sim%m_x4%delta_eta
          f1d(1:nc_x3+1)=sim%f4d_parx1(i1,i2,1:nc_x3+1,i4) 
          call sim%adv_x3%advect_1d_constant(&
            alpha, &
            dt, &
            f1d(1:nc_x3+1), &
            f1d(1:nc_x3+1))
            sim%f4d_parx1(i1,i2,1:nc_x3+1,i4)=f1d(1:nc_x3+1)
         enddo
      enddo
    enddo    
    
    SLL_DEALLOCATE_ARRAY(f1d,ierr)
    
  end subroutine advection_x3
  


  subroutine advection_x4( sim, dt )
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar) :: sim
    sll_real64,dimension(:), allocatable ::  f1d
    sll_real64,dimension(:), allocatable ::  f1d_new
    sll_real64, intent(in) :: dt
    sll_int32 :: nc_x1
    sll_int32 :: nc_x2
    sll_int32 :: nc_x3
    sll_int32 :: nc_x4
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: i3
    !sll_int32 :: i4
    sll_int32 :: ierr
    !sll_int32 :: global_indices(4)
    sll_real64 :: alpha
    sll_int32 :: loc_sz_x1
    sll_int32 :: loc_sz_x2
    sll_int32 :: loc_sz_x3
    sll_int32 :: loc_sz_x4
    

    nc_x1 = sim%m_x1%num_cells
    nc_x2 = sim%m_x2%num_cells
    nc_x3 = sim%m_x3%num_cells
    nc_x4 = sim%m_x4%num_cells


    
    SLL_ALLOCATE(f1d(nc_x4+1),ierr)  
    SLL_ALLOCATE(f1d_new(nc_x4+1),ierr)  
      
    call compute_local_sizes_4d( sim%layout4d_parx1, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4 )
    do i3 =1,loc_sz_x3
      do i2=1,loc_sz_x2
        do i1=1,loc_sz_x1
          alpha = sim%A3_parx1(i1,i2,i3)
          f1d(1:nc_x4+1)=sim%f4d_parx1(i1,i2,i3,1:nc_x4+1) 
          call sim%adv_x4%advect_1d_constant(&
            alpha, &
            dt, &
            f1d(1:nc_x4+1), &
            f1d_new(1:nc_x4+1))
            sim%f4d_parx1(i1,i2,i3,1:nc_x4+1)=f1d_new(1:nc_x4+1)
        enddo
      enddo      
    enddo    

    SLL_DEALLOCATE_ARRAY(f1d,ierr)
    SLL_DEALLOCATE_ARRAY(f1d_new,ierr)
    
  end subroutine advection_x4


  subroutine advection_x1x2( sim, dt )
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar) :: sim
    sll_real64,dimension(:,:), allocatable ::  f2d
    sll_real64,dimension(:,:), allocatable ::  f2d_new
    sll_real64, intent(in) :: dt
    sll_int32 :: nc_x1
    sll_int32 :: nc_x2
    sll_int32 :: nc_x3
    sll_int32 :: nc_x4
    !sll_int32 :: i1
    !sll_int32 :: i2
    sll_int32 :: i3
    sll_int32 :: i4
    sll_int32 :: ierr
    !sll_int32 :: global_indices(4)
    !sll_real64 :: alpha
    sll_int32 :: loc_sz_x1
    sll_int32 :: loc_sz_x2
    sll_int32 :: loc_sz_x3
    sll_int32 :: loc_sz_x4
    

    nc_x1 = sim%m_x1%num_cells
    nc_x2 = sim%m_x2%num_cells
    nc_x3 = sim%m_x3%num_cells
    nc_x4 = sim%m_x4%num_cells



    call apply_remap_4D( &
      sim%remap_plan_parx1_to_parx3x4, &
      sim%f4d_parx1, &
      sim%f4d_parx3x4 )



    
    SLL_ALLOCATE(f2d(nc_x1+1,nc_x2+1),ierr)  
    SLL_ALLOCATE(f2d_new(nc_x1+1,nc_x2+1),ierr)  
    !SLL_ALLOCATE(A1(nc_x1+1,nc_x2+1),ierr)  
    !SLL_ALLOCATE(A2(nc_x1+1,nc_x2+1),ierr)  
      
    call compute_local_sizes_4d( sim%layout4d_parx3x4, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4 )
    do i4 = 1,loc_sz_x4
      do i3 = 1,loc_sz_x3
        !A1 = 0._f64
        !A2 = 0._f64
        f2d(1:nc_x1+1,1:nc_x2+1)=sim%f4d_parx3x4(1:nc_x1+1,1:nc_x2+1,i3,i4) 
        call sim%adv_x1x2%advect_2d(&
          sim%A1_parx3(1:nc_x1+1,1:nc_x2+1,i3), &
          sim%A2_parx3(1:nc_x1+1,1:nc_x2+1,i3), &
          dt, &
          f2d(1:nc_x1+1,1:nc_x2+1), &
          f2d_new(1:nc_x1+1,1:nc_x2+1))
        sim%f4d_parx3x4(1:nc_x1+1,1:nc_x2+1,i3,i4)=f2d_new(1:nc_x1+1,1:nc_x2+1)
      enddo  
    enddo    

    SLL_DEALLOCATE_ARRAY(f2d,ierr)
    SLL_DEALLOCATE_ARRAY(f2d_new,ierr)
    !SLL_DEALLOCATE_ARRAY(A1,ierr)
    !SLL_DEALLOCATE_ARRAY(A2,ierr)


    call apply_remap_4D( &
      sim%remap_plan_parx3x4_to_parx1, &
      sim%f4d_parx3x4, &
      sim%f4d_parx1 )

    
  end subroutine advection_x1x2

  
  subroutine delete_dk4d_field_aligned_polar( sim )
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar) :: sim
    !sll_int32 :: ierr
    if(sll_get_collective_rank(sll_world_collective)==0)then    
      print *,'#delete_dk4d_polar not implemented'
      print *,sim%dt
    endif
    
  end subroutine delete_dk4d_field_aligned_polar
  
  
  subroutine initialize_profiles_analytic(sim)
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar), intent(inout) :: sim
    sll_int32 :: i,ierr,nc_x1
    sll_real64 :: x1,delta_x1,rpeak,tmp,x1_min,x1_max
    sll_real64 :: inv_Ln
    sll_real64 :: inv_LTi
    sll_real64 :: inv_LTe
    sll_real64 :: R0
    sll_real64 :: x3_min
    sll_real64 :: x3_max
    sll_real64 :: deltarn
    sll_real64 :: deltarTi
    sll_real64 :: deltarTe
    sll_real64 :: Lr
    sll_real64 :: Lz
    
    nc_x1 = sim%m_x1%num_cells
    delta_x1 = sim%m_x1%delta_eta
    x1_min = sim%m_x1%eta_min
    x1_max = sim%m_x1%eta_max
    x3_min = sim%m_x3%eta_min
    x3_max = sim%m_x3%eta_max
    
    Lr = x1_max-x1_min
    Lz = x3_max-x3_min
    
    R0 = Lz/(2._f64*sll_pi)
    inv_Ln = sim%kappan/R0
    inv_LTi = sim%kappaTi/R0
    inv_LTe = sim%kappaTe/R0
    deltarn = sim%deltarn * Lr
    deltarTi = sim%deltarTi * Lr
    deltarTe = sim%deltarTe * Lr
    
   
    SLL_ALLOCATE(sim%n0_r(nc_x1+1),ierr)
    SLL_ALLOCATE(sim%Ti_r(nc_x1+1),ierr)
    SLL_ALLOCATE(sim%Te_r(nc_x1+1),ierr)
    SLL_ALLOCATE(sim%dlog_density_r(nc_x1+1),ierr)
    SLL_ALLOCATE(sim%B_norm_r(nc_x1+1),ierr)
    SLL_ALLOCATE(sim%Bstar_par_v_r(nc_x1+1),ierr)
    SLL_ALLOCATE(sim%c_r(nc_x1+1),ierr)
    
    rpeak = x1_min+sim%rho_peak*(x1_max-x1_min)
    do i=1,nc_x1+1
      x1 = x1_min+real(i-1,f64)*delta_x1
      !sim%n0_r(i) = exp(-sim%kappan*sim%deltarn*tanh((x1-rpeak)/sim%deltarn))
      !sim%Ti_r(i)=exp(-sim%kappaTi*sim%deltarTi*tanh((x1-rpeak)/sim%deltarTi))    
      !sim%Te_r(i)=exp(-sim%kappaTe*sim%deltarTe*tanh((x1-rpeak)/sim%deltarTe))
      !sim%dlog_density_r(i) = -sim%kappan*cosh((x1-rpeak)/sim%deltarn)**(-2)    

      sim%n0_r(i) = exp(-inv_Ln*deltarn*tanh((x1-rpeak)/deltarn))
      sim%Ti_r(i)=exp(-inv_LTi*deltarTi*tanh((x1-rpeak)/deltarTi))    
      sim%Te_r(i)=exp(-inv_LTe*deltarTe*tanh((x1-rpeak)/deltarTe))
      sim%dlog_density_r(i) = -inv_Ln*cosh((x1-rpeak)/deltarn)**(-2)    
      !constant q case
      sim%c_r(i) = x1*sim%iota_r(i)/R0!x1/(R0*sim%q0)
      sim%B_norm_r(i) = (1._f64+sim%c_r(i)**2)**(sim%B_norm_exponent)
      !sim%Bstar_par_v_r = (2._f64*sim%c_r(i)-R0*sim%Dr_q0*sim%c_r(i)**2)/(1+sim%c_r(i)**2)
      sim%Bstar_par_v_r(i) = (2._f64*sim%c_r(i)-sim%Diota_r(i)*sim%c_r(i)**2)/(1+sim%c_r(i)**2)
      !print *,i,sim%Bstar_par_v_r(i)
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
  subroutine allocate_fdistribu4d_and_QN_DK_parx1( sim )
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar), intent(inout) :: sim

    sll_int32 :: ierr !, itemp
    sll_int32 :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4
    sll_int32 :: power2
    sll_int32 :: power2_x3
    sll_int32 :: k_min
    sll_int32 :: nproc3d_x1_parx1
    


    sim%nproc_x1 = sll_get_collective_size(sll_world_collective)
    sim%nproc_x2 = 1
    sim%nproc_x3 = 1
    sim%nproc_x4 = 1
   
    if(sll_get_collective_rank(sll_world_collective)==0)then
    
      print *,'#num_points=',sim%m_x1%num_cells+1, &
        sim%m_x2%num_cells+1, &
        sim%m_x3%num_cells+1, &
        sim%m_x4%num_cells+1
   
      print *,'#num_proc=',sll_get_collective_size(sll_world_collective)
   
      print *,'#num_proc_parx1:',sim%nproc_x1,sim%nproc_x2,sim%nproc_x3,sim%nproc_x4
    endif
    
    power2 = int(log(real(sll_get_collective_size(sll_world_collective)))/log(2.0))
    power2_x3 = int(log(real(sim%m_x3%num_cells+1))/log(2.0))
    power2_x3 = min(power2_x3,power2)

    !--> Initialization of parallel layout of f4d in x1 direction
    !-->  (x2,x3,x4) : sequential
    !-->  x1 : parallelized layout
    sim%layout4d_parx1  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells+1, & 
      sim%m_x3%num_cells+1, &
      sim%m_x4%num_cells+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%layout4d_parx1 )
    
    ! Allocate the array needed to store the local chunk 
    ! of the distribution function data. First compute the 
    ! local sizes. Since the remap operations
    ! are out-of-place, we will allocate two different arrays, 
    ! one for each layout.
    call compute_local_sizes_4d( sim%layout4d_parx1, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
    
    !print *,'#locsize',sim%my_rank,loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4
    
    
      
      
    SLL_ALLOCATE(sim%f4d_parx1(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)
    SLL_ALLOCATE(sim%phi3d_parx1(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)
    SLL_ALLOCATE(sim%rho3d_parx1(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)
    SLL_ALLOCATE(sim%A3_parx1(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)

    sim%layout2d_parx1  => new_layout_2D( sll_world_collective )
    call initialize_layout_with_distributed_2D_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells, & 
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%layout2d_parx1 )

    sim%layout2d_parx2  => new_layout_2D( sll_world_collective )
    call initialize_layout_with_distributed_2D_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells, & 
      sim%nproc_x2, &
      sim%nproc_x1, &
      sim%layout2d_parx2 )





    !--> Initialization of parallel layout of f4d in (x3,x4) directions
    !-->  (x3,x4) : parallelized layout
    !-->  (x1,x2) : sequential
    !--> we take the most important number of points in x3
    
    sim%nproc_x1 = 1
    sim%nproc_x2 = 1
    sim%nproc_x3 = 2**(power2_x3)
    sim%nproc_x4 = 2**(power2-power2_x3)
     

    sim%layout4d_parx3x4  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells+1, & 
      sim%m_x3%num_cells+1, &
      sim%m_x4%num_cells+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%layout4d_parx3x4 )
        
    call compute_local_sizes_4d( sim%layout4d_parx3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )    
    SLL_ALLOCATE(sim%f4d_parx3x4(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)
    SLL_ALLOCATE(sim%rho3d_parx3(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_parx3(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)
    SLL_ALLOCATE(sim%A1_parx3(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)
    SLL_ALLOCATE(sim%A2_parx3(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)
    
    
    sim%remap_plan_parx1_to_parx3x4 => NEW_REMAP_PLAN( &
      sim%layout4d_parx1, &
      sim%layout4d_parx3x4, &
      sim%f4d_parx1)
    sim%remap_plan_parx3x4_to_parx1 => NEW_REMAP_PLAN( &
      sim%layout4d_parx3x4, &
      sim%layout4d_parx1, &
      sim%f4d_parx3x4)

    if(sll_get_collective_rank(sll_world_collective)==0)then
      print *,'#num_proc_parx3x4:',sim%nproc_x1,sim%nproc_x2,sim%nproc_x3,sim%nproc_x4
    endif


    !col_f => get_layout_collective( layout_f )    
    !my_rank_f  = sll_get_collective_rank( col_f )
    !col_sz_f  = sll_get_collective_size( col_f )
    
    k_min = get_layout_k_min( &
      sim%layout4d_parx3x4, &
      sll_get_collective_rank(sll_world_collective) )
    !k_max = get_layout_k_max( layout_f, my_rank_f )

    !print *,'kmin for new=',k_min,k_max,my_rank_f
    sim%new_collective_per_locx3 => sll_new_collective( &
      sll_world_collective, &
      k_min, &
      sll_get_collective_rank(sll_world_collective))    

    sim%new_collective_per_locx4 => sll_new_collective( &
      sll_world_collective, &
      sll_get_collective_rank(sim%new_collective_per_locx3), &
      sll_get_collective_rank(sll_world_collective))    

    
    nproc3d_x1_parx1 = sll_get_collective_size(sll_world_collective)
    
    sim%layout3d_parx1  => layout_3D_from_layout_4D( sim%layout4d_parx1 )
    sim%layout3d_parx3  => layout_3D_from_layout_4D( sim%layout4d_parx3x4 )
    
    !new_layout_3D( sim%new_collective_per_locx4 )
!    call initialize_layout_with_distributed_3D_array( &
!      sim%m_x1%num_cells+1, & 
!      sim%m_x2%num_cells+1, & 
!      sim%m_x3%num_cells+1, &
!      nproc3d_x1_parx1, &
!      nproc3d_x2_parx1, &
!      nproc3d_x3_parx1, &
!      sim%layout3d_parx1 )

!   if(sll_get_collective_rank(sll_world_collective).eq. 0) then
!     do i=0,sll_get_collective_size(sll_world_collective)-1
!       print *,'layout_parx1', i, &
!         get_layout_3D_i_min(sim%layout3d_parx1,i), &
!         get_layout_3D_i_max(sim%layout3d_parx1,i), &
!         get_layout_3D_j_min(sim%layout3d_parx1,i), &
!         get_layout_3D_j_max(sim%layout3d_parx1,i), &
!         get_layout_3D_k_min(sim%layout3d_parx1,i), & 
!         get_layout_3D_k_max(sim%layout3d_parx1,i)
!     enddo
!   endif
!
!   if(sll_get_collective_rank(sll_world_collective).eq. 0) then
!     do i=0,sll_get_collective_size(sll_world_collective)-1
!       print *,'layout_parx3', i, &
!         get_layout_3D_i_min(sim%layout3d_parx3,i), &
!         get_layout_3D_i_max(sim%layout3d_parx3,i), &
!         get_layout_3D_j_min(sim%layout3d_parx3,i), &
!         get_layout_3D_j_max(sim%layout3d_parx3,i), &
!         get_layout_3D_k_min(sim%layout3d_parx3,i), & 
!         get_layout_3D_k_max(sim%layout3d_parx3,i)
!     enddo
!   endif


    sim%remap_plan_parx1_to_parx3 => NEW_REMAP_PLAN( &
      sim%layout3d_parx1, &
      sim%layout3d_parx3, &
      sim%phi3d_parx1)

    
  end subroutine allocate_fdistribu4d_and_QN_DK_parx1









  !----------------------------------------------------
  ! Initialization of the distribution function for
  !   drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine initialize_fdistribu4d_DK(sim,layout,f4d)
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar), intent(inout) :: sim
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
    
    call init_fequilibrium( &
      sim%m_x1%num_cells+1, &
      sim%m_x4%num_cells+1, &
      x1_node, &
      x4_node, &
      sim%n0_r, &
      sim%Ti_r, &
      sim%feq_x1x4 )
!    do i1=1,sim%m_x1%num_cells+1
!      do i4=1,sim%m_x4%num_cells+1 
!        sim%feq_x1x4(i1,i4) = compute_equil_analytic(sim,x1_node(i1),x4_node(i4))  
!      enddo
!    enddo 
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
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar), intent(in) :: sim
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
 
  subroutine solve_quasi_neutral_parx1(sim)
    class(sll_simulation_4d_drift_kinetic_field_aligned_polar), intent(inout) :: sim
      sll_int32 :: glob_ind(4)
      sll_int32 :: loc4d_sz_x1
      sll_int32 :: loc4d_sz_x2
      sll_int32 :: loc4d_sz_x3
      sll_int32 :: loc4d_sz_x4
      sll_int32 :: i
      sll_real64, dimension(:), allocatable :: tmp
      sll_int32 :: ierr

      call compute_local_sizes_4d( &
        sim%layout4d_parx1, &
        loc4d_sz_x1, &
        loc4d_sz_x2, &
        loc4d_sz_x3, &
        loc4d_sz_x4 )        
      
      SLL_ALLOCATE(tmp(loc4d_sz_x1),ierr)
      
      do i=1,loc4d_sz_x1
        glob_ind(:) = local_to_global_4D( &
          sim%layout4d_parx1, &
          (/i,1,1,1/))
        tmp(i) = 1._f64/sim%n0_r(glob_ind(1))  
      enddo
      do i=1,loc4d_sz_x1
        sim%rho3d_parx1(i,:,:) =  sim%rho3d_parx1(i,:,:)*tmp(i)-1._f64
      enddo
      call sim%poisson3d%compute_phi_from_rho( &
        sim%phi3d_parx1(:,1:loc4d_sz_x2-1,:), &
        sim%rho3d_parx1(:,1:loc4d_sz_x2-1,:) )
      
      sim%phi3d_parx1(:,loc4d_sz_x2,:) = sim%phi3d_parx1(:,1,:)  
      !print *,'#maxval rho=',maxval(sim%rho3d_parx1)
      !print *,'#minval rho=',minval(sim%rho3d_parx1)
      
      !print *,'#maxval phi=',maxval(sim%phi3d_parx1)
      !print *,'#minval phi=',minval(sim%phi3d_parx1)
      !sim%phi3d_parx1 = sll_get_collective_rank(sll_world_collective)
       
      call apply_remap_3D( &
          sim%remap_plan_parx1_to_parx3, &
          sim%phi3d_parx1, &
          sim%phi3d_parx3 ) 
          
      !print *,'#hello!'     
!      print *,'phi3d_parx1',maxval(sim%phi3d_parx1(:,1:loc4d_sz_x2-1,:)),minval(sim%phi3d_parx1(:,1:loc4d_sz_x2-1,:)),&
!        sll_get_collective_rank(sll_world_collective)
!      print *,'phi3d_parx3',maxval(sim%phi3d_parx3(:,1:loc4d_sz_x2-1,:)),minval(sim%phi3d_parx3(:,1:loc4d_sz_x2-1,:)),&
!        sll_get_collective_rank(sll_world_collective)
!       print *,'phi3d_parx1 Nt',maxval(sim%phi3d_parx1(:,loc4d_sz_x2,:)),minval(sim%phi3d_parx1(:,loc4d_sz_x2,:)),&
!        sll_get_collective_rank(sll_world_collective)
!      print *,'phi3d_parx3 Nt',maxval(sim%phi3d_parx3(:,loc4d_sz_x2,:)),minval(sim%phi3d_parx3(:,loc4d_sz_x2,:)),&
!        sll_get_collective_rank(sll_world_collective)
!     
!      call sll_halt_collective()
!      stop
        

      SLL_DEALLOCATE_ARRAY(tmp,ierr)
    
    
  end subroutine solve_quasi_neutral_parx1

 
 
  

#ifndef NOHDF5
!*********************
!*********************

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f_polar(iplot,f,m_x1,m_x2)
    use sll_xdmf
    use sll_hdf5_io_serial
    sll_int32 :: file_id
    sll_int32 :: error
    sll_real64, dimension(:,:), allocatable :: x1
    sll_real64, dimension(:,:), allocatable :: x2
    sll_int32 :: i, j
    sll_int32, intent(in) :: iplot
    character(len=4)      :: cplot
    sll_int32             :: nnodes_x1, nnodes_x2
    type(sll_logical_mesh_1d), pointer :: m_x1
    type(sll_logical_mesh_1d), pointer :: m_x2
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64 :: r
    sll_real64 :: theta
    sll_real64 :: rmin
    sll_real64 :: rmax
    sll_real64 :: dr
    sll_real64 :: dtheta
    
    
    nnodes_x1 = m_x1%num_cells+1
    nnodes_x2 = m_x2%num_cells+1
    rmin = m_x1%eta_min
    rmax = m_x1%eta_max
    dr = m_x1%delta_eta
    dtheta = m_x2%delta_eta
    
    !print *,'#maxf=',iplot,maxval(f),minval(f)
    

    
    if (iplot == 1) then

      SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
      SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
      do j = 1,nnodes_x2
        do i = 1,nnodes_x1
          r       = rmin+real(i-1,f32)*dr
          theta   = real(j-1,f32)*dtheta
          x1(i,j) = r*cos(theta)
          x2(i,j) = r*sin(theta)
        end do
      end do
      call sll_hdf5_file_create("polar_mesh-x1.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_hdf5_file_close(file_id, error)
      call sll_hdf5_file_create("polar_mesh-x2.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_hdf5_file_close(file_id, error)
      deallocate(x1)
      deallocate(x2)

    end if

    call int2string(iplot,cplot)
    call sll_xdmf_open("f_x1x2_"//cplot//".xmf","polar_mesh", &
      nnodes_x1,nnodes_x2,file_id,error)
    call sll_xdmf_write_array("f_x1x2_"//cplot,f,"values", &
      error,file_id,"Node")
    call sll_xdmf_close(file_id,error)
  end subroutine plot_f_polar

#endif

#ifndef NOHDF5
!*********************
!*********************

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f_cartesian(iplot,f,m_x1,m_x2)
    use sll_xdmf
    use sll_hdf5_io_serial
    sll_int32 :: file_id
    sll_int32 :: error
    sll_real64, dimension(:,:), allocatable :: x1
    sll_real64, dimension(:,:), allocatable :: x2
    sll_int32 :: i, j
    sll_int32, intent(in) :: iplot
    character(len=4)      :: cplot
    sll_int32             :: nnodes_x1, nnodes_x2
    type(sll_logical_mesh_1d), pointer :: m_x1
    type(sll_logical_mesh_1d), pointer :: m_x2
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64 :: r
    sll_real64 :: theta
    sll_real64 :: rmin
    sll_real64 :: rmax
    sll_real64 :: dr
    sll_real64 :: dtheta
    
    
    nnodes_x1 = m_x1%num_cells+1
    nnodes_x2 = m_x2%num_cells+1
    rmin = m_x1%eta_min
    rmax = m_x1%eta_max
    dr = m_x1%delta_eta
    dtheta = m_x2%delta_eta
    
    !print *,'#maxf=',iplot,maxval(f),minval(f)
    

    
    if (iplot == 1) then

      SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
      SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
      do j = 1,nnodes_x2
        do i = 1,nnodes_x1
          r       = rmin+real(i-1,f32)*dr
          theta   = real(j-1,f32)*dtheta
          x1(i,j) = r!*cos(theta)
          x2(i,j) = theta!*sin(theta)
        end do
      end do
      call sll_hdf5_file_create("cartesian_mesh-x1.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_hdf5_file_close(file_id, error)
      call sll_hdf5_file_create("cartesian_mesh-x2.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_hdf5_file_close(file_id, error)
      deallocate(x1)
      deallocate(x2)

    end if

    call int2string(iplot,cplot)
    call sll_xdmf_open("f_x2x3_"//cplot//".xmf","cartesian_mesh", &
      nnodes_x1,nnodes_x2,file_id,error)
    call sll_xdmf_write_array("f_x2x3_"//cplot,f,"values", &
      error,file_id,"Node")
    call sll_xdmf_close(file_id,error)
  end subroutine plot_f_cartesian

#endif



  

end module sll_simulation_4d_drift_kinetic_field_aligned_polar_module





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!here we give some routines for the computations of the fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
!!
!!  !>  compute second component of magnetic field
!  subroutine compute_B_x2( &
!    B_x2_array, &
!    x1_array, &
!    num_points_x1, &
!    R, &
!    a )
!    sll_real64, dimension(:), intent(out) :: B_x2_array
!    sll_real64, dimension(:), intent(in) :: x1_array
!    sll_int32, intent(in) :: num_points_x1
!    sll_real64, intent(in) :: R
!    sll_real64, intent(in) :: a
!    sll_real64 :: q
!    sll_int32 :: i
!    sll_real64 :: x1
!    
!    if(size(B_x2_array,1)<num_points_x1)then
!      print *,'#bad size of B_x2_array in compute_B_x2', &
!        size(B_x2_array,1), &
!        num_points_x1
!      stop
!    endif
!
!    if(size(x1_array,1)<num_points_x1)then
!      print *,'#bad size of x1_array in compute_B_x2', &
!        size(B_x2_array,1), &
!        num_points_x1
!      stop
!    endif
!    do i=1,num_points_x1    
!      x1 = x1_array(i)
!      q = 1+(x1/a)**2
!      B_x2_array(i) = 1._f64+(x1/(R*q))**2
!    enddo  
!  end subroutine compute_B_x2
!!
!!  !>  compute the norm of the magnetic field
!!  
!  subroutine compute_B_norm( &
!    B_norm_array, &
!    B_x2_array, &
!    num_points_x1)
!    sll_real64, dimension(:), intent(in) :: B_x2_array
!    sll_real64, dimension(:), intent(out) :: B_norm_array
!    sll_int32, intent(in) :: num_points_x1
!    sll_int32 :: i
!
!    if(size(B_x2_array,1)<num_points_x1)then
!      print *,'#bad size of B_x2_array in compute_B_norm', &
!        size(B_x2_array,1), &
!        num_points_x1
!      stop
!    endif
!    if(size(B_norm_array,1)<num_points_x1)then
!      print *,'#bad size of B_norm_array in compute_B_norm', &
!        size(B_norm_array,1), &
!        num_points_x1
!      stop
!    endif
!    do i=1,num_points_x1    
!      B_norm_array(i) = sqrt(1._f64+B_x2_array(i)**2)
!    enddo          
!    
!  end subroutine compute_B_norm
!!
!!  !>  compute the unit magnetic field
!!
!  subroutine compute_b_unit( &
!    b_unit_x2_array, &
!    b_unit_x3_array, &
!    B_x2_array, &
!    num_points_x1)
!    sll_real64, dimension(:), intent(in) :: B_x2_array
!    sll_real64, dimension(:), intent(out) :: b_unit_x2_array
!    sll_real64, dimension(:), intent(out) :: b_unit_x3_array
!    sll_int32, intent(in) :: num_points_x1
!    sll_int32 :: i
!
!    if(size(B_x2_array,1)<num_points_x1)then
!      print *,'#bad size of B_x2_array in compute_B_norm', &
!        size(B_x2_array,1), &
!        num_points_x1
!      stop
!    endif
!    if(size(b_unit_x2_array,1)<num_points_x1)then
!      print *,'#bad size of b_unit_x2_array in compute_b_unit', &
!        size(b_unit_x2_array,1), &
!        num_points_x1
!      stop
!    endif
!    if(size(b_unit_x3_array,1)<num_points_x1)then
!      print *,'#bad size of b_unit_x3_array in compute_b_unit', &
!        size(b_unit_x2_array,1), &
!        num_points_x1
!      stop
!    endif
!    do i=1,num_points_x1    
!      b_unit_x2_array(i) = B_x2_array(i)/sqrt(1._f64+B_x2_array(i)**2)
!      b_unit_x3_array(i) = 1._f64/sqrt(1._f64+B_x2_array(i)**2)
!    enddo          
!    
!  end subroutine compute_b_unit
!!> orthogonal coordinate system
!!> reference: http://en.wikipedia.org/wiki/Orthogonal_coordinates
!!> reference: http://en.wikipedia.org/wiki/Curvilinear_coordinates
!!> eta_1,eta_2,eta_3
!!> stand for q_1,q_2,q_3 in the reference
!!> x stands for r in the reference


!!> eta = (eta_1,eta_2,eta_3)
!!> \mathbf{x} = (x_1,x_2,x_3) = x_1e_1+x_2e_2+x_3e_3 cartesian basis
!!> \mathbf{e_1} = \partial_{x_1} \mathbf{x}
!!> \mathbf{e_2} = \partial_{x_2} \mathbf{x}
!!> \mathbf{e_3} = \partial_{x_3} \mathbf{x}
!!> the transformation is given by the change form the cartesian grid
!!> \mathbf{x} = \mathbf{x}(\mathbf{eta})
!!> that is
!!> x_1(eta_1,eta_2,eta_3) = ...
!!> x_2(eta_1,eta_2,eta_3) = ...
!!> x_3(eta_1,eta_2,eta_3) = ...
!!> we define
!!> \mathbf{h_1} = \partial_{eta_1} \mathbf{x}
!!> \mathbf{h_2} = \partial_{eta_2} \mathbf{x}
!!> \mathbf{h_3} = \partial_{eta_3} \mathbf{x}
!!> and (covariant normalized basis = contravariant normalized basis in orthogonal geometry)
!!> \hat{h_1} = \mathbf{h_1}/h_1, h_1 = |\mathbf{h_1}|
!!> \hat{h_2} = \mathbf{h_2}/h_2, h_2 = |\mathbf{h_2}|
!!> \hat{h_3} = \mathbf{h_3}/h_3, h_3 = |\mathbf{h_3}|

!!> gradient of a scalar field in orthogonal coordinate system
!!> \nabla\phi = (\hat{h_1}/h_1) \partial_{eta_1}\phi
!!>   +(\hat{h_2}/h_2) \partial_{eta_2}\phi
!!>   +(\hat{h_3}/h_3) \partial_{eta_3}\phi

!!> Vector field \mathbf{F}
!!> F_1 = \mathbf{F} \cdot \hat{h_1}
!!> F_2 = \mathbf{F} \cdot \hat{h_2}
!!> F_3 = \mathbf{F} \cdot \hat{h_3}

!!> divergence of a vector field in orthogonal coordinate system
!!> \nabla\cdot\mathbf{F} = 1/(h1h2h3) [ \partial_{eta_1}(F_1h_2h_3)
!!>   +\partial_{eta_2}(F_2h_3h_1)  
!!>   +\partial_{eta_3}(F_3h_1h_2) ]  

!!> curl of a vector field in orthogonal coordinate system
!!> \nabla \times \mathbf{F} = 1/(h1h2h3) 
!!  Det[ h_1\hat{h_1} h_2\hat{h_2} h_3\hat{h_3}
!!>   \partial_{eta_1} \partial_{eta_2} \partial_{eta_3}  
!!>   h_1F_1 h_2F_2 h_3F_3 ]  

!!> Laplacian of a scalar field in orthogonal coordinate system
!!> \nabla^2 \phi = 1/(h1h2h3) 
!!  Det[ h_1\hat{h_1} h_2\hat{h_2} h_3\hat{h_3}
!!>   \partial_{eta_1} \partial_{eta_2} \partial_{eta_3}  
!!>   h_1F_1 h_2F_2 h_3F_3 ]  


!!> specify domain for eta_1,eta_2,eta_3 renamed
!!> specify change from cartesian
!!> specify scale factors h_1,h_2,h_3




!!> example of cylindrical coordinates
!!> eta = (r,theta,z)
!!> x_1 = r*cos(theta)
!!> x_2 = r*sin(theta)
!!> x_3 = z
!!> h_1 = 1
!!> h_2 = r
!!> h_3 = 1

!!> \nabla \phi(r,theta,z) = (\partial_r \phi)\hat{r}
!!> +(\partial_theta \phi)/r\hat{theta} +\partial_z \phi\hat{z}

!!> alpha = iota /R, R=L/(2pi), L = z_max-z_min
!!> We suppose that the magnetic field writes
!!> B = B_norm b, b=b_theta hat_theta + b_z hat_z
!!> hat z = b/b_z - (b_theta/b_z) hat_theta
!!> as example, we have
!!> B_norm = (1+alpha^2*r^2)**(-1/2)
!!> b_theta = alpha*r/(1+alpha^2*r^2)**(1/2)
!!> b_z = 1/(1+alpha^2*r^2)**(1/2)
!!> b_theta/b_z = alpha*r
!!> b_theta^2+b_z^2 = 1
!!> Db_theta = alpha/(1+alpha^2*r^2)^(3/2)
!!> Db_z = -b_theta*Db_theta/b_z
!!> Db_z = -alpha^2*r/(1+alpha^2*r^2)^(3/2)
!!> curl_b_theta = -rDb_z = (alpha*r)^2/(1+alpha^2*r^2)^(3/2)
!!>   = r*(b_theta/b_z)*Db_theta
!!> curl_b_z = D(r*b_theta) = b_theta+r*Db_theta = alpha*r*(2+alpha^2*r^2)/(1+alpha ^2*r^2)^(3/2)
!!> curl_b_dot_b = r*(b_theta/b_z)*Db_theta*b_theta+(b_theta+r*Db_theta)*b_z
!!>   = b_z*(b_theta+rDb_theta*(1+(b_theta/b_z)**2)
!!>   = 2*alpha*r/(1+alpha^2*r^2)
!!> Bstar_par = B_norm+v*curl_b_dot_b
!!>   = (1+alpha^2*r^2)**(-1/2) + v*(2*alpha*r)/(1+alpha^2*r^2)
!! grad_phi_r = Dr_phi
!! grad_phi_theta = Dtheta_phi/r
!! grad_phi_z = Dz_phi

!!> Bstar_theta = B_norm*b_theta+v*curl_b_theta
!!>   = alpha*r/(1+alpha^2*r^2)+v*(alpha*r)^2/(1+alpha^2*r^2)^(3/2)
!!> Bstar_z = B_norm*b_z+v*curl_b_z
!!>   = 1/(1+alpha^2*r^2)+v*alpha*r*(2+alpha^2*r^2)/(1+alpha ^2*r^2)^(3/2)

!!> Bstar = Bstar_theta hat_theta+Bstar_z hat_z
!!>   = (Bstar_theta-(b_theta/b_z)Bstar_z)hat_theta+(Bstar_z/b_z)b
!!>   = v*(curl_b_theta-(b_theta/b_z)*curl_b_z)hat_theta+(Bstar_z/b_z)b


!!> bstar_dot_grad_phi = (B_norm*b_theta+v*curl_b_theta)*Dtheta_phi/r
!!>  +(B_norm*b_z+v*curl_b_z)*Dz_phi

!!> b_dot_grad_phi = b_theta*Dtheta_phi/r+b_z*Dz_phi

!!> bstar_dot_grad_phi = B_norm*b_dot_grad_phi+v*(TRUC)
!!> TRUC = curl_b_theta*Dtheta_phi/r+curl_b_z*Dz_phi
!!>   = curl_b_theta*Dtheta_phi/r+curl_b_z*(b_dot_grad_phi-b_theta*Dtheta_phi/r)/b_z
!!>   = (curl_b_z/b_z)*b_dot_grad_phi+(curl_b_theta-curl_b_z*b_theta/b_z)*Dtheta_phi/r
!!>   = alpha*r*(1+1/(1+alpha^2*r^2))*b_dot_grad_phi-(alpha^2*r^2)/(1+alpha^2*r^2)^(1/2)*Dtheta_phi/r

!
!  subroutine compute_curl_b_unit_cubic_splines( &
!    curl_b_unit_x2_array, &
!    curl_b_unit_x3_array, &
!    b_unit_x2_array, &
!    b_unit_x3_array, &
!    x1_array, &
!    num_points_x1)
!    sll_real64, dimension(:), intent(in) :: b_unit_x2_array
!    sll_real64, dimension(:), intent(in) :: b_unit_x3_array
!    sll_real64, dimension(:), intent(out) :: curl_b_unit_x2_array
!    sll_real64, dimension(:), intent(out) :: curl_b_unit_x3_array
!    sll_real64, dimension(:), intent(in) :: x1_array
!    sll_int32, intent(in) :: num_points_x1
!    sll_int32 :: i
!    
!    
!    
!  end subroutine compute_curl_b_unit_cubic_splines
!
!
!
!
!
!  !>  compute B star parallel in an array
!  subroutine compute_B_star_parallel( &
!    B_star_parallel_array, &
!    B_x2_array, &
!    x1_array, &
!    num_points_x1, &
!    R, &
!    a ) &
!    result(res)
!    sll_real64, dimension(:), intent(in) :: B_x2_array
!    sll_real64, dimension(:), intent(in) :: x1_array
!    sll_int32, intent(in) :: num_points_x1
!    sll_real64, intent(in) :: R
!    sll_real64, intent(in) :: a
!    sll_real64 :: q
!    sll_int32 :: i
!    sll_real64 :: x1
!    
!    if(size(B_x2_array,1)<num_points_x1)then
!      print *,'#bad size of B_x2_array in compute_B', &
!        size(B_x2_array,1), &
!        num_points_x1
!      print *,'#in subroutine compute_B'
!      stop
!    endif
!    do i=1,num_points_x1    
!      x1 = x1_array(i)
!      q = 1+(x1/a)**2
!      B_x2_array(i) = 1._f64+(x1/(R*q))**2
!    enddo  
!  end subroutine compute_B_star_parallel
!
!
!
!



