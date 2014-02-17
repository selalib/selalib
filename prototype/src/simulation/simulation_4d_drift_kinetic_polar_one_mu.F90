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
!> Christophe Steiner
!> Michel Mehrenberger (mehrenbe@math.unistra.fr)
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


module sll_simulation_4d_drift_kinetic_polar_one_mu_module
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
  use sll_module_gyroaverage_2d_polar_hermite_solver
  use sll_module_gyroaverage_2d_polar_splines_solver
  use sll_module_gyroaverage_2d_polar_pade_solver


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




  
  type, extends(sll_simulation_base_class) :: &
    sll_simulation_4d_drift_kinetic_polar_one_mu

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
     sll_int32  :: freq_diag_time
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
     !--> Gyroaverage
     sll_real64  :: mu 

     !--> 4D logical mesh (r,theta,phi,vpar)
     !type(sll_logical_mesh_4d), pointer :: logical_mesh4d


     !--> Density and temperature profiles
     sll_real64, dimension(:)  , pointer :: n0_r
     sll_real64, dimension(:)  , pointer :: Ti_r
     sll_real64, dimension(:)  , pointer :: Te_r
     sll_real64, dimension(:)  , pointer :: dlog_density_r

     !--> Equilibrium distribution function
     sll_real64, dimension(:,:), pointer :: feq_x1x4


     !--> 4D distribution function 
     !----> sequential in (x1,x2,x4) and parallel in (x3)
     type(layout_4D), pointer :: layout4d_seqx1x2x4
     sll_real64, dimension(:,:,:,:), pointer :: f4d_seqx1x2x4 
     !----> parallel in (x3) and sequential in (x1,x2,x4) 
     type(layout_4D), pointer :: layout4d_seqx3
     sll_real64, dimension(:,:,:,:), pointer :: f4d_seqx3
     !----> definition of remap
     type(remap_plan_4D_real64), pointer ::remap_plan_seqx1x2x4_to_seqx3
     type(remap_plan_4D_real64), pointer ::remap_plan_seqx3_to_seqx1x2x4
     

     !--> 3D charge density and 3D electric potential
     !----> sequential in (x1,x2)
     type(layout_3D), pointer :: layout3d_seqx1x2
     sll_real64, dimension(:,:,:), pointer :: rho3d_seqx1x2 
     sll_real64, dimension(:,:,:), pointer :: phi3d_seqx1x2 
     sll_real64, dimension(:,:,:), pointer :: A1_seqx1x2 
     sll_real64, dimension(:,:,:), pointer :: A2_seqx1x2 
     sll_real64, dimension(:,:,:), pointer :: A3_seqx1x2 
     !----> sequential in x3
     type(layout_3D), pointer :: layout3d_seqx3
     sll_real64, dimension(:,:,:), pointer :: rho3d_seqx3
     sll_real64, dimension(:,:,:), pointer :: phi3d_seqx3
     sll_real64, dimension(:,:,:), pointer :: A3_seqx3
     !----> definition of remap
     type(remap_plan_3D_real64), pointer ::remap_plan_seqx1x2_to_seqx3
     type(remap_plan_3D_real64), pointer ::remap_plan_seqx3_to_seqx1x2

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
    
    class(sll_gyroaverage_2d_base), pointer :: gyroaverage

    class(sll_poisson_2d_base), pointer   :: poisson2d
    class(sll_poisson_2d_base), pointer   :: poisson2d_mean


    !for computing advection field from phi
    class(sll_interpolator_2d_base), pointer   :: phi_interp_x1x2
    class(sll_interpolator_1d_base), pointer   :: phi_interp_x3


     !--> temporary structures that are used in CG_polar
     !type(sll_SL_polar), pointer :: plan_sl_polar


   contains
     procedure, pass(sim) :: run => run_dk4d_polar
     procedure, pass(sim) :: init_from_file => init_dk4d_polar
  end type sll_simulation_4d_drift_kinetic_polar_one_mu

  interface delete
     module procedure delete_dk4d_polar
  end interface delete

contains

!we should not give directly the file here
!but a long list of parameters that would be initialized with
!a read_from_file routine

  subroutine init_dk4d_polar( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_drift_kinetic_polar_one_mu), intent(inout) :: sim
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
    sll_real64 :: eta_min_gyro(2)
    sll_real64 :: eta_max_gyro(2)
    sll_int32 :: Nc_gyro(2)
    sll_int32 :: interp_degree_gyro(2)


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
    !sll_int32  :: QN_case
    !--> Pertubation
    sll_int32  :: perturb_choice
    sll_int32  :: mmode
    sll_int32  :: nmode
    sll_real64 :: eps_perturb   
    !--> Gyroaverage
    sll_real64              :: mu
    character(len=256)      :: gyroaverage_case
    sll_int32               :: gyroaverage_N_points
    sll_int32               :: gyroaverage_interp_degree_x1
    sll_int32               :: gyroaverage_interp_degree_x2
    !--> Algorithm
    sll_real64 :: dt
    sll_int32  :: number_iterations
    sll_int32  :: freq_diag_time
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
      poisson2d_BC_rmax
    namelist /perturbation/ &
      perturb_choice, &
      mmode, &
      nmode, &
      eps_perturb
    namelist /sim_params/ &
      dt, & 
      number_iterations, &
      freq_diag_time, &
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
      poisson2d_case, &
      gyroaverage_case, &
      mu, &
      gyroaverage_N_points, &
      gyroaverage_interp_degree_x1, &
      gyroaverage_interp_degree_x2
      
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
    
    SLL_ALLOCATE(tmp_r(num_cells_x1+1,2),ierr)
    

    
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
      print *,'##gyroaverage'
      print *,'#gyroaverage_case=',gyroaverage_case
      print *,'#mu=',mu
      print *,'#gyroaverage_N_points=',gyroaverage_N_points
      print *,'#gyroaverage_interp_degree=',gyroaverage_interp_degree_x1,gyroaverage_interp_degree_x2   
    endif
    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)

    
    call initialize_profiles_analytic(sim)    
    call allocate_fdistribu4d_DK(sim)
    call allocate_QN_DK( sim )
    
    
    call initialize_eta1_node_1d(sim%m_x1,sim%x1_node)
    call initialize_eta1_node_1d(sim%m_x2,sim%x2_node)
    call initialize_eta1_node_1d(sim%m_x3,sim%x3_node)
    call initialize_eta1_node_1d(sim%m_x4,sim%x4_node)
    
    
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

          
      case default
        print *,'#bad poisson2d_case',poisson2d_case
        print *,'#not implemented'
        print *,'#in init_dk4d_polar'
        stop
    end select
    
    !--> gyroaverage
    
    sim%mu = mu
    
    eta_min_gyro(1) = sim%m_x1%eta_min
    eta_max_gyro(1) = sim%m_x1%eta_max
    eta_min_gyro(2) = 0._f64
    eta_max_gyro(2) = 2._f64*sll_pi
    
    print *,"eta_min_gyro=",eta_min_gyro
    print *,"eta_max_gyro=",eta_max_gyro
    
    Nc_gyro(1)=sim%m_x1%num_cells
    Nc_gyro(2)=sim%m_x2%num_cells
    interp_degree_gyro(1)=gyroaverage_interp_degree_x1
    interp_degree_gyro(2)=gyroaverage_interp_degree_x2
  
    select case (gyroaverage_case)
      case ("HERMITE")       

        sim%gyroaverage => new_gyroaverage_2d_polar_hermite_solver( &
          eta_min_gyro, &
          eta_max_gyro, &
          Nc_gyro, &
          gyroaverage_N_points, &
          interp_degree_gyro, &
          1)
          
      case ("HERMITE_C1")       

        sim%gyroaverage => new_gyroaverage_2d_polar_hermite_solver( &
          eta_min_gyro, &
          eta_max_gyro, &
          Nc_gyro, &
          gyroaverage_N_points, &
          interp_degree_gyro, &
          2)
          
       case ("SPLINES")       

        sim%gyroaverage => new_gyroaverage_2d_polar_splines_solver( &
          eta_min_gyro, &
          eta_max_gyro, &
          Nc_gyro, &
          gyroaverage_N_points, &
          1)
          
       case ("PADE")       

        sim%gyroaverage => new_gyroaverage_2d_polar_pade_solver( &
          eta_min_gyro, &
          eta_max_gyro, &
          Nc_gyro)
          
      case default
        print *,'#bad gyroaverage_case',gyroaverage_case
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
        


     
    

  end subroutine init_dk4d_polar

  subroutine run_dk4d_polar(sim)
    class(sll_simulation_4d_drift_kinetic_polar_one_mu), intent(inout) :: sim
    !--> For initial profile HDF5 saving
    integer                      :: file_err
    sll_int32                    :: file_id
    character(len=12), parameter :: filename_prof = "init_prof.h5"
    sll_real64,dimension(:,:,:,:), allocatable :: f4d_store
    sll_int32 :: loc4d_sz_x1
    sll_int32 :: loc4d_sz_x2
    sll_int32 :: loc4d_sz_x3
    sll_int32 :: loc4d_sz_x4
    sll_int32 :: iter
    sll_int32 :: nc_x1
    sll_int32 :: nc_x2
    sll_int32 :: nc_x3
    sll_int32 :: nc_x4
    !sll_int32 :: i1
    !sll_int32 :: i2
    !sll_int32 :: i3
    !sll_int32 :: i4
    sll_int32 :: ierr
    sll_real64 :: dt
    sll_int32 :: th_diag_id 
    
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

      
    end if


    call compute_local_sizes_4d( sim%layout4d_seqx1x2x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
    SLL_ALLOCATE(f4d_store(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)

    if(sll_get_collective_rank(sll_world_collective)==0) then
      call sll_ascii_file_create('thdiag.dat', th_diag_id, ierr)
    endif

    call initialize_fdistribu4d_DK(sim,sim%layout4d_seqx1x2x4,sim%f4d_seqx1x2x4)

        
    do iter=1,sim%num_iterations    


      call compute_local_sizes_4d( sim%layout4d_seqx1x2x4, &
        loc4d_sz_x1, &
        loc4d_sz_x2, &
        loc4d_sz_x3, &
        loc4d_sz_x4 )
    
    !print *,'#locsize it',iter,sim%my_rank,loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4

   
      
      call compute_rho_dk(sim)    

      if(sll_get_collective_rank(sll_world_collective)==0) then
        print*,'#iteration=',iter
      
        if(iter==1)then
          call sll_gnuplot_write( &
            sim%rho3d_seqx1x2(1:nc_x1+1,1,1)/sim%n0_r(1:nc_x1+1)-1._f64, &
            'rho_0_init', &
            ierr)
        endif
        if(iter==2)then
          call sll_gnuplot_write( &
            sim%rho3d_seqx1x2(1:nc_x1+1,1,1)/sim%n0_r(1:nc_x1+1)-1._f64, &
            'rho_1_init', &
            ierr)
        endif
        if(iter==2)then
          call sll_gnuplot_write( &
            sim%phi3d_seqx1x2(1:nc_x1+1,1,1), &
            'phi_1', &
            ierr)
        endif
      endif

          
      !call solve_quasi_neutral( sim )
      call solve_quasi_neutral_with_gyroaverage( sim )
      call compute_field_dk( sim )
      call gyroaverage_field_dk( sim )


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
          f4d_store = sim%f4d_seqx1x2x4
          call advection_x3( sim, 0.5_f64*dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x1x2( sim, 0.5_f64*dt )
          call compute_rho_dk(sim)    
          call solve_quasi_neutral( sim )
          call compute_field_dk( sim )
          
          !correction
          sim%f4d_seqx1x2x4 = f4d_store
          call advection_x3( sim, 0.5_f64*dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x1x2( sim, dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x3( sim, 0.5_f64*dt )
 
        case default
          print *,'#sim%time_case=',sim%time_case
          print *, '#not implemented'
          print *,'#in run_dk4d_polar'
          stop
      end select          
    
    enddo
    

  end subroutine run_dk4d_polar

  subroutine time_history_diagnostic_dk_polar( &
    sim, &
    file_id, &    
    step)
    class(sll_simulation_4d_drift_kinetic_polar_one_mu) :: sim
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
      sim%phi3d_seqx1x2(:,:,1)**2, &
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
        call sll_gnuplot_write(sim%phi3d_seqx1x2(:,1,1),'phi_0',ierr)
        call sll_gnuplot_write(sim%rho3d_seqx1x2(:,1,1)/sim%n0_r(:)-1._f64,'rho_0',ierr)
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

  
  
  subroutine compute_field_dk( sim )
    class(sll_simulation_4d_drift_kinetic_polar_one_mu) :: sim
    sll_int32 :: loc_sz_x1
    sll_int32 :: loc_sz_x2
    sll_int32 :: loc_sz_x3
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: i3
    sll_int32 :: nc_x1
    sll_int32 :: nc_x2
    sll_int32 :: nc_x3
    
    
    nc_x1 = sim%m_x1%num_cells
    nc_x2 = sim%m_x2%num_cells
    nc_x3 = sim%m_x3%num_cells
    
    
    call compute_local_sizes_3d( &
      sim%layout3d_seqx1x2, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3 )
    
    do i3 = 1,loc_sz_x3
      call compute_field_from_phi_polar( &
        sim%phi3d_seqx1x2(1:nc_x1+1,1:nc_x2+1,i3), &
        sim%m_x1, &
        sim%m_x2, &
        sim%A1_seqx1x2(1:nc_x1+1,1:nc_x2+1,i3), &
        sim%A2_seqx1x2(1:nc_x1+1,1:nc_x2+1,i3), &
        sim%phi_interp_x1x2)
    enddo

    call compute_local_sizes_3d( &
      sim%layout3d_seqx3, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3 )

    do i2=1, loc_sz_x2
      do i1=1, loc_sz_x1
        call compute_field_from_phi_cartesian_1d( &
          sim%phi3d_seqx3(i1,i2,1:nc_x3+1), &
          sim%m_x3, &
          sim%A3_seqx3(i1,i2,1:nc_x3+1), &
          sim%phi_interp_x3)
        sim%A3_seqx3(i1,i2,1:nc_x3+1)=-sim%A3_seqx3(i1,i2,1:nc_x3+1)   
      enddo
    enddo

    call apply_remap_3D( &
      sim%remap_plan_seqx3_to_seqx1x2, &
      sim%A3_seqx3, &
      sim%A3_seqx1x2 )  




  end subroutine compute_field_dk


  subroutine compute_rho_dk( sim )
    class(sll_simulation_4d_drift_kinetic_polar_one_mu) :: sim
    sll_int32 :: loc_sz_x1
    sll_int32 :: loc_sz_x2
    sll_int32 :: loc_sz_x3
    sll_int32 :: loc_sz_x4
    
    
    
    
     call compute_local_sizes_4d( sim%layout4d_seqx1x2x4, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4 )

    !print *,sim%my_rank,loc_sz_x1,loc_sz_x2,loc_sz_x3,loc_sz_x4



    call compute_reduction_4d_to_3d_direction4(&
      sim%f4d_seqx1x2x4, &
      sim%rho3d_seqx1x2, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4, &
      sim%m_x4%delta_eta)

 
 
    call apply_remap_3D( &
      sim%remap_plan_seqx1x2_to_seqx3, &
      sim%rho3d_seqx1x2, &
      sim%rho3d_seqx3 )

    
    

  end subroutine compute_rho_dk





  subroutine advection_x3( sim, dt )
    class(sll_simulation_4d_drift_kinetic_polar_one_mu) :: sim
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
      
      
      
    call apply_remap_4D( &
      sim%remap_plan_seqx1x2x4_to_seqx3, &
      sim%f4d_seqx1x2x4, &
      sim%f4d_seqx3 )
    call compute_local_sizes_4d( sim%layout4d_seqx3, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4 )
     
    do i2=1,loc_sz_x2
      do i1=1,loc_sz_x1
        do i4=1,loc_sz_x4
          global_indices(1:4) = local_to_global_4D( &
            sim%layout4d_seqx3, &
            (/i1, i2, 1, i4/) )
          alpha = sim%m_x4%eta_min+real(global_indices(4)-1,f64)*sim%m_x4%delta_eta
          f1d(1:nc_x3+1)=sim%f4d_seqx3(i1,i2,1:nc_x3+1,i4) 
          call sim%adv_x3%advect_1d_constant(&
            alpha, &
            dt, &
            f1d(1:nc_x3+1), &
            f1d(1:nc_x3+1))
            sim%f4d_seqx3(i1,i2,1:nc_x3+1,i4)=f1d(1:nc_x3+1)
         enddo
      enddo
    enddo    
    call apply_remap_4D( &
        sim%remap_plan_seqx3_to_seqx1x2x4, &
        sim%f4d_seqx3, &
        sim%f4d_seqx1x2x4 )
    
    SLL_DEALLOCATE_ARRAY(f1d,ierr)
    
  end subroutine advection_x3
  


  subroutine advection_x4( sim, dt )
    class(sll_simulation_4d_drift_kinetic_polar_one_mu) :: sim
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
      
    call compute_local_sizes_4d( sim%layout4d_seqx1x2x4, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4 )
    do i3 =1,loc_sz_x3
      do i2=1,loc_sz_x2
        do i1=1,loc_sz_x1
          alpha = sim%A3_seqx1x2(i1,i2,i3)
          f1d(1:nc_x4+1)=sim%f4d_seqx1x2x4(i1,i2,i3,1:nc_x4+1) 
          call sim%adv_x4%advect_1d_constant(&
            alpha, &
            dt, &
            f1d(1:nc_x4+1), &
            f1d_new(1:nc_x4+1))
            sim%f4d_seqx1x2x4(i1,i2,i3,1:nc_x4+1)=f1d_new(1:nc_x4+1)
        enddo
      enddo      
    enddo    

    SLL_DEALLOCATE_ARRAY(f1d,ierr)
    SLL_DEALLOCATE_ARRAY(f1d_new,ierr)
    
  end subroutine advection_x4


  subroutine advection_x1x2( sim, dt )
    class(sll_simulation_4d_drift_kinetic_polar_one_mu) :: sim
    sll_real64,dimension(:,:), allocatable ::  f2d
    sll_real64,dimension(:,:), allocatable ::  f2d_new
    sll_real64,dimension(:,:), allocatable ::  A1
    sll_real64,dimension(:,:), allocatable ::  A2
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


    
    SLL_ALLOCATE(f2d(nc_x1+1,nc_x2+1),ierr)  
    SLL_ALLOCATE(f2d_new(nc_x1+1,nc_x2+1),ierr)  
    SLL_ALLOCATE(A1(nc_x1+1,nc_x2+1),ierr)  
    SLL_ALLOCATE(A2(nc_x1+1,nc_x2+1),ierr)  
      
    call compute_local_sizes_4d( sim%layout4d_seqx1x2x4, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4 )
    do i4 = 1, nc_x4+1
      do i3 =1,loc_sz_x3
        A1 = 0._f64
        A2 = 0._f64
        f2d(1:nc_x1+1,1:nc_x2+1)=sim%f4d_seqx1x2x4(1:nc_x1+1,1:nc_x2+1,i3,i4) 
        call sim%adv_x1x2%advect_2d(&
          sim%A1_seqx1x2(1:nc_x1+1,1:nc_x2+1,i3), &
          sim%A2_seqx1x2(1:nc_x1+1,1:nc_x2+1,i3), &
          dt, &
          f2d(1:nc_x1+1,1:nc_x2+1), &
          f2d_new(1:nc_x1+1,1:nc_x2+1))
        sim%f4d_seqx1x2x4(1:nc_x1+1,1:nc_x2+1,i3,i4)=f2d_new(1:nc_x1+1,1:nc_x2+1)
      enddo  
    enddo    

    SLL_DEALLOCATE_ARRAY(f2d,ierr)
    SLL_DEALLOCATE_ARRAY(A1,ierr)
    SLL_DEALLOCATE_ARRAY(A2,ierr)
    
  end subroutine advection_x1x2

  
  
  
  
  

  subroutine delete_dk4d_polar( sim )
    class(sll_simulation_4d_drift_kinetic_polar_one_mu) :: sim
    !sll_int32 :: ierr
    
    print *,'#delete_dk4d_polar not implemented'
    print *,sim%dt
    
  end subroutine delete_dk4d_polar
  
  
  subroutine initialize_profiles_analytic(sim)
    class(sll_simulation_4d_drift_kinetic_polar_one_mu), intent(inout) :: sim
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
    class(sll_simulation_4d_drift_kinetic_polar_one_mu), intent(inout) :: sim

    sll_int32 :: ierr !, itemp
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
    sim%layout4d_seqx1x2x4  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells+1, & 
      sim%m_x3%num_cells+1, &
      sim%m_x4%num_cells+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%layout4d_seqx1x2x4 )
    
    ! Allocate the array needed to store the local chunk 
    ! of the distribution function data. First compute the 
    ! local sizes. Since the remap operations
    ! are out-of-place, we will allocate two different arrays, 
    ! one for each layout.
    call compute_local_sizes_4d( sim%layout4d_seqx1x2x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
    
    !print *,'#locsize',sim%my_rank,loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4
    
    
      
      
      
    SLL_ALLOCATE(sim%f4d_seqx1x2x4(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)

    !--> Initialization of parallel layout of f4d in (x1,x2,x4) directions
    !-->  (x1,x2,x4) : parallelized layout
    !-->  (x3) : sequential
    
    sim%nproc_x1 = 2**(sim%power2/3)
    sim%nproc_x2 = 2**(sim%power2/3)
    sim%nproc_x3 = 1
    sim%nproc_x4 = 2**(sim%power2-2*(sim%power2/3))
     

    sim%layout4d_seqx3  => new_layout_4D( sll_world_collective )
    call initialize_layout_with_distributed_4D_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells+1, & 
      sim%m_x3%num_cells+1, &
      sim%m_x4%num_cells+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%layout4d_seqx3 )
        
    call compute_local_sizes_4d( sim%layout4d_seqx3, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )    
    SLL_ALLOCATE(sim%f4d_seqx3(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)
    
    
    sim%remap_plan_seqx1x2x4_to_seqx3 => NEW_REMAP_PLAN( &
      sim%layout4d_seqx1x2x4, &
      sim%layout4d_seqx3, &
      sim%f4d_seqx1x2x4)
    sim%remap_plan_seqx3_to_seqx1x2x4 => NEW_REMAP_PLAN( &
      sim%layout4d_seqx3, &
      sim%layout4d_seqx1x2x4, &
      sim%f4d_seqx3)

    
    
  end subroutine allocate_fdistribu4d_DK


  !----------------------------------------------------
  ! Initialization of the distribution function for
  !   drift-kinetic 4D simulation
  !----------------------------------------------------
  subroutine initialize_fdistribu4d_DK(sim,layout,f4d)
    class(sll_simulation_4d_drift_kinetic_polar_one_mu), intent(inout) :: sim
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
    class(sll_simulation_4d_drift_kinetic_polar_one_mu), intent(in) :: sim
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
    class(sll_simulation_4d_drift_kinetic_polar_one_mu), intent(inout) :: sim

    !type(sll_logical_mesh_2d), pointer :: logical_mesh2d
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
    sim%layout3d_seqx1x2  => new_layout_3D( sll_world_collective )
    nproc3d_x3 = sim%nproc_x3*sim%nproc_x4
    call initialize_layout_with_distributed_3D_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells+1, & 
      sim%m_x3%num_cells+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      nproc3d_x3, &
      sim%layout3d_seqx1x2 )
    call compute_local_sizes_3d( &
      sim%layout3d_seqx1x2, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    SLL_ALLOCATE(sim%rho3d_seqx1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_seqx1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%A1_seqx1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%A2_seqx1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%A3_seqx1x2(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
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
        
    sim%layout3d_seqx3  => new_layout_3D( sll_world_collective )
    call initialize_layout_with_distributed_3D_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells+1, & 
      sim%m_x3%num_cells+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%layout3d_seqx3 )
    call compute_local_sizes_3d( &
      sim%layout3d_seqx3, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    SLL_ALLOCATE(sim%rho3d_seqx3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_seqx3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    SLL_ALLOCATE(sim%A3_seqx3(loc3d_sz_x1,loc3d_sz_x2,loc3d_sz_x3),ierr)
    
    
    sim%remap_plan_seqx1x2_to_seqx3 => NEW_REMAP_PLAN( &
      sim%layout3d_seqx1x2, &
      sim%layout3d_seqx3, &
      sim%rho3d_seqx1x2)
    sim%remap_plan_seqx3_to_seqx1x2 => NEW_REMAP_PLAN( &
      sim%layout3d_seqx3, &
      sim%layout3d_seqx1x2, &
      sim%rho3d_seqx3)

    
   end subroutine allocate_QN_DK

 
 
  
  subroutine solve_quasi_neutral(sim)
    class(sll_simulation_4d_drift_kinetic_polar_one_mu), intent(inout) :: sim
    sll_int32 :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3
    sll_int32 :: iloc1, iloc2, iloc3
    !sll_int32 :: i1, i2
    sll_int32 :: i3
    sll_real64 :: tmp
    sll_int32 :: glob_ind(3)    



    
    select case (sim%QN_case)
      case (SLL_NO_QUASI_NEUTRAL)
      ! no quasi neutral solver as in CRPP-CONF-2001-069
        call compute_local_sizes_3d( &
          sim%layout3d_seqx3, &
          loc3d_sz_x1, &
          loc3d_sz_x2, &
          loc3d_sz_x3 )        
        if((loc3d_sz_x3).ne.(sim%m_x3%num_cells+1))then
          print *,'#Problem of parallelization dimension in solve_quasi_neutral'
          print *,'#sll_simulation_4d_drift_kinetic_polar type simulation'
          stop
        endif        
        do iloc2 = 1, loc3d_sz_x2
          do iloc1 = 1, loc3d_sz_x1          
            tmp = sum(sim%rho3d_seqx3(iloc1,iloc2,1:sim%m_x3%num_cells))&
              /real(sim%m_x3%num_cells,f64)
            SLL_ASSERT(loc3d_sz_x3==sim%m_x3%num_cells+1)
            do i3 = 1,sim%m_x3%num_cells+1
              glob_ind(:) = local_to_global_3D(sim%layout3d_seqx3, &
                (/iloc1,iloc2,i3/))                        
              sim%phi3d_seqx3(iloc1,iloc2,i3) = (sim%rho3d_seqx3(iloc1,iloc2,i3)-tmp)&
                *sim%Te_r(glob_ind(1))/sim%n0_r(glob_ind(1))
            enddo    
          enddo
        enddo  
        call apply_remap_3D( &
          sim%remap_plan_seqx3_to_seqx1x2, &
          sim%phi3d_seqx3, &
          sim%phi3d_seqx1x2 )  
      case (SLL_QUASI_NEUTRAL_WITHOUT_ZONAL_FLOW)
        call compute_local_sizes_3d( &
          sim%layout3d_seqx1x2, &
          loc3d_sz_x1, &
          loc3d_sz_x2, &
          loc3d_sz_x3 )
                  
        do iloc2=1, loc3d_sz_x2
          do iloc1=1, loc3d_sz_x1
            sim%phi3d_seqx1x2(iloc1,iloc2,:) = &
              sim%rho3d_seqx1x2(iloc1,iloc2,:)/sim%n0_r(iloc1)-1._f64
          enddo
        enddo
        do iloc3=1, loc3d_sz_x3
          call sim%poisson2d%compute_phi_from_rho( &
            sim%phi3d_seqx1x2(:,:,iloc3), &
            sim%phi3d_seqx1x2(:,:,iloc3) )
        enddo
        call apply_remap_3D( &
          sim%remap_plan_seqx1x2_to_seqx3, &
          sim%phi3d_seqx1x2, &
          sim%phi3d_seqx3 )            
      case (SLL_QUASI_NEUTRAL_WITH_ZONAL_FLOW)
        print *,'#SLL_QUASI_NEUTRAL_WITH_ZONAL_FLOW'
        print *,'#not implemented yet '
        stop      
      case default
        print *,'#bad value for sim%QN_case'
        stop  
    end select        
  
  end subroutine solve_quasi_neutral
  
  
  
  
  subroutine solve_quasi_neutral_with_gyroaverage(sim)
    class(sll_simulation_4d_drift_kinetic_polar_one_mu), intent(inout) :: sim
    sll_int32 :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3
    sll_int32 :: iloc1, iloc2, iloc3
    !sll_int32 :: i1, i2
    sll_int32 :: i3
    sll_real64 :: tmp
    sll_int32 :: glob_ind(3)    
    sll_int32 :: nc_x1, nc_x2
            
    nc_x1 = sim%m_x1%num_cells
    nc_x2 = sim%m_x2%num_cells

    select case (sim%QN_case)
      case (SLL_NO_QUASI_NEUTRAL)
      ! no quasi neutral solver as in CRPP-CONF-2001-069
        call compute_local_sizes_3d( &
          sim%layout3d_seqx3, &
          loc3d_sz_x1, &
          loc3d_sz_x2, &
          loc3d_sz_x3 )        
        if((loc3d_sz_x3).ne.(sim%m_x3%num_cells+1))then
          print *,'#Problem of parallelization dimension in solve_quasi_neutral'
          print *,'#sll_simulation_4d_drift_kinetic_polar type simulation'
          stop
        endif        
        do iloc2 = 1, loc3d_sz_x2
          do iloc1 = 1, loc3d_sz_x1          
            tmp = sum(sim%rho3d_seqx3(iloc1,iloc2,1:sim%m_x3%num_cells))&
              /real(sim%m_x3%num_cells,f64)
            SLL_ASSERT(loc3d_sz_x3==sim%m_x3%num_cells+1)
            do i3 = 1,sim%m_x3%num_cells+1
              glob_ind(:) = local_to_global_3D(sim%layout3d_seqx3, &
                (/iloc1,iloc2,i3/))                        
              sim%phi3d_seqx3(iloc1,iloc2,i3) = (sim%rho3d_seqx3(iloc1,iloc2,i3)-tmp)&
                *sim%Te_r(glob_ind(1))/sim%n0_r(glob_ind(1))
            enddo    
          enddo
        enddo  
        call apply_remap_3D( &
          sim%remap_plan_seqx3_to_seqx1x2, &
          sim%phi3d_seqx3, &
          sim%phi3d_seqx1x2 )  
      case (SLL_QUASI_NEUTRAL_WITHOUT_ZONAL_FLOW)
        call compute_local_sizes_3d( &
          sim%layout3d_seqx1x2, &
          loc3d_sz_x1, &
          loc3d_sz_x2, &
          loc3d_sz_x3 )
            
        do iloc3 = 1,loc3d_sz_x3    
          call sim%gyroaverage%compute_gyroaverage( &
          sim%mu, &
          sim%rho3d_seqx1x2(1:nc_x1+1,1:nc_x2+1,iloc3))   
        enddo    
                  
        do iloc2=1, loc3d_sz_x2
          do iloc1=1, loc3d_sz_x1
            sim%phi3d_seqx1x2(iloc1,iloc2,:) = &
              sim%rho3d_seqx1x2(iloc1,iloc2,:)/sim%n0_r(iloc1)-1._f64
          enddo
        enddo
        do iloc3=1, loc3d_sz_x3
          call sim%poisson2d%compute_phi_from_rho( &
            sim%phi3d_seqx1x2(:,:,iloc3), &
            sim%phi3d_seqx1x2(:,:,iloc3) )
        enddo
        call apply_remap_3D( &
          sim%remap_plan_seqx1x2_to_seqx3, &
          sim%phi3d_seqx1x2, &
          sim%phi3d_seqx3 )            
      case (SLL_QUASI_NEUTRAL_WITH_ZONAL_FLOW)
        print *,'#SLL_QUASI_NEUTRAL_WITH_ZONAL_FLOW'
        print *,'#not implemented yet '
        stop      
      case default
        print *,'#bad value for sim%QN_case'
        stop  
    end select        
  
  end subroutine solve_quasi_neutral_with_gyroaverage
  



subroutine gyroaverage_field_dk(sim)
    class(sll_simulation_4d_drift_kinetic_polar_one_mu), intent(inout) :: sim
    sll_int32 :: loc_sz_x1
    sll_int32 :: loc_sz_x2
    sll_int32 :: loc_sz_x3
    sll_int32 :: i3
    sll_int32 :: nc_x1, nc_x2
            
    nc_x1 = sim%m_x1%num_cells
    nc_x2 = sim%m_x2%num_cells

    call compute_local_sizes_3d( &
      sim%layout3d_seqx1x2, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3 )

    do i3 = 1,loc_sz_x3
      call sim%gyroaverage%compute_gyroaverage( &
        sqrt(2*sim%mu), &
        sim%A1_seqx1x2(1:nc_x1+1,1:nc_x2+1,i3))
      call sim%gyroaverage%compute_gyroaverage( &
        sqrt(2*sim%mu), &
        sim%A2_seqx1x2(1:nc_x1+1,1:nc_x2+1,i3))
      call sim%gyroaverage%compute_gyroaverage( &
        sqrt(2*sim%mu), &
        sim%A3_seqx1x2(1:nc_x1+1,1:nc_x2+1,i3))
    enddo          


  end subroutine gyroaverage_field_dk

  

end module sll_simulation_4d_drift_kinetic_polar_one_mu_module



