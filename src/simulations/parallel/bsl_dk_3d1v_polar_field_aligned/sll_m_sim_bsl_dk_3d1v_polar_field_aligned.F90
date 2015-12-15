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
!> Yaman Güçlü (yaman.guclu@gmail.com)
!> @brief 
!> Simulation class to solve drift kinetic equation in polar coordinates
!> (3d space (x1=r,x2=theta,x3=z) 1d velocity (x4=v))
!> it is an extension of the simulation_4d_drift_kinetic_polar.F90
!> we explore the use of field aligned interpolation
!> using the strategy introduced by Ottaviani
!> we follow here exactly the direction of the magnetic field
!> the model is simpler that the screw pinch test case
!> as a consequence, splitting is still valid
!> we have still constant advection in x4 and along the direction of the magnetic field
!> the magnetic field is monitored by iota
!> which can be given as an array
!> for the moment we suppose constant or linear iota of the form
!> iota = iota0+Diota0*r
!> where   iota0 = 0.8 Dr_iota0 = 0. are examples
!> and can be changed in the namelist file
!> @details
!> Example of use in test program: see unit_test_4d_dk_polar.F90 file
!> 
!> \code
!>
!>  use sll_simulation_4d_drift_kinetic_polar
!>  type(sll_simulation_4d_vp_polar)    :: simulation
!>  call simulation%init_from_file(trim(filename))
!>  call simulation%run()
!>  call sll_o_delete(simulation)
!> \endcode

module sll_m_sim_bsl_dk_3d1v_polar_field_aligned

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_advection_1d_base, only: &
    sll_c_advection_1d_base

  use sll_m_advection_1d_periodic, only: &
    sll_f_new_periodic_1d_advector

  use sll_m_advection_2d_base, only: &
    sll_c_advection_2d_base

  use sll_m_advection_2d_bsl, only: &
    sll_f_new_bsl_2d_advector

  use sll_m_advection_2d_oblic, only: &
    sll_f_new_oblic_2d_advector, &
    sll_t_oblic_2d_advector, &
    sll_s_oblic_advect_2d_constant

  use sll_m_ascii_io, only: &
    sll_s_ascii_file_create

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_hermite, &
    sll_p_neumann, &
    sll_p_neumann_mode_0, &
    sll_p_periodic, &
    sll_p_set_to_limit

  use sll_m_cartesian_meshes, only: &
    sll_o_get_node_positions, &
    sll_f_new_cartesian_mesh_1d, &
    sll_t_cartesian_mesh_1d

  use sll_m_characteristics_2d_base, only: &
    sll_c_characteristics_2d_base

  use sll_m_characteristics_2d_explicit_euler, only: &
    sll_f_new_explicit_euler_2d_charac

  use sll_m_characteristics_2d_verlet, only: &
    sll_f_new_verlet_2d_charac

  use sll_m_collective, only: &
    sll_s_collective_barrier, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_cubic_spline_interpolator_1d, only: &
    sll_f_new_cubic_spline_interpolator_1d

  use sll_m_cubic_spline_interpolator_2d, only: &
    sll_f_new_cubic_spline_interpolator_2d

  use sll_m_derivative_2d_oblic, only: &
    sll_s_compute_oblic_derivative_2d, &
    sll_f_new_oblic_2d_derivative, &
    sll_t_oblic_2d_derivative

  use sll_m_fdistribu4d_dk, only: &
    sll_s_init_fequilibrium

  use sll_m_gnuplot, only: &
    sll_o_gnuplot_1d

  use sll_m_hdf5_io_serial, only: &
    sll_o_hdf5_file_close, &
    sll_o_hdf5_file_create, &
    sll_o_hdf5_write_array

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  use sll_m_periodic_interp, only: &
    sll_p_lagrange, &
    sll_p_spline

  use sll_m_poisson_2d_base, only: &
    sll_c_poisson_2d_base

  use sll_m_poisson_3d_base, only: &
    sll_c_poisson_3d_base

  use sll_m_qn_solver_3d_polar_parallel_x1_wrapper, only: &
    sll_f_new_qn_solver_3d_polar_parallel_x1_wrapper

  use sll_m_reduction, only: &
    sll_s_compute_reduction_2d_to_0d, &
    sll_s_compute_reduction_4d_to_3d_direction4

  use sll_m_remapper, only: &
    sll_o_apply_remap_3d, &
    sll_o_apply_remap_4d, &
    sll_o_compute_local_sizes, &
    sll_o_global_to_local, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_t_layout_2d, &
    sll_t_layout_3d, &
    sll_t_layout_4d, &
    sll_o_local_to_global, &
    sll_f_new_layout_2d, &
    sll_f_new_layout_3d_from_layout_4d, &
    sll_f_new_layout_4d, &
    sll_o_new_remap_plan, &
    sll_t_remap_plan_3d_real64, &
    sll_t_remap_plan_4d_real64

  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  use sll_m_utilities, only: &
    sll_s_int2string

  use sll_xdmf_io, only: &
    sll_t_hdf5_serial

  use sll_xdmf_io_parallel, only: &
    sll_t_xdmf_parallel_file

  implicit none

  public :: &
    sll_o_delete, &
    sll_t_simulation_4d_drift_kinetic_field_aligned_polar

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !! choice of QNS solver
  !! TODO: should be elsewhere
  sll_int32, parameter :: SLL_NO_QUASI_NEUTRAL = 0
  sll_int32, parameter :: SLL_QUASI_NEUTRAL_WITH_ZONAL_FLOW = 1
  sll_int32, parameter :: SLL_QUASI_NEUTRAL_WITHOUT_ZONAL_FLOW = 2

  !! choice of time scheme solver
  !! TODO: should be elsewhere
  sll_int32, parameter :: SLL_TIME_LOOP_EULER = 0
  sll_int32, parameter :: SLL_TIME_LOOP_PREDICTOR_CORRECTOR = 1
  sll_int32, parameter :: SLL_TIME_LOOP_PREDICTOR2_CORRECTOR = 2

  !============================================================================
  ! SIMULATION TYPE: Drift-kinetic, 3D-1V, field-aligned, cylindrical coords.
  !============================================================================
  type, extends(sll_c_simulation_base_class) :: &
    sll_t_simulation_4d_drift_kinetic_field_aligned_polar

    ! Parallel environment parameters
    sll_int32  :: world_size
    sll_int32  :: my_rank
    sll_int32  :: power2 ! 2^power2 = number of processes available
    ! Processor mesh sizes
    sll_int32  :: nproc_x1
    sll_int32  :: nproc_x2
    sll_int32  :: nproc_x3
    sll_int32  :: nproc_x4
!    type(sll_t_collective_t), pointer :: new_collective_per_locx3
!    type(sll_t_collective_t), pointer :: new_collective_per_locx4

    ! Mesh parameters
    type(sll_t_cartesian_mesh_1d), pointer :: m_x1
    type(sll_t_cartesian_mesh_1d), pointer :: m_x2
    type(sll_t_cartesian_mesh_1d), pointer :: m_x3
    type(sll_t_cartesian_mesh_1d), pointer :: m_x4
    
    ! Instantaneous time value (updated at the end of the time-step)
    sll_real64 :: time

    ! Physics/numerical parameters
    sll_real64 :: dt
    sll_int32  :: num_iterations
    sll_int32  :: freq_diag_time
    sll_int32  :: freq_diag
    sll_int32  :: time_case
    sll_int32  :: charac_case
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

    !--> Density and temperature profiles
    sll_real64, pointer :: n0_r(:)
    sll_real64, pointer :: Ti_r(:)
    sll_real64, pointer :: Te_r(:)
    sll_real64, pointer :: dlog_density_r(:)

    !--> Magnetic field
    !sll_real64 :: q0
    !sll_real64 :: Dr_q0
    sll_real64, pointer :: iota_r(:)
    !sll_real64, pointer :: Diota_r(:)
    !sll_real64 :: B_norm_exponent
    !sll_real64, pointer :: B_norm_r(:)
    !sll_real64, pointer :: Bstar_par_v_r(:)
    sll_real64 :: B0
    sll_real64, pointer :: b_unit_x2(:)
    sll_real64, pointer :: b_unit_x3(:)

    !variables that permit to compare the code when non field alignement is used
    !the convenience, her is tu use the same code
    !it maybe temporary
    logical :: use_field_aligned_derivative
    logical :: use_field_aligned_interpolation

    !--> Equilibrium distribution function, 1D-1V
    sll_real64, pointer :: feq_x1x4(:,:)

    !----> parallel in x1
    type(sll_t_layout_4d), pointer :: layout4d_parx1
    sll_real64,      pointer ::      f4d_parx1(:,:,:,:)

    !----> parallel in (x3,x4)
    type(sll_t_layout_4d), pointer :: layout4d_parx3x4
    sll_real64,      pointer ::      f4d_parx3x4(:,:,:,:)

    !----> definition of remap
    type(sll_t_remap_plan_4d_real64), pointer :: remap_plan_parx1_to_parx3x4
    type(sll_t_remap_plan_4d_real64), pointer :: remap_plan_parx3x4_to_parx1

    !----> parallel in x1
    type(sll_t_layout_3d), pointer ::       layout3d_parx1
    sll_real64,      pointer ::          rho3d_parx1(:,:,:)
    sll_real64,      pointer ::          phi3d_parx1(:,:,:)
    sll_real64,      pointer :: Daligned_phi3d_parx1(:,:,:)
    sll_real64,      pointer ::             A3_parx1(:,:,:)

    !----> parallel in x3
    type(sll_t_layout_3d), pointer :: layout3d_parx3
    sll_real64,      pointer ::    rho3d_parx3(:,:,:)
    sll_real64,      pointer ::    phi3d_parx3(:,:,:)
    sll_real64,      pointer ::       A1_parx3(:,:,:)
    sll_real64,      pointer ::       A2_parx3(:,:,:)
    type(sll_t_remap_plan_3d_real64), pointer ::remap_plan_parx1_to_parx3

    !----> for Poisson
    type(sll_t_layout_2d), pointer :: layout2d_parx1
    type(sll_t_layout_2d), pointer :: layout2d_parx2

    !--> cubic splines interpolation
    !type(sll_t_cubic_spline_2d), pointer :: interp_x1x2
    !type(sll_t_cubic_spline_1d), pointer :: interp_x3
    !type(sll_t_cubic_spline_1d), pointer :: interp_x4

    sll_real64, pointer :: x1_node(:)
    sll_real64, pointer :: x2_node(:)
    sll_real64, pointer :: x3_node(:)
    sll_real64, pointer :: x4_node(:)

    class(sll_c_advection_2d_base), pointer :: adv_x1x2
    !class(sll_c_interpolator_2d), pointer :: interp_x1x2
    class(sll_c_characteristics_2d_base), pointer :: charac_x1x2
    class(sll_c_advection_1d_base), pointer :: adv_x3
    class(sll_c_advection_1d_base), pointer :: adv_x4
    type(sll_t_oblic_2d_advector),      pointer :: adv_x2x3
    
    class(sll_c_poisson_2d_base), pointer :: poisson2d
    class(sll_c_poisson_2d_base), pointer :: poisson2d_mean
    class(sll_c_poisson_3d_base), pointer :: poisson3d

    !for computing advection field from phi
    class(sll_c_interpolator_2d), pointer :: phi_interp_x1x2
    class(sll_c_interpolator_1d), pointer :: phi_interp_x3
    class(sll_c_advection_1d_base),    pointer :: adv_x2
    type(sll_t_oblic_2d_derivative),       pointer :: deriv
    !class(sll_c_interpolator_1d), pointer   :: phi_interp_fa !for field aligned interp.
    !should replace phi_interp_x3 in future

   contains

     procedure :: run            =>  run_dk4d_field_aligned_polar
     procedure :: init_from_file => init_dk4d_field_aligned_polar

  end type sll_t_simulation_4d_drift_kinetic_field_aligned_polar

  !============================================================================
  ! sll_o_delete SUBROUTINE (overloaded)
  !============================================================================
  interface sll_o_delete
     module procedure delete_dk4d_field_aligned_polar
  end interface sll_o_delete

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  ! Initialize simulation from input file
  subroutine init_dk4d_field_aligned_polar( sim, filename )
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar), &
                      intent(inout) :: sim
    character(len=*), intent(in)    :: filename

    intrinsic :: trim
    character( len=256 ) :: err_msg ! used by SLL_ERROR and SLL_WARNING

    sll_int32            :: IO_stat
    sll_int32, parameter :: input_file = 99  ! TODO: use Pierre's I/O
    class(sll_c_characteristics_2d_base), pointer :: charac2d   ! computation of characteristics

    !> 2D interpolator (in poloidal plane) for r component of adv. field
    class(sll_c_interpolator_2d), pointer   :: A1_interp2d 

    !> 2D interpolator (in poloidal plane) for theta component of adv. field
    class(sll_c_interpolator_2d), pointer   :: A2_interp2d

    !> 2D interpolator (in poloidal plane) for distribution function
    class(sll_c_interpolator_2d), pointer   :: f_interp2d

    !> 1D interpolators (along r) for (r,theta) components of adv. field
    class(sll_c_interpolator_1d), pointer   :: A1_interp1d_x1
    class(sll_c_interpolator_1d), pointer   :: A2_interp1d_x1

    sll_real64 :: charac2d_tol     !< Tolerance for fixed point iteration
    sll_int32  :: charac2d_maxiter !< Max no. of fixed point iterations

    !--> Mesh
    sll_int32  :: num_cells_x1  !< r
    sll_int32  :: num_cells_x2  !< theta
    sll_int32  :: num_cells_x3  !< z
    sll_int32  :: num_cells_x4  !< v
    sll_real64 :: r_min
    sll_real64 :: r_max
    sll_real64 :: z_min
    sll_real64 :: z_max
    sll_real64 :: v_min
    sll_real64 :: v_max

    !--> Equilibrium
    sll_real64  :: tau0
    sll_real64  :: rho_peak
    sll_real64  :: kappan
    sll_real64  :: deltarn
    sll_real64  :: kappaTi
    sll_real64  :: deltarTi
    sll_real64  :: kappaTe
    sll_real64  :: deltarTe
    sll_real64  :: iota0
    sll_real64  :: Dr_iota0
    character(len=256) :: iota_file
    sll_int32  :: size_iota_file
    logical    :: is_iota_file
    sll_real64 :: B0

    ! TODO: check if following quantities can be used in the future
    !sll_int32 :: size_Diota_file
    !logical :: is_Diota_file
    !character(len=256) :: Diota_file
    !sll_real64 :: B_norm_exponent
    !sll_int32  :: QN_case
    
    !--> Pertubation
    sll_int32  :: perturb_choice  !< 1 mode or multiple modes
    sll_int32  :: mmode           !< theta mode
    sll_int32  :: nmode           !< z mode
    sll_real64 :: eps_perturb     !< amplitude

    !--> Algorithm
    sll_real64 :: dt
    sll_int32  :: number_iterations
    sll_int32  :: freq_diag_time  !< print scalar diagnostics every freq_diag_t
    sll_int32  :: freq_diag       !< print slices of distr. func. every freq_di

    ! Options (TODO: see if descriptors could be used here)
    character(len=256) :: advect2d_case 
    character(len=256) :: charac2d_case
    character(len=256) :: time_loop_case 
    character(len=256) :: poisson2d_case
    character(len=256) :: QN_case
    character(len=256) :: advector_x2
    character(len=256) :: advector_x3
    character(len=256) :: advector_x4
    character(len=256) :: interp_x1x2
    character(len=256) :: phi_interp_x1x2
    character(len=256) :: phi_interp_x3
    character(len=256) :: poisson2d_BC_rmin
    character(len=256) :: poisson2d_BC_rmax
    
    ! Order of 1D periodic advectors
    sll_int32          :: order_x2
    sll_int32          :: order_x3
    sll_int32          :: order_x4

    ! TODO: remove stencil info from namelists, we should calculate them
    sll_int32          :: lagrange_stencil_left
    sll_int32          :: lagrange_stencil_right
    sll_int32          :: deriv_stencil_left
    sll_int32          :: deriv_stencil_right
    sll_int32          :: poisson2d_BC(2)
    sll_int32          :: ierr
    !sll_int32          :: spline_degree
    
    ! Options (True/False)
    logical :: use_field_aligned_derivative    ! For electric field from phi
    logical :: use_field_aligned_interpolation ! FCISL general idea

    ! Namelists read from input file
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
      !Diota_file, &
      size_iota_file, &
      !size_Diota_file, &
      is_iota_file, &
      !is_Diota_file, &
      B0
            
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
      advector_x2, &
      advector_x3, &
      advector_x4, &
      order_x2, &
      order_x3, &
      order_x4, &
      lagrange_stencil_left, &
      lagrange_stencil_right, &
      deriv_stencil_left, &
      deriv_stencil_right, &
      poisson2d_case, &
      use_field_aligned_derivative, &
      use_field_aligned_interpolation

      !, spline_degree
    
    !==========================================================================
    !default parameters
    !==========================================================================
    
    iota0           = 0._f64
    Dr_iota0        = 0._f64
    iota_file       = "no_q_file.dat"
    !Diota_file      = "no_q_file.dat"
    is_iota_file    = .false.
    !is_Diota_file  = .false.
    size_iota_file  = 0
    !size_Diota_file = 0
    num_cells_x2    = 32
    
    use_field_aligned_derivative    = .false.
    use_field_aligned_interpolation = .false.
    
    !B_norm_exponent = -0.5_f64
    B0 = 1._f64
    
    !LAG5
    lagrange_stencil_left  = -2
    lagrange_stencil_right =  3
    !PPM1
    deriv_stencil_left  = -2
    deriv_stencil_right =  2
    
    ! TODO: some defaults missing
    
    !==========================================================================
    ! Read parameters from input file
    ! TODO: use Pierre's I/O
    !==========================================================================

    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
!       print *, '#init_dk4d_polar() failed to open file ', filename
       err_msg = 'failed to open file '// filename 
       SLL_ERROR( 'init_dk4d_field_aligned_polar', trim( err_msg ) )
    end if
    read(input_file,mesh)
    read(input_file,equilibrium)
    read(input_file,perturbation)
    read(input_file,sim_params)
    close(input_file)

    !==========================================================================
    ! Set initial time value (it may be != 0 if the simulation is restarted)
    !==========================================================================

    sim%time = 0.0_f64

    !==========================================================================
    ! Create various data structures and operators in simulation
    !==========================================================================

    !--> Mesh
    sim%m_x1 => sll_f_new_cartesian_mesh_1d(num_cells_x1,eta_min=r_min,eta_max=r_max)
    sim%m_x2 => sll_f_new_cartesian_mesh_1d(num_cells_x2,&
      eta_min=0._f64,eta_max=2._f64*sll_p_pi)
    sim%m_x3 => sll_f_new_cartesian_mesh_1d(num_cells_x3,eta_min=z_min,eta_max=z_max)
    sim%m_x4 => sll_f_new_cartesian_mesh_1d(num_cells_x4,eta_min=v_min,eta_max=v_max)
    
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
    sim%B0 = B0
    !sim%B_norm_exponent = B_norm_exponent

    ! TODO: use descriptors 
    ! type( SLL_BOUNDARY_CONDITION_DESCRIPTOR ) :: poisson2d_BC(2)
    ! ...
    ! call poisson2d_BC(1)%parse( poisson2d_BC_rmin )
    ! call poisson2d_BC(2)%parse( poisson2d_BC_rmax )

    select case (poisson2d_BC_rmin)
      case ("SLL_DIRICHLET")
        poisson2d_BC(1) = sll_p_dirichlet
      case ("SLL_NEUMANN")
        poisson2d_BC(1) = sll_p_neumann
      case ("SLL_NEUMANN_MODE_0")
        poisson2d_BC(1) = sll_p_neumann_mode_0      
      case default
        print *,'#bad choice for poisson2d_BC_rmin'
        print *,'#in init_dk4d_polar'
        stop
    end select

    select case (poisson2d_BC_rmax)
      case ("SLL_DIRICHLET")
        poisson2d_BC(2) = sll_p_dirichlet
      case ("SLL_NEUMANN")
        poisson2d_BC(2) = sll_p_neumann
      case ("SLL_NEUMANN_MODE_0")
        poisson2d_BC(2) = sll_p_neumann_mode_0
      case default
        print *,'#bad choice for poisson2d_BC_rmax'
        print *,'#in i nit_dk4d_polar'
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

    !--> Perturbation
    sim%perturb_choice = perturb_choice
    sim%mmode          = mmode
    sim%nmode          = nmode
    sim%eps_perturb    = eps_perturb

    !--> Algorithm
    sim%dt                 = dt
    sim%num_iterations     = number_iterations
    sim%freq_diag_time     = freq_diag_time
    sim%freq_diag     = freq_diag
    sim%use_field_aligned_derivative = use_field_aligned_derivative
    sim%use_field_aligned_interpolation = use_field_aligned_interpolation

    ! MPI info: no. of tasks and processor rank
    sim%world_size = sll_f_get_collective_size(sll_v_world_collective)
    sim%my_rank    = sll_f_get_collective_rank(sll_v_world_collective)

    ! Master prints info to standard output
    ! TODO: add subroutine: sim%print_info()
    if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
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
      print *,'#QN_case=',trim(QN_case)
      print *,'##perturbation'
      print *,'#perturb_choice=',perturb_choice
      print *,'#mmode=',mmode
      print *,'#nmode=',nmode
      print *,'#eps_perturb=',eps_perturb
      print *,'#dt=',dt
      print *,'#number_iterations=',number_iterations
      print *,'#time_loop_case=',trim(time_loop_case)
      print *,'#charac2d_case=',trim(charac2d_case)
      print *,'#use_field_aligned_derivative=',use_field_aligned_derivative
      print *,'#use_field_aligned_interpolation=',use_field_aligned_interpolation
    endif

    ! TODO: check if equivalent to
    ! call sim%m_x1%get_node_positions( sim%x1_node )
    call sll_o_get_node_positions(sim%m_x1,sim%x1_node)
    call sll_o_get_node_positions(sim%m_x2,sim%x2_node)
    call sll_o_get_node_positions(sim%m_x3,sim%x3_node)
    call sll_o_get_node_positions(sim%m_x4,sim%x4_node)

    ! Allocate iota_r
    SLL_ALLOCATE(sim%iota_r(num_cells_x1+1),ierr)
    !SLL_ALLOCATE(sim%Diota_r(num_cells_x1+1),ierr)
 
    ! Compute iota_r
    call initialize_iota_profile( &
      iota0, &    
      Dr_iota0, &
      iota_file, &
      !Diota_file, &
      size_iota_file, &
      !size_Diota_file, &
      is_iota_file, &
      !is_Diota_file, &
      sim%m_x1%num_cells+1, &
      sim%x1_node, &
      sim%iota_r) !, &
      !sim%Diota_r )

    ! Allocate B unit vector
    SLL_ALLOCATE(sim%b_unit_x2(num_cells_x1+1),ierr)
    SLL_ALLOCATE(sim%b_unit_x3(num_cells_x1+1),ierr)

    ! Compute B unit vector
    call initialize_b_unit_from_iota_profile( &
      sim%iota_r, &
      sim%x1_node, &
      sim%m_x1%num_cells+1, &
      sim%m_x3%eta_max-sim%m_x3%eta_min, &
      sim%b_unit_x2, &
      sim%b_unit_x3)

    ! Initialize equilibrium profiles
    call initialize_profiles_analytic(sim)

    ! Allocate f and fields with various layouts
    ! TODO: data used for initialization should not be inside sim
    call allocate_fdistribu4d_and_QN_DK_parx1(sim)

    ! Create (=allocate+initialize) Poisson's solver
    select case (poisson2d_case)
      case ("POLAR_FFT")
        ! Wrapper around a loop of 2D Poisson's solvers, with 3D interface
        ! Solves a quasi-neutral equation (of Poisson-type)
        ! Input/output data is parallelized in x1=r
        ! Internally it also uses parallelization in x2=theta
        ! Needs both par-x1 and par-x2 layouts
        sim%poisson3d => sll_f_new_qn_solver_3d_polar_parallel_x1_wrapper( &
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
          inv_Te=1._f64/sim%Te_r )

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

    ! Create interpolators in (r,theta) plane for f distribution and E field
    ! Only option available: tensor-product cubic splines
    ! Periodic BCs in theta, Hermite BC in r (~Neumann: 1st derivative = 0) 
    ! TODO: check if linear interpolation on uniform grid is available in
    !       - interpolation/sll_lagrange_interpolator_1d.F90, or
    !       - sll_m_lagrange_interpolation/sll_lagrange_interpolation.F90
    select case (interp_x1x2)
      case ("SLL_CUBIC_SPLINES")
        f_interp2d => sll_f_new_cubic_spline_interpolator_2d( &
          sim%m_x1%num_cells+1, &
          sim%m_x2%num_cells+1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sim%m_x2%eta_min, &
          sim%m_x2%eta_max, &
          sll_p_hermite, &
          sll_p_periodic, &
          const_eta1_min_slope = 0._f64, & !to prevent problem on the boundary
          const_eta1_max_slope = 0._f64)
        A1_interp1d_x1 => sll_f_new_cubic_spline_interpolator_1d( &
          sim%m_x1%num_cells+1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sll_p_hermite, &
          slope_left  = 0._f64, &
          slope_right = 0._f64 )
        A2_interp1d_x1 => sll_f_new_cubic_spline_interpolator_1d( &
          sim%m_x1%num_cells+1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sll_p_hermite, &
          slope_left  = 0._f64, &
          slope_right = 0._f64 )
        A1_interp2d => sll_f_new_cubic_spline_interpolator_2d( &
          sim%m_x1%num_cells+1, &
          sim%m_x2%num_cells+1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sim%m_x2%eta_min, &
          sim%m_x2%eta_max, &
          sll_p_hermite, &
          sll_p_periodic, &
          const_eta1_min_slope = 0._f64, &
          const_eta1_max_slope = 0._f64)
        A2_interp2d => sll_f_new_cubic_spline_interpolator_2d( &
          sim%m_x1%num_cells+1, &
          sim%m_x2%num_cells+1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sim%m_x2%eta_min, &
          sim%m_x2%eta_max, &
          sll_p_hermite, &
          sll_p_periodic, &
          const_eta1_min_slope = 0._f64, &
          const_eta1_max_slope = 0._f64)
      case default
        print *,'#bad interp_x1x2',interp_x1x2
        print *,'#not implemented'
        print *,'#in init_dk4d_polar'
        stop
    end select

    ! Create (r,theta) interpolator for phi, used ONLY to obtain E = -dphi/dx
    ! For better mass conservation, here we should use the same interpolator
    ! used for f
    ! NOTE: BC at r=r_min only consistent with sll_p_neumann_mode_0 option
    select case (phi_interp_x1x2)
      case ("SLL_CUBIC_SPLINES")
         sim%phi_interp_x1x2 => sll_f_new_cubic_spline_interpolator_2d( &
          sim%m_x1%num_cells+1, &
          sim%m_x2%num_cells+1, &
          sim%m_x1%eta_min, &
          sim%m_x1%eta_max, &
          sim%m_x2%eta_min, &
          sim%m_x2%eta_max, &
          sll_p_hermite, &
          sll_p_periodic, &
          const_eta1_min_slope = 0._f64, & !to prevent problem on the boundary
          const_eta1_max_slope = 0._f64)
      case default
        print *,'#bad phi_interp_x1x2',phi_interp_x1x2
        print *,'#not implemented'
        print *,'#in init_dk4d_polar'
        stop
    end select

    ! Interpolator of phi along z.  This is used to compute Ez = -dphi/dz,
    ! when field-aligned interpolation option is not chosen.
    select case (phi_interp_x3)
      case ("SLL_CUBIC_SPLINES")
        sim%phi_interp_x3 => sll_f_new_cubic_spline_interpolator_1d( &
          sim%m_x3%num_cells+1, &
          sim%m_x3%eta_min, &
          sim%m_x3%eta_max, &
          sll_p_periodic)
      case default
        print *,'#bad phi_interp_x3',phi_interp_x3
        print *,'#not implemented'
        print *,'#in init_dk4d_polar'
        stop
    end select

    ! Computation of 2D characteristics
    ! Given (r,theta) grid, gives initial position of each point 
    ! TODO: implement wrapper to explicit RK integrators
    select case (charac2d_case)
      case("SLL_CHARAC_EULER")
        charac2d => sll_f_new_explicit_euler_2d_charac(&
          sim%m_x1%num_cells+1, &
          sim%m_x2%num_cells+1, &
          bc_type_1= sll_p_set_to_limit, &
          bc_type_2= sll_p_periodic, &
          eta1_min = sim%m_x1%eta_min, &
          eta1_max = sim%m_x1%eta_max, &
          eta2_min = sim%m_x2%eta_min, &
          eta2_max = sim%m_x2%eta_max)
      case("SLL_CHARAC_VERLET")
        charac2d => sll_f_new_verlet_2d_charac(&
          sim%m_x1%num_cells+1, &
          sim%m_x2%num_cells+1, &
          A1_interp2d, &
          A2_interp2d, &
          A1_interp1d_x1, &
          A2_interp1d_x1, &
          bc_type_1= sll_p_set_to_limit, &
          bc_type_2= sll_p_periodic, &
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

    ! 2D advector in (r,theta): given f and poloidal A=(A1,A2), advance by dt
    select case(advect2d_case)
      case ("SLL_BSL")
      sim%adv_x1x2 => sll_f_new_bsl_2d_advector(&
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

    ! 1D advector in theta: needed by field-aligned algorithm
    select case (advector_x2)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%adv_x2 => sll_f_new_periodic_1d_advector( &
          sim%m_x2%num_cells, &
          sim%m_x2%eta_min, &
          sim%m_x2%eta_max, &
          sll_p_spline, & 
          order_x2) 
      case("SLL_LAGRANGE") ! arbitrary order sll_p_lagrange periodic interpolation
        sim%adv_x2 => sll_f_new_periodic_1d_advector( &
          sim%m_x2%num_cells, &
          sim%m_x2%eta_min, &
          sim%m_x2%eta_max, &
          sll_p_lagrange, & 
          order_x2) 
       case default
         print*,'#advector in x2', advector_x2, ' not implemented'
         stop 
    end select

    ! 1D advector along z: not for field-aligned
    select case (advector_x3)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%adv_x3 => sll_f_new_periodic_1d_advector( &
          sim%m_x3%num_cells, &
          sim%m_x3%eta_min, &
          sim%m_x3%eta_max, &
          sll_p_spline, & 
          order_x3) 
      case("SLL_LAGRANGE") ! arbitrary order sll_p_lagrange periodic interpolation
        sim%adv_x3 => sll_f_new_periodic_1d_advector( &
          sim%m_x3%num_cells, &
          sim%m_x3%eta_min, &
          sim%m_x3%eta_max, &
          sll_p_lagrange, & 
          order_x3) 
       case default
         print*,'#advector in x3', advector_x3, ' not implemented'
         stop 
    end select

    ! 1D advector along v
    select case (advector_x4)
      case ("SLL_SPLINES") ! arbitrary order periodic splines
        sim%adv_x4 => sll_f_new_periodic_1d_advector( &
          sim%m_x4%num_cells, &
          sim%m_x4%eta_min, &
          sim%m_x4%eta_max, &
          sll_p_spline, & 
          order_x4) 
      case("SLL_LAGRANGE") ! arbitrary order sll_p_lagrange periodic interpolation
        sim%adv_x4 => sll_f_new_periodic_1d_advector( &
          sim%m_x4%num_cells, &
          sim%m_x4%eta_min, &
          sim%m_x4%eta_max, &
          sll_p_lagrange, & 
          order_x4) 
       case default
         print*,'#advector in x4', advector_x4, ' not implemented'
         stop 
    end select
    
    ! Create object for computing field aligned derivatives
    sim%deriv => sll_f_new_oblic_2d_derivative( &
      sim%m_x2%num_cells, &
      sim%adv_x2, &
      sim%m_x3%num_cells, &
      sim%m_x3%eta_min, &
      sim%m_x3%eta_max, &
      deriv_stencil_left, &
      deriv_stencil_right )

    ! 2D advector, field aligned
    ! TODO: add Hermite interpolation option
    sim%adv_x2x3 => sll_f_new_oblic_2d_advector( &
      sim%m_x2%num_cells, &
      sim%adv_x2, &
      sim%m_x3%num_cells, &
      sim%m_x3%eta_min, &
      sim%m_x3%eta_max, &
      lagrange_stencil_left, &
      lagrange_stencil_right )

  end subroutine init_dk4d_field_aligned_polar

!==============================================================================
! RUN
!==============================================================================

  subroutine run_dk4d_field_aligned_polar(sim)
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar), intent(inout) :: sim

    !--> For initial profile HDF5 saving
    !sll_int32                    :: file_err
    !sll_int32                    :: file_id
    character(len=12), parameter :: filename_prof = "init_prof.h5"
    sll_real64, allocatable      :: f4d_store(:,:,:,:)  ! for field prediction
    sll_int32  :: loc4d_sz_x1
    sll_int32  :: loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3
    sll_int32  :: loc4d_sz_x4
    sll_int32  :: loc4d(4)
    sll_int32  :: iter
    sll_int32  :: nc_x1
    sll_int32  :: nc_x2
    sll_int32  :: nc_x3
    sll_int32  :: nc_x4
    sll_int32  :: ierr
    sll_real64 :: dt
    sll_int32  :: th_diag_id
    sll_int32  :: i_plot 

    ! Output
    character(len=4)        :: cplot
    character(len=*), parameter :: filetype="HDF"
    character(len=64)       :: hdf5_file_name, xml_file_name, field_name
    !sll_int32               :: hdf5_file_id
    !sll_int32               :: error
    type( sll_t_hdf5_serial ) :: hdf5_file
    type( sll_t_xdmf_parallel_file ) :: xdmf_file
    logical                 :: to_file
    character(len=256)      :: field_path
    character(len=256)      :: dataset_x1_polar, dataset_x2_polar
    character(len=256)      :: dataset_x2_cart , dataset_x3_cart
    sll_int32               :: dims_x1x2(2), dims_x2x3(2)
    sll_int32               ::  gid_x1x2   ,  gid_x2x3

    ! Shortcuts 
    dt = sim%dt    
    nc_x1 = sim%m_x1%num_cells ! GLOBAL number of cells along each direction
    nc_x2 = sim%m_x2%num_cells
    nc_x3 = sim%m_x3%num_cells
    nc_x4 = sim%m_x4%num_cells
    
    !*** Saving of the radial profiles in HDF5 file ***
    if (sll_f_get_collective_rank(sll_v_world_collective)==0) then
      call hdf5_file%init( filename_prof )
      call hdf5_file%create()
      call hdf5_file%write_array( sim%n0_r, 'n0_r' )
      call hdf5_file%write_array( sim%Ti_r, 'Ti_r' )
      call hdf5_file%write_array( sim%Te_r, 'Te_r' )
      call hdf5_file%delete()

      ierr = 1
      call sll_o_gnuplot_1d(sim%n0_r,'n0_r_init',ierr)
      call sll_o_gnuplot_1d(sim%Ti_r,'Ti_r_init',ierr)
      call sll_o_gnuplot_1d(sim%Te_r,'Te_r_init',ierr)
      call sll_o_gnuplot_1d(sim%iota_r,'iota_r_init',ierr)
      call sll_o_gnuplot_1d(sim%b_unit_x2,'b_theta_init',ierr)
      call sll_o_gnuplot_1d(sim%b_unit_x3,'b_z_init',ierr)
    end if

    ! Allocate temp. storage, distributed in r
    call sll_o_compute_local_sizes( sim%layout4d_parx1, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
    SLL_ALLOCATE(f4d_store(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)

    ! Master MPI task writes time diagnostics (scalar quantities) to file
    if(sll_f_get_collective_rank(sll_v_world_collective)==0) then
      call sll_s_ascii_file_create('thdiag.dat', th_diag_id, ierr)
    endif

    ! Initialize all parallel arrays
    call initialize_fdistribu4d_DK(sim,sim%layout4d_parx1,sim%f4d_parx1)

    !--------------------------------------------------------------------------
    ! Write mesh data to file (it does not change with time)
    !--------------------------------------------------------------------------

    ! (r,theta) slice: poloidal plane
    dims_x1x2 = [sim%m_x1%num_cells+1, sim%m_x2%num_cells+1]
    if(sll_f_get_collective_rank(sll_v_world_collective)==0) then
      call write_mesh_x1x2_polar( sim%m_x1, sim%m_x2, &
        dataset_x1_polar, dataset_x2_polar )
    else     
      call write_mesh_x1x2_polar_light( sim%m_x1, sim%m_x2, &
        dataset_x1_polar, dataset_x2_polar )
    endif
    

    ! (theta,z) slice: flux surface
    dims_x2x3 = [sim%m_x2%num_cells+1, sim%m_x3%num_cells+1]
    if(sll_f_get_collective_rank(sll_v_world_collective)==0) then
      call write_mesh_x2x3_cart( sim%m_x2, sim%m_x3, &
        dataset_x2_cart , dataset_x3_cart  )
    else
      call write_mesh_x2x3_cart_light( sim%m_x2, sim%m_x3, &
        dataset_x2_cart , dataset_x3_cart  )        
    endif
    !--------------------------------------------------------------------------
    ! Time cycle
    !--------------------------------------------------------------------------
    i_plot = 0
    do iter=1,sim%num_iterations

      ! Compute 3D charge density from 4D distribution function
      call sll_o_compute_local_sizes( sim%layout4d_parx1, &
        loc4d_sz_x1, &
        loc4d_sz_x2, &
        loc4d_sz_x3, &
        loc4d_sz_x4 )
      call compute_rho_dk(sim)    

      if(modulo(iter,sim%freq_diag)==0)then
        if(sll_f_get_collective_rank(sll_v_world_collective)==0) then
          print*,'#iteration=',iter      
        endif
      endif

      ! Solve QN Poisson's eq. and compute E field from phi
      call solve_quasi_neutral_parx1( sim )  ! updates sim%phi3d_parx1
      call compute_field_dk_parx1( sim )     ! updates sim%A1_parx3, A2_parx3, A3_parx1

      ! Print time diagnostics
      if(modulo(iter,sim%freq_diag_time)==0)then
        call time_history_diagnostic_dk_polar( sim, th_diag_id, iter-1 )
      endif

      ! Evolve full solution by one time-step (operator splitting is used)
      select case (sim%time_case)

        ! Option a) Lie splitting, 1st order
        case (SLL_TIME_LOOP_EULER)  
          call advection_x3( sim, dt )
          call advection_x4( sim, dt )
          call advection_x1x2( sim, dt )

        ! Option b) Predictor-corrector plus Strang splitting, 2nd order
        case (SLL_TIME_LOOP_PREDICTOR_CORRECTOR)
          !prediction (Lie, 1st order)
          f4d_store = sim%f4d_parx1
          call advection_x3( sim, 0.5_f64*dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x1x2( sim, 0.5_f64*dt )
          call compute_rho_dk(sim)  
          call solve_quasi_neutral_parx1( sim )
          call compute_field_dk_parx1( sim )          
          !correction (Strang, 2nd order)
          sim%f4d_parx1 = f4d_store
          call advection_x3( sim, 0.5_f64*dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x1x2( sim, dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x3( sim, 0.5_f64*dt )

        ! Option c) As above, but with 2nd order prediction (more expensive)
        case (SLL_TIME_LOOP_PREDICTOR2_CORRECTOR)
          !prediction (Strang, 2nd order)
          f4d_store = sim%f4d_parx1
          call advection_x3( sim, 0.25_f64*dt )
          call advection_x4( sim, 0.25_f64*dt )
          call advection_x1x2( sim, 0.5_f64*dt )
          call advection_x4( sim, 0.25_f64*dt )
          call advection_x3( sim, 0.25_f64*dt )
          call compute_rho_dk(sim)  
          call solve_quasi_neutral_parx1( sim )
          call compute_field_dk_parx1( sim )          
          !correction (Strang, 2nd order)
          sim%f4d_parx1 = f4d_store
          call advection_x3( sim, 0.5_f64*dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x1x2( sim, dt )
          call advection_x4( sim, 0.5_f64*dt )
          call advection_x3( sim, 0.5_f64*dt )

        case default
          call sll_s_halt_collective()
          print *,'#sim%time_case=',sim%time_case
          print *, '#not implemented'
          print *,'#in run_dk4d_polar'
          stop
      end select          

      !------------------------------------------------------------------------
      ! Update time value
      !------------------------------------------------------------------------

      ! TODO: time should be incremented, but roundoff could be an issue
      sim%time = iter*dt

      !------------------------------------------------------------------------
      ! Write "heavy" diagnostics
      !------------------------------------------------------------------------
      if (modulo(iter,sim%freq_diag)==0) then

        i_plot = i_plot+1

        ! Determine file names
        call sll_s_int2string( i_plot, cplot )
        hdf5_file_name = 'diag2d_'//cplot//'.h5'  ! HDF5 file (heavy data)
        xml_file_name  = 'diag2d_'//cplot//'.xmf' ! XML  file (light data)

        ! Initialize HDF5 wrapper
        call hdf5_file%init( hdf5_file_name )

        ! MASTER: Create new HDF5 file and close it
        if (sll_f_get_collective_rank(sll_v_world_collective)==0) then
          call hdf5_file%create()
          call hdf5_file%close()
        end if

        ! Wait for HDF5 file to be available to other processors
        call sll_s_collective_barrier( sll_v_world_collective )

        ! Initialize parallel XDMF file
        call xdmf_file%init( sim%time, sll_v_world_collective )

        ! Add cartesian grid (x2,x3) to XDMF file
        call xdmf_file%add_grid( &
          'mesh_x2x3_cart', &
          dataset_x2_cart, &
          dataset_x3_cart, &
          dims_x2x3, &
          gid_x2x3 )

        ! Add polar grid (x1,x2) to XDMF file
        call xdmf_file%add_grid( &
          'mesh_x1x2_polar', &
          dataset_x1_polar, &
          dataset_x2_polar, &
          dims_x1x2, &
          gid_x1x2 )

        !----------------------------------------------------------------------
        ! Flux surface plot: f_x2x3
        !----------------------------------------------------------------------

        ! Get global indices of point that has:
        !  . r = (r_min+r_max)/2  (center value)
        !  . theta = 0
        !  . z = 0
        !  . v = 0 (center value)
        loc4d(1:4)  = sll_o_global_to_local( &
          sim%layout4d_parx1, &
          [nc_x1/2+1,1,1,nc_x4/2+1] )

        ! If point is in local domain, print slice of f on flux surface
        if(loc4d(1) > 0) then
          field_name = 'f_x2x3'
          call hdf5_file%open()
          call hdf5_file%write_array( sim%f4d_parx1(loc4d(1),:,:,loc4d(4)), &
                                      '/'//field_name )
          call hdf5_file%close()
          to_file    = .true.
          field_path = trim( hdf5_file_name )//':/'//trim( field_name )
        else
          to_file    = .false.
          field_path = ''
        end if

        ! Add field to appropriate grid in XDMF file
        call xdmf_file%add_field( gid_x2x3, field_name, field_path, to_file )

        !----------------------------------------------------------------------
        ! Polar plot: f_x1x2
        !----------------------------------------------------------------------

        ! Copy f distributed in x1 into f distributed in (x3,x4)
        call sll_o_apply_remap_4d( &
          sim%remap_plan_parx1_to_parx3x4, &
          sim%f4d_parx1, &
          sim%f4d_parx3x4 )

        ! Get global indices of point that has:
        !  . r = r_min
        !  . theta = 0
        !  . z = 0
        !  . v = 0 (center value)
        loc4d(1:4) = sll_o_global_to_local( sim%layout4d_parx3x4, [1,1,1,nc_x4/2+1] )

        ! If point is in local domain, print slice of f on poloidal plane z=0
        if(loc4d(3) > 0) then
#ifndef NOHDF5
          field_name = 'f_x1x2'
          call hdf5_file%open()
          call hdf5_file%write_array( sim%f4d_parx3x4(:,:,loc4d(3),loc4d(4)), &
                                      '/'//field_name )
          call hdf5_file%close()
          to_file    = .true.
          field_path = trim( hdf5_file_name )//':/'//trim( field_name )
#endif
        else
          to_file    = .false.
          field_path = ''
        endif

        ! Add field to appropriate grid in XDMF file
        call xdmf_file%add_field( gid_x1x2, field_name, field_path, to_file )

        !----------------------------------------------------------------------
        ! Write XML document, then sll_o_delete XDMF file contents
        !----------------------------------------------------------------------
        call xdmf_file%write( xml_file_name )
        call xdmf_file%delete()

        ! sll_o_delete HDF5 wrapper (will be initialized again)
        call hdf5_file%delete()

      endif !if (modulo(iter,sim%freq_diag)==0)

    enddo !do iter=1,sim%num_iterations
   
  end subroutine run_dk4d_field_aligned_polar

!==============================================================================

!> initialize q_profile
!> if is_q_file=.false. and is_Dq_file=.false., takes q(r) = q0_r+r*Dq0_r
!> if is_q_file=.false. and is_Dq_file=.true., error
!> if is_q_file=.true. and is_Dq_file=.false., takes q(r) from q_file Dq(r) by interp of q
!> if is_q_file=.true. and is_Dq_file=.false., takes q(r) from q_file Dq(r) from Dq_file
!> when q (or Dq) are taken from q_file, we use cubic sll_p_spline interpolation for getting
!> values on grid points
  subroutine initialize_iota_profile( &
    iota0_r, &
    Diota0_r, &
    iota_file, &
    !Diota_file, &
    size_iota_file, &
    !size_Diota_file, &
    is_iota_file, &
    !is_Diota_file, &
    num_points_r, &
    r_array, &
    iota) !, &
    !Diota )

    sll_real64,         intent(in)  :: iota0_r
    sll_real64,         intent(in)  :: Diota0_r
    character(len=256), intent(in)  :: iota_file
    !character(len=256), intent(in)  :: Diota_file
    sll_int32,          intent(in)  :: size_iota_file
    !sll_int32,          intent(in)  :: size_Diota_file
    logical,            intent(in)  :: is_iota_file
    !logical,            intent(in)  :: is_Diota_file
    sll_int32,          intent(in)  :: num_points_r
    sll_real64,         intent(in)  :: r_array(:)
    sll_real64,         intent(out) :: iota(:)
    !sll_real64,         intent(out) :: Diota(:)

    !local variables
    sll_int32 :: i
    
    ! some checking
    if (size(iota)<num_points_r) then
      call sll_s_halt_collective()
      print *,'#bad size for iota in initialize_iota_profile'
      stop
    endif
    !if(size(Diota)<num_points_r)then
    !  call sll_s_halt_collective()
    !  print *,'#bad size for Diota in initialize_iota_profile'
    !  stop
    !endif
    
    if ((is_iota_file .eqv. .false.)) then ! .and. (is_Diota_file .eqv. .false.)) then
      if(size(r_array)<num_points_r)then
        call sll_s_halt_collective()
        print *,'#bad size for r_array in initialize_iota_profile'
        stop
      endif
      do i=1,num_points_r
        iota(i) = iota0_r+Diota0_r*r_array(i)
        !Diota(i) = Diota0_r        
      enddo
    endif

    !if ((is_iota_file .eqv. .false.) .and. (is_Diota_file .eqv. .true.)) then
    !  call sll_s_halt_collective()
    !  print *,'#bad value for is_iota_file and is_Diota_file in initialize_iota_profile'
    !  stop      
    !endif

    !if ((is_iota_file .eqv. .true.) .and. (is_Diota_file .eqv. .false.))then
    !  call sll_s_halt_collective()
    !  print *,'#not implemented for the moment'
    !  print *,'#in initialize_iota_profile'
    !  stop
    !endif

    if ((is_iota_file .eqv. .true.) ) then !.and. (is_Diota_file .eqv. .true.))then
      call sll_s_halt_collective()
      print *,'#not implemented for the moment'
      print *,'#in initialize_iota_profile'
      stop
    endif
    
  end subroutine initialize_iota_profile
 

!  subroutine initialize_c_from_iota_profile( &
!    iota, &
!    r_array, &
!    num_points_r, &
!    L, &
!    c_r)
!    sll_real64, dimension(:), intent(in) :: iota
!    sll_real64, dimension(:), intent(in) :: r_array
!    sll_int32, intent(in) :: num_points_r
!    sll_real64, intent(in) :: L
!    sll_real64, dimension(:), intent(out) :: c_r
!    !local variables
!    sll_int32 :: i
!    sll_real64 :: big_R
!
!    ! some checking
!    if(size(iota)<num_points_r)then
!      call sll_s_halt_collective()
!      print *,'#bad size for iota in initialize_c_from_iota_profile'
!      stop
!    endif
!
!    if(size(r_array)<num_points_r)then
!      call sll_s_halt_collective()
!      print *,'#bad size for r_array in initialize_c_from_iota_profile'
!      stop
!    endif
!
!    if(size(c_r)<num_points_r)then
!      call sll_s_halt_collective()
!      print *,'#bad size for c_r in initialize_c_from_iota_profile'
!      stop
!    endif
!
!    big_R = L/(2._f64*sll_p_pi)
!
!    do i=1,num_points_r
!      c_r(i) = r_array(i)*iota(i)/big_R
!    enddo
!
!  end subroutine initialize_c_from_iota_profile


  subroutine initialize_b_unit_from_iota_profile( &
    iota, &
    r_array, &
    num_points_r, &
    L, &
    b_unit_x2, &
    b_unit_x3)
    
    sll_real64, intent(in)  :: iota(:)
    sll_real64, intent(in)  :: r_array(:)
    sll_int32,  intent(in)  :: num_points_r
    sll_real64, intent(in)  :: L
    sll_real64, intent(out) :: b_unit_x2(:)
    sll_real64, intent(out) :: b_unit_x3(:)

    !local variables
    sll_int32  :: i
    sll_real64 :: big_R
    sll_real64 :: c

    ! some checking
    if(size(iota)<num_points_r)then
      call sll_s_halt_collective()
      print *,'#bad size for iota in initialize_b_unit_from_iota_profile'
      stop
    endif

    if(size(r_array)<num_points_r)then
      call sll_s_halt_collective()
      print *,'#bad size for r_array in initialize_b_unit_from_iota_profile'
      stop
    endif

    if(size(b_unit_x2)<num_points_r)then
      call sll_s_halt_collective()
      print *,'#bad size for b_unit_x2 in initialize_b_unit_from_iota_profile'
      stop
    endif

    if(size(b_unit_x3)<num_points_r)then
      call sll_s_halt_collective()
      print *,'#bad size for b_unit_x3 in initialize_b_unit_from_iota_profile'
      stop
    endif

    big_R = L/(2._f64*sll_p_pi)
    
    do i=1,num_points_r
      c = r_array(i)*iota(i)/big_R
      b_unit_x2(i) = (c/sqrt(1._f64+c**2))/r_array(i)
      b_unit_x3(i) = 1._f64/sqrt(1._f64+c**2)
    enddo
    
  end subroutine initialize_b_unit_from_iota_profile


!  subroutine initialize_iota_modif( &
!    Nc_x1, &
!    Nc_x2, &
!    iota, &
!    num_points_r, &
!    spaghetti_size_guess, &
!    spaghetti_size, &
!    shift_r)
!    sll_int32, intent(in) :: Nc_x1
!    sll_int32, intent(in) :: Nc_x2
!    sll_real64, dimension(:), intent(in) :: iota
!    sll_int32, intent(in) :: num_points_r
!    sll_int32, intent(in) :: spaghetti_size_guess
!    sll_int32, intent(out) :: spaghetti_size
!    sll_int32, dimension(:), intent(out) :: shift_r
!    !local variables
!    sll_int32 :: i
!    sll_real64 :: big_R
!    sll_int32 :: spaghetti_size0
!
!    if(size(iota)<num_points_r)then
!      call sll_s_halt_collective()
!      print *,'#bad size for iota in initialize_iota_modif'
!      stop
!    endif
!
!    if(size(shift_r)<num_points_r)then
!      call sll_s_halt_collective()
!      print *,'#bad size for shift_R in initialize_iota_modif'
!      stop
!    endif
!
!    do i=1,num_points_r
!      call sll_s_compute_spaghetti_and_shift_from_guess( &
!        Nc_x1, &
!        Nc_x2, &
!        iota(i), &
!        spaghetti_size_guess, &
!        shift_r(i), &
!        spaghetti_size)
!      if(i==1)then
!        spaghetti_size0 = spaghetti_size
!      endif
!      if(spaghetti_size .ne. spaghetti_size0)then
!        print *,'#bad spaghetti size in initialize_iota_modif'
!        print *,'#we want to have same spaghetti_size for all the r'
!      endif  
!    enddo
!    
!
!  end subroutine initialize_iota_modif


  subroutine time_history_diagnostic_dk_polar( sim, file_id, step )
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar), &
               intent(in) :: sim
    sll_int32, intent(in) :: file_id
    sll_int32, intent(in) :: step

    sll_int32  :: Nc_x1
    sll_int32  :: Nc_x2
    sll_real64 :: nrj
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_int32  :: ierr

    Nc_x1  = sim%m_x1%num_cells
    Nc_x2  = sim%m_x2%num_cells
    delta1 = sim%m_x1%delta_eta
    delta2 = sim%m_x2%delta_eta

    nrj = 0._f64
    call sll_s_compute_reduction_2d_to_0d(&
      sim%phi3d_parx3(:,:,1)**2, &
      nrj, &
      Nc_x1+1, &
      Nc_x2+1, &
      delta1, &    
      delta2)

    if (sll_f_get_collective_rank(sll_v_world_collective)==0) then
      write(file_id,'(f12.5,2g20.12)') real(step,f64)*sim%dt, nrj
      if (step==0) then
        ierr = 1
        call sll_o_gnuplot_1d(sim%phi3d_parx3(:,1,1),'phi_0',ierr)
        !call sll_o_gnuplot_1d(sim%rho3d_seqx1x2(:,1,1)/sim%n0_r(:)-1._f64,'rho_0',ierr)
        call sll_o_gnuplot_1d(sim%Ti_r(:),'Ti_r',ierr)
        call sll_o_gnuplot_1d(sim%Te_r(:),'Te_r',ierr)
        call sll_o_gnuplot_1d(sim%n0_r(:),'n0_r',ierr)
      endif
    endif
 
  end subroutine time_history_diagnostic_dk_polar


  subroutine compute_field_from_phi_polar(phi,mesh1,mesh2,A1,A2,interp2d,B0)
    sll_real64,                      intent(in)  :: phi(:,:)
    type(sll_t_cartesian_mesh_1d),     pointer     :: mesh1
    type(sll_t_cartesian_mesh_1d),     pointer     :: mesh2
    sll_real64,                      intent(out) :: A1(:,:)
    sll_real64,                      intent(out) :: A2(:,:)
    class(sll_c_interpolator_2d), pointer     :: interp2d
    sll_real64,                      intent(in)  :: B0

    sll_int32  :: Nc_x1
    sll_int32  :: Nc_x2
    sll_real64 :: x1_min
    sll_real64 :: x2_min
    sll_real64 :: delta_x1
    sll_real64 :: delta_x2
    sll_real64 :: x1
    sll_real64 :: x2
    sll_int32  :: i1
    sll_int32  :: i2
    
    Nc_x1    = mesh1%num_cells
    Nc_x2    = mesh2%num_cells
    x1_min   = mesh1%eta_min
    x2_min   = mesh2%eta_min
    delta_x1 = mesh1%delta_eta
    delta_x2 = mesh2%delta_eta

    call interp2d%compute_interpolants(phi)

    do i2=1,Nc_x2+1
      x2=x2_min+real(i2-1,f64)*delta_x2
      do i1=1,Nc_x1+1
        x1=x1_min+real(i1-1,f64)*delta_x1
        A1(i1,i2)= interp2d%interpolate_from_interpolant_derivative_eta2(x1,x2)/(x1*B0)
        A2(i1,i2)=-interp2d%interpolate_from_interpolant_derivative_eta1(x1,x2)/(x1*B0)
      end do
    end do
    
  end subroutine compute_field_from_phi_polar


  subroutine compute_field_from_phi_cartesian_1d( phi, mesh, A, interp )
    sll_real64,                      intent(in)  :: phi(:)
    type(sll_t_cartesian_mesh_1d),     pointer     :: mesh
    sll_real64,                      intent(out) :: A(:)
    class(sll_c_interpolator_1d), pointer     :: interp

    sll_int32  :: Nc_x1
    sll_real64 :: x1_min
    sll_real64 :: delta_x1
    sll_real64 :: x1
    sll_int32  :: i1
    
    Nc_x1    = mesh%num_cells
    x1_min   = mesh%eta_min
    delta_x1 = mesh%delta_eta

    call interp%compute_interpolants(phi)

    do i1=1,Nc_x1+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      A(i1) = interp%interpolate_from_interpolant_derivative_eta1(x1)
    end do

  end subroutine compute_field_from_phi_cartesian_1d

  !----------------------------------------------------------------------------
  ! Input:             Output:
  !  - phi3d_parx3       - A1_parx3, A2_parx3
  !  - phi3d_parx1       - A3_parx1
  !----------------------------------------------------------------------------
  subroutine compute_field_dk_parx1( sim )
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar) :: sim

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
    sll_int32 :: glob_ind(4)
    
    nc_x1 = sim%m_x1%num_cells
    nc_x2 = sim%m_x2%num_cells
    nc_x3 = sim%m_x3%num_cells
    
    ! Use layout parallel in (z,v)
    call sll_o_compute_local_sizes( &
      sim%layout4d_parx3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
    
    !print *,'#phi3d_parx3=',maxval(sim%phi3d_parx3),minval(sim%phi3d_parx3), &
    !  sll_f_get_collective_rank(sll_v_world_collective)

    ! Compute A1_parx3 and A2_parx3 from phi3d_parx3
    ! 2D gradient of phi in poloidal plane (r,theta)
    do i3 = 1,loc4d_sz_x3
      call compute_field_from_phi_polar( &
        sim%phi3d_parx3(1:nc_x1+1,1:nc_x2+1,i3), &
        sim%m_x1, &
        sim%m_x2, &
        sim%A1_parx3(1:nc_x1+1,1:nc_x2+1,i3), &
        sim%A2_parx3(1:nc_x1+1,1:nc_x2+1,i3), &
        sim%phi_interp_x1x2, &
        sim%B0)
    enddo

    ! Use layout parallel in r (we now work on 2D flux surfaces)
    call sll_o_compute_local_sizes( &
      sim%layout4d_parx1, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )

    ! Compute A3_parx1 from phi3d_parx1
    ! Option a) Field-aligned algorithm on (theta,z) flux surface
    if (sim%use_field_aligned_derivative .eqv. .true.) then
      do i1=1, loc4d_sz_x1
        glob_ind(1:4) = sll_o_local_to_global( sim%layout4d_parx1, [i1,1,1,1] )
        call sll_s_compute_oblic_derivative_2d( &
          sim%deriv, &
          sim%b_unit_x2(glob_ind(1)), &
          sim%b_unit_x3(glob_ind(1)), &
          sim%phi3d_parx1(i1,1:nc_x2+1,1:nc_x3+1), &
          sim%A3_parx1(i1,1:nc_x2+1,1:nc_x3+1))
        sim%A3_parx1(i1,1:nc_x2+1,1:nc_x3+1) = -sim%A3_parx1(i1,1:nc_x2+1,1:nc_x3+1)
      enddo
      
    ! Option b) Just take 1D derivative of phi along z
    else
      do i2=1, loc4d_sz_x2
        do i1=1, loc4d_sz_x1
          call compute_field_from_phi_cartesian_1d( &
            sim%phi3d_parx1(i1,i2,1:nc_x3+1), &
            sim%m_x3, &
            sim%A3_parx1(i1,i2,1:nc_x3+1), &
            sim%phi_interp_x3)
          sim%A3_parx1(i1,i2,1:nc_x3+1) = -sim%A3_parx1(i1,i2,1:nc_x3+1)
        enddo
      enddo
    endif

  end subroutine compute_field_dk_parx1

!------------------------------------------------------------------------------
! Compute 3D charge density from 4D distribution function
!------------------------------------------------------------------------------
  subroutine compute_rho_dk( sim )
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar) :: sim

    sll_int32 :: loc_sz_x1
    sll_int32 :: loc_sz_x2
    sll_int32 :: loc_sz_x3
    sll_int32 :: loc_sz_x4
    
    call sll_o_compute_local_sizes( sim%layout4d_parx1, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4 )

    call sll_s_compute_reduction_4d_to_3d_direction4(&
      sim%f4d_parx1, &
      sim%rho3d_parx1, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4, &
      sim%m_x4%delta_eta)

  end subroutine compute_rho_dk

!------------------------------------------------------------------------------

  subroutine advection_x3( sim, dt )
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar), &
                intent(inout) :: sim
    sll_real64, intent(in)    :: dt

    sll_real64, allocatable :: f1d(:)
    sll_real64, allocatable :: f2d(:,:)
    sll_int32  :: nc_x1
    sll_int32  :: nc_x2
    sll_int32  :: nc_x3
    sll_int32  :: nc_x4
    sll_int32  :: i1
    sll_int32  :: i2
    !sll_int32 :: i3
    sll_int32  :: i4
    sll_int32  :: ierr
    sll_int32  :: glob_ind(4)
    sll_real64 :: alpha
    sll_int32  :: loc_sz_x1
    sll_int32  :: loc_sz_x2
    sll_int32  :: loc_sz_x3
    sll_int32  :: loc_sz_x4
    
    nc_x1 = sim%m_x1%num_cells
    nc_x2 = sim%m_x2%num_cells
    nc_x3 = sim%m_x3%num_cells
    nc_x4 = sim%m_x4%num_cells

    call sll_o_compute_local_sizes( sim%layout4d_parx1, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4 )

    if (sim%use_field_aligned_interpolation .eqv. .true.) then
      SLL_ALLOCATE(f2d(nc_x2+1,nc_x3+1),ierr)
      do i1=1,loc_sz_x1
        do i4=1,loc_sz_x4
          glob_ind(1:4) = sll_o_local_to_global( sim%layout4d_parx1, [i1,1,1,i4] )
          alpha = sim%m_x4%eta_min+real(glob_ind(4)-1,f64)*sim%m_x4%delta_eta
          f2d(1:nc_x2+1,1:nc_x3+1) = sim%f4d_parx1(i1,1:nc_x2+1,1:nc_x3+1,i4)
          call sll_s_oblic_advect_2d_constant( &
            sim%adv_x2x3, &
            sim%b_unit_x2(glob_ind(1)), &
            sim%b_unit_x3(glob_ind(1)), &
            alpha*dt, &
            f2d(1:nc_x2+1,1:nc_x3+1), &
            sim%f4d_parx1(i1,1:nc_x2+1,1:nc_x3+1,i4))
        enddo
      enddo
      SLL_DEALLOCATE_ARRAY(f2d,ierr)

    else
      SLL_ALLOCATE(f1d(nc_x3+1),ierr)
      do i2=1,loc_sz_x2
        do i1=1,loc_sz_x1
          do i4=1,loc_sz_x4
            glob_ind(1:4) = sll_o_local_to_global( sim%layout4d_parx1, [i1,i2,1,i4] )
            alpha = sim%m_x4%eta_min+real(glob_ind(4)-1,f64)*sim%m_x4%delta_eta
            f1d(1:nc_x3+1) = sim%f4d_parx1(i1,i2,1:nc_x3+1,i4)
            call sim%adv_x3%advect_1d_constant(&
              alpha*sim%b_unit_x3(glob_ind(1)), &
              dt, &
              f1d(1:nc_x3+1), &
              f1d(1:nc_x3+1))
            sim%f4d_parx1(i1,i2,1:nc_x3+1,i4) = f1d(1:nc_x3+1)
           enddo
        enddo
      enddo
      SLL_DEALLOCATE_ARRAY(f1d,ierr)

    endif

  end subroutine advection_x3
  

  subroutine advection_x4( sim, dt )
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar), &
                intent(inout) :: sim
    sll_real64, intent(in)    :: dt

    sll_real64, allocatable ::  f1d(:)
    sll_real64, allocatable ::  f1d_new(:)
    sll_int32  :: nc_x1
    sll_int32  :: nc_x2
    sll_int32  :: nc_x3
    sll_int32  :: nc_x4
    sll_int32  :: i1
    sll_int32  :: i2
    sll_int32  :: i3
    !sll_int32 :: i4
    sll_int32  :: ierr
    sll_real64 :: alpha
    sll_int32  :: loc_sz_x1
    sll_int32  :: loc_sz_x2
    sll_int32  :: loc_sz_x3
    sll_int32  :: loc_sz_x4
    
    nc_x1 = sim%m_x1%num_cells
    nc_x2 = sim%m_x2%num_cells
    nc_x3 = sim%m_x3%num_cells
    nc_x4 = sim%m_x4%num_cells

    SLL_ALLOCATE(f1d(nc_x4+1),ierr)
    SLL_ALLOCATE(f1d_new(nc_x4+1),ierr)
      
    call sll_o_compute_local_sizes( sim%layout4d_parx1, &
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
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar), &
                intent(inout) :: sim
    sll_real64, intent(in)    :: dt

    sll_real64, allocatable :: f2d(:,:)
    sll_real64, allocatable :: f2d_new(:,:)
    sll_real64, allocatable :: A2(:,:)
    sll_real64, allocatable :: v_array(:)
    sll_int32 :: nc_x1
    sll_int32 :: nc_x2
    sll_int32 :: nc_x3
    sll_int32 :: nc_x4
    sll_int32 :: i1
    sll_int32 :: i3
    sll_int32 :: i4
    sll_int32 :: ierr
    sll_int32 :: glob_ind(4)
    sll_int32 :: loc_sz_x1
    sll_int32 :: loc_sz_x2
    sll_int32 :: loc_sz_x3
    sll_int32 :: loc_sz_x4

    nc_x1 = sim%m_x1%num_cells
    nc_x2 = sim%m_x2%num_cells
    nc_x3 = sim%m_x3%num_cells
    nc_x4 = sim%m_x4%num_cells

    call sll_o_apply_remap_4d( &
      sim%remap_plan_parx1_to_parx3x4, &
      sim%f4d_parx1, &
      sim%f4d_parx3x4 )

    SLL_ALLOCATE(f2d(nc_x1+1,nc_x2+1),ierr)
    SLL_ALLOCATE(f2d_new(nc_x1+1,nc_x2+1),ierr)
    SLL_ALLOCATE(A2(nc_x1+1,nc_x2+1),ierr)
    SLL_ALLOCATE(v_array(nc_x4+1),ierr)

    call sll_o_compute_local_sizes( sim%layout4d_parx3x4, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3, &
      loc_sz_x4 )

    v_array(1:nc_x4+1) = 0._f64
    if (sim%use_field_aligned_interpolation .eqv. .false.) then
      do i4=1,nc_x4+1
        v_array(i4) = sim%m_x4%eta_min+real(i4-1,f64)*sim%m_x4%delta_eta
      enddo
    endif

    do i4 = 1,loc_sz_x4
      do i3 = 1,loc_sz_x3
        glob_ind(:) = sll_o_local_to_global( sim%layout4d_parx3x4, [1,1,i3,i4] )
        do i1=1,nc_x1+1
          A2(i1,1:nc_x2+1) = sim%A2_parx3(i1,1:nc_x2+1,i3) &
            -v_array(glob_ind(4))*sim%b_unit_x2(i1)
        enddo
        f2d(1:nc_x1+1,1:nc_x2+1) = sim%f4d_parx3x4(1:nc_x1+1,1:nc_x2+1,i3,i4) 
        call sim%adv_x1x2%advect_2d(&
          sim%A1_parx3(1:nc_x1+1,1:nc_x2+1,i3), &
          A2(1:nc_x1+1,1:nc_x2+1), &
          !sim%A2_parx3(1:nc_x1+1,1:nc_x2+1,i3), &
          dt, &
          f2d(1:nc_x1+1,1:nc_x2+1), &
          f2d_new(1:nc_x1+1,1:nc_x2+1))
        sim%f4d_parx3x4(1:nc_x1+1,1:nc_x2+1,i3,i4)=f2d_new(1:nc_x1+1,1:nc_x2+1)
      enddo  
    enddo    

    SLL_DEALLOCATE_ARRAY(f2d,ierr)
    SLL_DEALLOCATE_ARRAY(f2d_new,ierr)
    SLL_DEALLOCATE_ARRAY(A2,ierr)
    SLL_DEALLOCATE_ARRAY(v_array,ierr)

    call sll_o_apply_remap_4d( &
      sim%remap_plan_parx3x4_to_parx1, &
      sim%f4d_parx3x4, &
      sim%f4d_parx1 )
    
  end subroutine advection_x1x2

  
  subroutine delete_dk4d_field_aligned_polar( sim )
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar), &
      intent(inout) :: sim

    if (sll_f_get_collective_rank(sll_v_world_collective)==0) then
      print *,'#delete_dk4d_polar not implemented'
      print *,sim%dt
    endif
    
  end subroutine delete_dk4d_field_aligned_polar
  

  !< Initialize equilibrium profiles with analytical functions 
  ! TODO: maybes these profiles should not be part of sim object
  subroutine initialize_profiles_analytic(sim)
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar), &
      intent(inout) :: sim

    sll_int32  :: i,ierr,nc_x1
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
    
    nc_x1    = sim%m_x1%num_cells
    delta_x1 = sim%m_x1%delta_eta
    x1_min   = sim%m_x1%eta_min
    x1_max   = sim%m_x1%eta_max
    x3_min   = sim%m_x3%eta_min
    x3_max   = sim%m_x3%eta_max
    
    Lr = x1_max-x1_min
    Lz = x3_max-x3_min
    
    R0 = Lz/(2._f64*sll_p_pi)
    inv_Ln   = sim%kappan  / R0
    inv_LTi  = sim%kappaTi / R0
    inv_LTe  = sim%kappaTe / R0
    deltarn  = sim%deltarn  * Lr
    deltarTi = sim%deltarTi * Lr
    deltarTe = sim%deltarTe * Lr
    
    SLL_ALLOCATE(sim%n0_r(nc_x1+1),ierr)
    SLL_ALLOCATE(sim%Ti_r(nc_x1+1),ierr)
    SLL_ALLOCATE(sim%Te_r(nc_x1+1),ierr)
    SLL_ALLOCATE(sim%dlog_density_r(nc_x1+1),ierr)
    !SLL_ALLOCATE(sim%B_norm_r(nc_x1+1),ierr)
    !SLL_ALLOCATE(sim%Bstar_par_v_r(nc_x1+1),ierr)
    !SLL_ALLOCATE(sim%c_r(nc_x1+1),ierr)
    
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
      !sim%c_r(i) = x1*sim%iota_r(i)/R0!x1/(R0*sim%q0)
      !sim%B_norm_r(i) = (1._f64+sim%c_r(i)**2)**(sim%B_norm_exponent)
      !sim%Bstar_par_v_r = (2._f64*sim%c_r(i)-R0*sim%Dr_q0*sim%c_r(i)**2)/(1+sim%c_r(i)**2)
      !sim%Bstar_par_v_r(i) = (2._f64*sim%c_r(i)-sim%Diota_r(i)*sim%c_r(i)**2)/(1+sim%c_r(i)**2)
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
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar), &
      intent(inout) :: sim

    sll_int32 :: ierr !, itemp
    sll_int32 :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4
    sll_int32 :: power2
    sll_int32 :: power2_x3
!    sll_int32 :: k_min
 
    sim%nproc_x1 = sll_f_get_collective_size(sll_v_world_collective)
    sim%nproc_x2 = 1
    sim%nproc_x3 = 1
    sim%nproc_x4 = 1
   
    if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
    
      print *,'#num_points=',sim%m_x1%num_cells+1, &
        sim%m_x2%num_cells+1, &
        sim%m_x3%num_cells+1, &
        sim%m_x4%num_cells+1
   
      print *,'#num_proc=',sll_f_get_collective_size(sll_v_world_collective)
   
      print *,'#num_proc_parx1:',sim%nproc_x1,sim%nproc_x2,sim%nproc_x3,sim%nproc_x4
    endif
    
    power2 = int(log(real(sll_f_get_collective_size(sll_v_world_collective)))/log(2.0))
    power2_x3 = int(log(real(sim%m_x3%num_cells+1))/log(2.0))
    power2_x3 = min(power2_x3,power2)

    !--> Initialization of parallel layout of f4d in x1 direction (r)
    !-->  (x2,x3,x4) : sequential
    !-->  x1 : parallelized layout
    sim%layout4d_parx1  => sll_f_new_layout_4d( sll_v_world_collective )
    call sll_o_initialize_layout_with_distributed_array( &
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
    call sll_o_compute_local_sizes( sim%layout4d_parx1, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
 
    SLL_ALLOCATE(sim%f4d_parx1(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)
    SLL_ALLOCATE(sim%rho3d_parx1(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_parx1(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)
    SLL_ALLOCATE(sim%A3_parx1(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)

    !--> Initialization of parallel layout for various arrays in poloidal plane
    !--> x1 (r)     : parallel
    !--> x2 (theta) : sequential
    sim%layout2d_parx1  => sll_f_new_layout_2d( sll_v_world_collective )
    call sll_o_initialize_layout_with_distributed_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells, &     ! TODO: check if we should add "+1"
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%layout2d_parx1 )

    ! As above, but parallelized in theta
    sim%layout2d_parx2  => sll_f_new_layout_2d( sll_v_world_collective )
    call sll_o_initialize_layout_with_distributed_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells, &    ! TODO: check if we should add "+1" 
      sim%nproc_x2, &
      sim%nproc_x1, &
      sim%layout2d_parx2 )

    !--> Initialization of parallel layout of f4d in (x3,x4) directions
    !-->  (x3,x4) : parallelized layout
    !-->  (x1,x2) : sequential
    !--> we take the largest possible number of processors in x3
    sim%nproc_x1 = 1
    sim%nproc_x2 = 1
    sim%nproc_x3 = 2**(power2_x3)
    sim%nproc_x4 = 2**(power2-power2_x3)
     
    sim%layout4d_parx3x4  => sll_f_new_layout_4d( sll_v_world_collective )
    call sll_o_initialize_layout_with_distributed_array( &
      sim%m_x1%num_cells+1, & 
      sim%m_x2%num_cells+1, & 
      sim%m_x3%num_cells+1, &
      sim%m_x4%num_cells+1, &
      sim%nproc_x1, &
      sim%nproc_x2, &
      sim%nproc_x3, &
      sim%nproc_x4, &
      sim%layout4d_parx3x4 )
        
    ! Compute size of local arrays, then allocate them (A3 not needed)
    call sll_o_compute_local_sizes( sim%layout4d_parx3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )

    ! f4d is parallelized in (x3,x4)
    ! all other 3D arrays are only parallelized in (x3)
    ! NOTE: look into implications of this
    SLL_ALLOCATE(sim%f4d_parx3x4(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr)
    SLL_ALLOCATE(sim%rho3d_parx3(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)
    SLL_ALLOCATE(sim%phi3d_parx3(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)
    SLL_ALLOCATE(sim%A1_parx3(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)
    SLL_ALLOCATE(sim%A2_parx3(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3),ierr)
    
    !--------------------------------------------------------------------------
    ! Create objects that will perform "remap" between 2 layouts
    ! MPI communication will happen when "remap" is called
    !--------------------------------------------------------------------------

    sim%remap_plan_parx1_to_parx3x4 => sll_o_new_remap_plan( &
      sim%layout4d_parx1, &
      sim%layout4d_parx3x4, &
      sim%f4d_parx1)

    sim%remap_plan_parx3x4_to_parx1 => sll_o_new_remap_plan( &
      sim%layout4d_parx3x4, &
      sim%layout4d_parx1, &
      sim%f4d_parx3x4)

    ! Get 3D layouts from 4D layout, by discarding x4 dimension (v)
    ! If there was parallelization in v, after conversion multiple tasks will
    ! end up working on the same array chunk!
    sim%layout3d_parx1  => sll_f_new_layout_3d_from_layout_4d( sim%layout4d_parx1 )
    sim%layout3d_parx3  => sll_f_new_layout_3d_from_layout_4d( sim%layout4d_parx3x4 )
    
    sim%remap_plan_parx1_to_parx3 => sll_o_new_remap_plan( &
      sim%layout3d_parx1, &
      sim%layout3d_parx3, &
      sim%phi3d_parx1)

    !--------------------------------------------------------------------------

    ! Print new domain decomposition
    if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
      print *,'#num_proc_parx3x4:', sim%nproc_x1, sim%nproc_x2, &
                                    sim%nproc_x3, sim%nproc_x4
    endif

    !--------------------------------------------------------------------------
    ! Use MPI colors to deal with redundant data
    !--------------------------------------------------------------------------

!    k_min = sll_o_get_layout_k_min( &
!      sim%layout4d_parx3x4, &
!      sll_f_get_collective_rank(sll_v_world_collective) )
!
!    sim%new_collective_per_locx3 => sll_f_new_collective( &
!      sll_v_world_collective, &
!      k_min, &
!      sll_f_get_collective_rank(sll_v_world_collective))    
!
!    sim%new_collective_per_locx4 => sll_f_new_collective( &
!      sll_v_world_collective, &
!      sll_f_get_collective_rank(sim%new_collective_per_locx3), &
!      sll_f_get_collective_rank(sll_v_world_collective))    

  end subroutine allocate_fdistribu4d_and_QN_DK_parx1


  !----------------------------------------------------
  ! Initialization of the distribution function for
  !   drift-kinetic 4D simulation
  ! Different parallelization layouts may be used, by changing the input
  ! arguments 'layout' and 'f4d' in a consistent way
  ! TODO: use more specific input arguments instead of 'sim'
  !----------------------------------------------------
  subroutine initialize_fdistribu4d_DK(sim,layout,f4d)
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar), intent(inout) :: sim
    type(sll_t_layout_4d), pointer :: layout
    sll_real64     , pointer :: f4d(:,:,:,:)

    sll_int32  :: ierr
    sll_int32  :: i1, i2, i3, i4
    sll_int32  :: iloc1, iloc2, iloc3, iloc4
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4
    sll_int32  :: glob_ind(1:4)
    sll_real64, pointer :: x1_node(:)
    sll_real64, pointer :: x2_node(:)
    sll_real64, pointer :: x3_node(:)
    sll_real64, pointer :: x4_node(:)
    sll_real64 :: rpeak,k_x2,k_x3
    sll_real64 :: tmp_mode,tmp
    sll_real64 :: x1_min,x1_max
          
    !--> Initialization of the equilibrium distribution function
    SLL_ALLOCATE(sim%feq_x1x4(sim%m_x1%num_cells+1,sim%m_x4%num_cells+1),ierr)
    
    x1_min = sim%m_x1%eta_min
    x1_max = sim%m_x1%eta_max
        
    call sll_o_get_node_positions(sim%m_x1,x1_node)
    call sll_o_get_node_positions(sim%m_x2,x2_node)
    call sll_o_get_node_positions(sim%m_x3,x3_node)
    call sll_o_get_node_positions(sim%m_x4,x4_node)
    
    call sll_s_init_fequilibrium( &
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
    call sll_o_compute_local_sizes( layout, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
   
    k_x2  = 2._f64*sll_p_pi/(sim%m_x2%eta_max - sim%m_x2%eta_min)
    k_x3  = 2._f64*sll_p_pi/(sim%m_x3%eta_max - sim%m_x3%eta_min)
      
    rpeak = x1_min+sim%rho_peak*(x1_max-x1_min) 
    
    do iloc4 = 1,loc4d_sz_x4
      do iloc3 = 1,loc4d_sz_x3
        do iloc2 = 1,loc4d_sz_x2
          do iloc1 = 1,loc4d_sz_x1
            glob_ind(:) = sll_o_local_to_global( layout, [iloc1,iloc2,iloc3,iloc4] )
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
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar), &
                intent(in) :: sim
    sll_real64, intent(in) :: x1, x4

    sll_real64 :: compute_equil_analytic
    sll_real64 :: tmp(2), rpeak, x1_min, x1_max

    x1_min = sim%m_x1%eta_min
    x1_max = sim%m_x1%eta_max
    rpeak  = x1_min+sim%rho_peak*(x1_max-x1_min)
    tmp(1) = sim%n0_at_rpeak*exp(-sim%kappan*sim%deltarn*tanh((x1-rpeak)/sim%deltarn))
    tmp(2) = exp(-sim%kappaTi*sim%deltarTi*tanh((x1-rpeak)/sim%deltarTi))
    compute_equil_analytic = tmp(1)/sqrt(2._f64*sll_p_pi*tmp(2))*exp(-0.5_f64*x4**2/tmp(2))
  
  end function compute_equil_analytic


!------------------------------------------------------------------------------
! Quasi-neutral solver
!------------------------------------------------------------------------------
  subroutine solve_quasi_neutral_parx1(sim)
    class(sll_t_simulation_4d_drift_kinetic_field_aligned_polar), intent(inout) :: sim

    sll_int32 :: glob_ind(4)
    sll_int32 :: loc4d_sz_x1
    sll_int32 :: loc4d_sz_x2
    sll_int32 :: loc4d_sz_x3
    sll_int32 :: loc4d_sz_x4
    sll_int32 :: i
    sll_real64 :: n0

    select case (sim%QN_case)
      case (SLL_NO_QUASI_NEUTRAL)
        call sll_o_compute_local_sizes( &
          sim%layout4d_parx1, &
          loc4d_sz_x1, &
          loc4d_sz_x2, &
          loc4d_sz_x3, &
          loc4d_sz_x4 )        
        sim%phi3d_parx1(:,1:loc4d_sz_x2,:) = 0._f64
      case (SLL_QUASI_NEUTRAL_WITHOUT_ZONAL_FLOW)
        call sll_o_compute_local_sizes( &
          sim%layout4d_parx1, &
          loc4d_sz_x1, &
          loc4d_sz_x2, &
          loc4d_sz_x3, &
          loc4d_sz_x4 )        

        ! Compute RHS for QNP solver (divide by n0 and subtract 1)
        ! Store result in place!! (rho3d_parx1)
        do i=1,loc4d_sz_x1
          glob_ind(:) = sll_o_local_to_global( sim%layout4d_parx1, [i,1,1,1] )
          n0 = sim%n0_r(glob_ind(1))
          sim%rho3d_parx1(i,:,:) = sim%rho3d_parx1(i,:,:) / n0 - 1._f64
        enddo
 
        ! Solver
        call sim%poisson3d%compute_phi_from_rho( &
          sim%phi3d_parx1(:,1:loc4d_sz_x2-1,:), &
          sim%rho3d_parx1(:,1:loc4d_sz_x2-1,:) )
 
        ! Enforce periodic boundary conditions in theta
        sim%phi3d_parx1(:,loc4d_sz_x2,:) = sim%phi3d_parx1(:,1,:)  
       
        ! Copy phi3d_parx1 into phi3d_parx3
        call sll_o_apply_remap_3d( &
          sim%remap_plan_parx1_to_parx3, &
          sim%phi3d_parx1, &
          sim%phi3d_parx3 )           
!
!      case default
!        print *,'#bad choice for QN_case', sim%QN_case
!        print *,'#in init_dk4d_polar'
!        stop

    end select

  end subroutine solve_quasi_neutral_parx1

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine write_mesh_x1x2_polar( mesh_x1, mesh_x2, dataset_x1, dataset_x2 )
    type(sll_t_cartesian_mesh_1d), intent(in   ) :: mesh_x1
    type(sll_t_cartesian_mesh_1d), intent(in   ) :: mesh_x2
    character(len=*)           , intent(  out) :: dataset_x1
    character(len=*)           , intent(  out) :: dataset_x2

    sll_real64, allocatable :: x1(:,:), x2(:,:)
    sll_real64              :: r, dr, rmin, rmax, theta, dtheta
    sll_int32               :: i, j, nnodes_x1, nnodes_x2
    sll_int32               :: error
    type(sll_t_hdf5_serial) :: hdf5_file

    nnodes_x1 = mesh_x1%num_cells+1
    nnodes_x2 = mesh_x2%num_cells+1
    rmin      = mesh_x1%eta_min
    rmax      = mesh_x1%eta_max
    dr        = mesh_x1%delta_eta
    dtheta    = mesh_x2%delta_eta

    SLL_ALLOCATE( x1(nnodes_x1,nnodes_x2), error )
    SLL_ALLOCATE( x2(nnodes_x1,nnodes_x2), error )

    do j = 1,nnodes_x2
      do i = 1,nnodes_x1
        r       = rmin+real(i-1,f64)*dr
        theta   = real(j-1,f64)*dtheta
        x1(i,j) = r*cos(theta)
        x2(i,j) = r*sin(theta)
      end do
    end do

    call hdf5_file%init( "mesh_x1x2_polar.h5" )
    call hdf5_file%create()
    call hdf5_file%write_array( x1, "/x1" )
    call hdf5_file%write_array( x2, "/x2" )
    call hdf5_file%delete()

    dataset_x1 = "mesh_x1x2_polar.h5:/x1"
    dataset_x2 = "mesh_x1x2_polar.h5:/x2"

    deallocate( x1 )
    deallocate( x2 )

  end subroutine


  subroutine write_mesh_x1x2_polar_light( mesh_x1, mesh_x2, dataset_x1, dataset_x2 )
    type(sll_t_cartesian_mesh_1d), intent(in   ) :: mesh_x1
    type(sll_t_cartesian_mesh_1d), intent(in   ) :: mesh_x2
    character(len=*)           , intent(  out) :: dataset_x1
    character(len=*)           , intent(  out) :: dataset_x2

    !sll_real64, allocatable :: x1(:,:), x2(:,:)
    !sll_real64              :: r, dr, rmin, rmax, theta, dtheta
    !sll_int32               :: i, j, nnodes_x1, nnodes_x2
    !sll_int32               :: file_id, error
    !type(sll_t_hdf5_serial) :: hdf5_file

    !nnodes_x1 = mesh_x1%num_cells+1
    !nnodes_x2 = mesh_x2%num_cells+1
    !rmin      = mesh_x1%eta_min
    !rmax      = mesh_x1%eta_max
    !dr        = mesh_x1%delta_eta
    !dtheta    = mesh_x2%delta_eta

    !SLL_ALLOCATE( x1(nnodes_x1,nnodes_x2), error )
    !SLL_ALLOCATE( x2(nnodes_x1,nnodes_x2), error )

    !do j = 1,nnodes_x2
    !  do i = 1,nnodes_x1
    !    r       = rmin+real(i-1,f64)*dr
    !    theta   = real(j-1,f64)*dtheta
    !    x1(i,j) = r*cos(theta)
    !    x2(i,j) = r*sin(theta)
    !  end do
    !end do

    !call hdf5_file%init( "mesh_x1x2_polar.h5" )
    !call hdf5_file%create()
    !call hdf5_file%write_array( x1, "/x1" )
    !call hdf5_file%write_array( x2, "/x2" )
    !call hdf5_file%delete()

    dataset_x1 = "mesh_x1x2_polar.h5:/x1"
    dataset_x2 = "mesh_x1x2_polar.h5:/x2"

    !deallocate( x1 )
    !deallocate( x2 )

  end subroutine


  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine write_mesh_x2x3_cart( mesh_x2, mesh_x3, dataset_x2, dataset_x3 )
    type(sll_t_cartesian_mesh_1d), intent(in   ) :: mesh_x2
    type(sll_t_cartesian_mesh_1d), intent(in   ) :: mesh_x3
    character(len=*)           , intent(  out) :: dataset_x2
    character(len=*)           , intent(  out) :: dataset_x3

    sll_real64, allocatable :: x2(:,:), x3(:,:)
    sll_int32               :: i, j, nnodes_x2, nnodes_x3
    sll_int32               :: file_id, error

    nnodes_x2 = mesh_x2%num_cells+1
    nnodes_x3 = mesh_x3%num_cells+1

    SLL_ALLOCATE( x2(nnodes_x2,nnodes_x3), error )
    SLL_ALLOCATE( x3(nnodes_x2,nnodes_x3), error )

    do i = 1,nnodes_x2
      x2(i,:) = mesh_x2%eta1_node( i )
    end do

    do j = 1,nnodes_x3
      x3(:,j) = mesh_x3%eta1_node( j )
    end do

    call sll_o_hdf5_file_create( "mesh_x2x3_cart.h5", file_id, error )
    call sll_o_hdf5_write_array( file_id, x2, "/x2", error )
    call sll_o_hdf5_write_array( file_id, x3, "/x3", error )
    call sll_o_hdf5_file_close ( file_id, error )

    dataset_x2 = "mesh_x2x3_cart.h5:/x2"
    dataset_x3 = "mesh_x2x3_cart.h5:/x3"

    deallocate( x2 )
    deallocate( x3 )

  end subroutine



  subroutine write_mesh_x2x3_cart_light( mesh_x2, mesh_x3, dataset_x2, dataset_x3 )
    type(sll_t_cartesian_mesh_1d), intent(in   ) :: mesh_x2
    type(sll_t_cartesian_mesh_1d), intent(in   ) :: mesh_x3
    character(len=*)           , intent(  out) :: dataset_x2
    character(len=*)           , intent(  out) :: dataset_x3

    !sll_real64, allocatable :: x2(:,:), x3(:,:)
    !sll_int32               :: i, j, nnodes_x2, nnodes_x3
    !sll_int32               :: file_id, error

    !nnodes_x2 = mesh_x2%num_cells+1
    !nnodes_x3 = mesh_x3%num_cells+1

    !SLL_ALLOCATE( x2(nnodes_x2,nnodes_x3), error )
    !SLL_ALLOCATE( x3(nnodes_x2,nnodes_x3), error )

    !do i = 1,nnodes_x2
    !  x2(i,:) = mesh_x2%eta1_node( i )
    !end do

    !do j = 1,nnodes_x3
    !  x3(:,j) = mesh_x3%eta1_node( j )
    !end do

    !call sll_o_hdf5_file_create( "mesh_x2x3_cart.h5", file_id, error )
    !call sll_o_hdf5_write_array( file_id, x2, "/x2", error )
    !call sll_o_hdf5_write_array( file_id, x3, "/x3", error )
    !call sll_o_hdf5_file_close ( file_id, error )

    dataset_x2 = "mesh_x2x3_cart.h5:/x2"
    dataset_x3 = "mesh_x2x3_cart.h5:/x3"

    !deallocate( x2 )
    !deallocate( x3 )

  end subroutine


!==============================================================================

end module sll_m_sim_bsl_dk_3d1v_polar_field_aligned





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
!    sll_int32, inent(in) :: num_points_x1
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
!    do           i=1,num_points_x1    
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



