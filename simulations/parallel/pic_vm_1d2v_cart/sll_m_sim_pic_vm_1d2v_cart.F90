! Simulation of 1d2v Vlasov-Maxwell with simple PIC method, periodic boundary conditions, Weibel instability. FEM with splines, degree 3 for B and 2 for E

! author: Katharina Kormann, IPP

module sll_m_sim_pic_vm_1d2v_cart

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"


  use sll_m_ascii_io, only: &
       sll_s_ascii_file_create

  use sll_m_binomial_filter, only: &
       sll_t_binomial_filter

  use sll_m_collective, only: &
       sll_o_collective_allreduce, &
       sll_s_collective_reduce_real64, &
       sll_f_get_collective_rank, &
       sll_f_get_collective_size, &
       sll_v_world_collective

  use sll_m_constants, only: &
       sll_p_twopi

  use sll_m_control_variate, only: &
       sll_t_control_variate, &
       sll_t_control_variates

  use sll_m_filter_base_1d, only: &
       sll_c_filter_base_1d

  use sll_m_fft_filter_1d, only: &
       sll_t_fft_filter_1d

  use sll_m_gauss_legendre_integration, only: &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_initial_distribution, only : &
       sll_c_distribution_params, &
       sll_s_initial_distribution_new

  use sll_m_io_utilities, only : &
       sll_s_read_data_real_array, &
       sll_s_concatenate_filename_and_path

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_maxwell_1d_base, only: &
       sll_c_maxwell_1d_base

  use sll_m_maxwell_1d_fem, only: &
       sll_t_maxwell_1d_fem

  use sll_m_maxwell_1d_fem_sm, only: &
       sll_t_maxwell_1d_fem_sm

  use sll_m_maxwell_1d_trafo, only:&
       sll_t_maxwell_1d_trafo

  use sll_m_maxwell_1d_ps, only: &
       sll_t_maxwell_1d_ps

  use sll_m_maxwell_clamped_1d_fem_sm, only: &
       sll_t_maxwell_clamped_1d_fem_sm

  use sll_m_maxwell_clamped_1d_trafo, only:&
       sll_t_maxwell_clamped_1d_trafo

  use sll_mpi, only: &
       mpi_sum

  use sll_m_particle_mesh_coupling_base_1d, only: &
       sll_p_galerkin, &
       sll_c_particle_mesh_coupling_1d

  use sll_m_particle_mesh_coupling_spline_1d, only: &
       sll_t_particle_mesh_coupling_spline_1d, &
       sll_s_new_particle_mesh_coupling_spline_1d, &
       sll_s_new_particle_mesh_coupling_spline_1d_ptr

  use sll_m_particle_mesh_coupling_spline_cl_1d, only: &
       sll_t_particle_mesh_coupling_spline_cl_1d, &
       sll_s_new_particle_mesh_coupling_spline_cl_1d, &
       sll_s_new_particle_mesh_coupling_spline_cl_1d_ptr

  use sll_m_particle_mesh_coupling_spline_strong_1d, only: &
       sll_t_particle_mesh_coupling_spline_strong_1d, &
       sll_s_new_particle_mesh_coupling_spline_strong_1d

  use sll_m_particle_mesh_coupling_spline_smooth_1d, only: &
       sll_t_particle_mesh_coupling_spline_smooth_1d, &
       sll_s_new_particle_mesh_coupling_spline_smooth_1d, &
       sll_s_new_particle_mesh_coupling_spline_smooth_1d_ptr

  use sll_m_particle_group_1d2v, only: &
       sll_t_particle_group_1d2v

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base, &
       sll_t_particle_array

  use sll_m_particle_sampling, only: &
       sll_t_particle_sampling

  use sll_m_sim_base, only: &
       sll_c_simulation_base_class

  use sll_m_splines_pp, only: &
       sll_t_spline_pp_1d, &
       sll_s_spline_pp_init_1d, &
       sll_s_spline_pp_free_1d, &
       sll_f_spline_pp_horner_1d, &
       sll_p_boundary_periodic, &
       sll_p_boundary_clamped, &
       sll_p_boundary_clamped_square, &
       sll_p_boundary_clamped_cubic, &
       sll_p_boundary_clampeddiri

  use sll_m_time_propagator_base, only: &
       sll_c_time_propagator_base

  use sll_m_time_propagator_pic_vm_1d2v_hs, only: &
       sll_s_new_time_propagator_pic_vm_1d2v_hs, &
       sll_t_time_propagator_pic_vm_1d2v_hs

  use sll_m_time_propagator_pic_vm_1d2v_momentum, only: &
       sll_s_new_time_propagator_pic_vm_1d2v_momentum, &
       sll_t_time_propagator_pic_vm_1d2v_momentum

  use sll_m_time_propagator_pic_vm_1d2v_boris, only: &
       sll_t_time_propagator_pic_vm_1d2v_boris

  use sll_m_time_propagator_pic_vm_1d2v_disgradE, only: &
       sll_t_time_propagator_pic_vm_1d2v_disgradE

  use sll_m_time_propagator_pic_vm_1d2v_disgradEC, only: &
       sll_t_time_propagator_pic_vm_1d2v_disgradEC

  use sll_m_time_propagator_pic_vm_1d2v_disgradEC_sub, only: &
       sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub

  use sll_m_time_propagator_pic_vm_1d2v_hs_trafo, only: &
       sll_t_time_propagator_pic_vm_1d2v_hs_trafo

  use sll_m_time_propagator_pic_vm_1d2v_disgradE_trafo, only: &
       sll_t_time_propagator_pic_vm_1d2v_disgradE_trafo

  use sll_m_time_propagator_pic_vm_1d2v_disgradEC_trafo, only: &
       sll_t_time_propagator_pic_vm_1d2v_disgradEC_trafo

  use sll_m_time_propagator_pic_vm_3d3v_cl_helper, only: &
       sll_p_boundary_particles_periodic, &
       sll_p_boundary_particles_singular, &
       sll_p_boundary_particles_reflection, &
       sll_p_boundary_particles_absorption

  use sll_m_time_propagator_pic_vm_1d2v_cef, only: &
       sll_t_time_propagator_pic_vm_1d2v_cef

  use sll_m_time_propagator_pic_vm_1d2v_ecsim, only: &
       sll_t_time_propagator_pic_vm_1d2v_ecsim

  use sll_m_time_propagator_pic_vm_1d2v_ecsim2o, only: &
       sll_t_time_propagator_pic_vm_1d2v_ecsim2o

  use sll_m_time_propagator_pic_vm_1d2v_subcyc, only: &
       sll_t_time_propagator_pic_vm_1d2v_subcyc

  use sll_m_time_propagator_pic_vm_1d2v_zigsub, only: &
       sll_t_time_propagator_pic_vm_1d2v_zigsub

  use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_t_time_mark

  use sll_m_utilities, only : &
       sll_s_int2string


  implicit none

  public :: &
       sll_t_sim_pic_vm_1d2v_cart

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: sll_p_splitting_hs = 0
  sll_int32, parameter :: sll_p_splitting_boris = 1
  sll_int32, parameter :: sll_p_splitting_disgradE = 2
  sll_int32, parameter :: sll_p_splitting_cef = 4
  sll_int32, parameter :: sll_p_splitting_ecsim = 6
  sll_int32, parameter :: sll_p_splitting_ecsim2o =7 
  sll_int32, parameter :: sll_p_splitting_disgradEC = 10
  sll_int32, parameter :: sll_p_splitting_disgradEC_sub = 13
  sll_int32, parameter :: sll_p_splitting_hs_trafo = 14
  sll_int32, parameter :: sll_p_splitting_disgradE_trafo = 15
  sll_int32, parameter :: sll_p_splitting_disgradEC_trafo = 16
  sll_int32, parameter :: sll_p_splitting_subcyc = 17
sll_int32, parameter :: sll_p_splitting_zigsub = 19 !< subcycling with zig zagging
  sll_int32, parameter :: sll_p_splitting_momentum = 20


  sll_int32, parameter :: sll_p_onegaussian = 0
  sll_int32, parameter :: sll_p_twogaussian = 1

  sll_int32, parameter :: sll_p_bfield_cos = 0
  sll_int32, parameter :: sll_p_bfield_sin = 1
  sll_int32, parameter :: sll_p_bfield_constant = 2
  sll_int32, parameter :: sll_p_bfield_cos_const = 3

  sll_int32, parameter :: sll_p_strang_splitting=0
  sll_int32, parameter :: sll_p_splitting_fourth=1
  sll_int32, parameter :: sll_p_lie_splitting=2
  sll_int32, parameter :: sll_p_lie_splitting_back=3
  sll_int32, parameter :: sll_p_splitting_fourth_10steps=4
  sll_int32, parameter :: sll_p_splitting_second_4steps=5

  type, extends(sll_c_simulation_base_class) :: sll_t_sim_pic_vm_1d2v_cart

     ! Abstract particle group
     class(sll_t_particle_array), pointer :: particle_group !< Particle group

     ! Fields
     sll_real64, pointer :: phi_dofs(:) !< DoFs describing the scalar potential
     sll_real64, pointer :: efield_dofs(:,:) !< DoFs describing the two components of the electric field
     sll_real64, pointer :: efield_dofs_n(:,:) !< DoFs describing the two components of the electric field
     sll_real64, allocatable :: bfield_dofs(:) !< DoFs describing the third component of the magnetic field
     sll_real64, allocatable :: rhob(:) !< charge at the boundary


     ! Cartesian mesh
     sll_real64 :: delta_x !< Grid spacing

     sll_real64 :: time !< time

     ! Maxwell solver 
     ! Abstract 
     class(sll_c_maxwell_1d_base), allocatable :: maxwell_solver !< Maxwell solver

     ! Abstract kernel smoothers
     class(sll_c_particle_mesh_coupling_1d), allocatable :: kernel_smoother_0 !< Particle mesh coupling
     class(sll_c_particle_mesh_coupling_1d), allocatable :: kernel_smoother_1 !< Particle mesh coupling   


     ! Specific operator splitting
     class(sll_c_time_propagator_base), allocatable :: propagator !< time propagator object 
     sll_int32 :: splitting_case !< time splitting case
     sll_int32 :: splitting_type !< time splitting type

     ! Control variate
     type(sll_t_control_variates) :: control_variate !< control variate
     sll_int32  :: no_weights !< weight number

     ! Physical parameters
     class(sll_c_distribution_params), allocatable :: init_distrib_params !< distribution parameter
     sll_real64 :: beta !< magnitude of initial magnetic field
     sll_real64 :: domain(3) !< domain length: x_min, x_max, Lx
     type(sll_t_particle_sampling) :: sampler !< particle sampler
     sll_real64 :: plasma_betar(3) !< reciprocal of plasma beta
     sll_real64 :: force_sign !< sign of particle force
     logical :: adiabatic_electrons = .false. !< true for run with adiabatic electrons

     ! Simulation parameters
     sll_real64 :: delta_t !< time step
     sll_int32  :: n_time_steps !< amount of time steps
     sll_int32  :: n_particles !< amount of particles per mpi process
     sll_int32  :: n_total_particles !< total number of particles
     sll_int32  :: degree_smoother !< spline degree
     sll_int32  :: degree_fem !< spline degree
     sll_int32  :: n_gcells !< numer of gridcells
     sll_int32  :: n_dofs0 !< number of Dofs for 0-form
     sll_int32  :: n_dofs1 !< number of Dofs for 1-form
     logical    :: boundary = .false. !< true for non periodic field boundary
     sll_int32  :: boundary_fields = 100 !< field boundary conditions
     sll_int32  :: boundary_particles = 100 !< particle boundary conditions

     ! Parameters for MPI
     sll_int32  :: rank !< mpi rank
     sll_int32  :: world_size !< mpi world size

     ! Case definitions
     sll_int32  :: initial_bfield !< case for intial magnetic field
     sll_real64 :: charge !< charge of particle species

     sll_real64, allocatable :: background_charge(:) !< background distribution of neutralizing ions

     ! Output
     character(len=256)   :: file_prefix !< name of diagnostic file
     logical              :: output_fields !< logical for print out fields
     logical              :: output_particles !< logical for print out particles

     ! For ctest
     logical    :: ctest_passed = .false. !< logical for ctest
     logical    :: make_ctest = .false. !< logical for ctest
     character(len=256)   :: ctest_ref_file_rho, ctest_ref_file_thdiag !< names of ctest reference files

     !coordinate transformation
     type(sll_t_mapping_3d) :: map !< coordinate transformation
     logical                :: ct = .false. !< true if coordinate transformation is used

     !spline_pp
     type(sll_t_spline_pp_1d) :: spline0_pp !< pp spline

     logical :: strong_ampere !< true for strong ampere formulation

     ! Filter
     class(sll_c_filter_base_1d), allocatable :: filter !< filter
     type(sll_t_binomial_filter) :: bfilter !< binomial filter

   contains
     procedure :: init_from_file => init_pic_vm_1d2v !< Initialization
     procedure :: run => run_pic_vm_1d2v !< Simulation
     procedure :: delete => delete_pic_vm_1d2v !< Finalization

  end type sll_t_sim_pic_vm_1d2v_cart


contains


  !------------------------------------------------------------------------------!
  ! Read in the simulation parameters from input file
  subroutine init_pic_vm_1d2v (sim, filename)
    class(sll_t_sim_pic_vm_1d2v_cart), intent(inout) :: sim !< Singlespecies simulation
    character(len=*),                  intent(in)    :: filename !< filename
    !local variables
    sll_int32   :: io_stat
    sll_int32   :: input_file, file_id 
    sll_int32   :: ierr
    type(sll_t_time_mark) :: start_init, end_init
    sll_real64         :: delta_t = 0._f64
    sll_int32          :: n_time_steps = 0
    sll_real64         :: beta = 0.0_f64
    character(len=256) :: initial_distrib = "none"
    character(len=256) :: initial_bfield = "none"
    sll_real64         :: charge = -1._f64
    sll_real64         :: mass = 1._f64
    sll_real64         :: plasma_beta(3)  = 1._f64
    character(len=256) :: particle_force = "repulsion"
    logical            :: electrostatic = .false.
    character(len=256) :: file_prefix = "none"
    logical            :: output_fields = .false.
    logical            :: output_particles = .false.
    sll_int32          :: ng_x = 0
    sll_real64         :: x1_min = 0._f64, x1_max = 0._f64
    logical            :: jmean = .false.
    sll_int32          :: degree_fem = -2
    sll_int32          :: n_particles = 0
    character(len=256) :: sampling_case = "none"
    character(len=256) :: splitting_case = "none"
    sll_int32          :: spline_degree = 0
    character(len=256) :: splitting_type = "none"
    logical            :: strong_ampere = .false.
    logical            :: eval_grid_points
    character(len=256) :: boundary_fields = "none"
    character(len=256) :: boundary_particles = "none"
    logical            :: smoothing = .false.
    character(len=256) :: ctest_case = "none"
    character(len=256) :: filtering = "none"
    sll_int32          :: filter_iter = 0
    sll_int32          :: mode = 2
    logical            :: with_control_variate = .false.
    sll_int32          :: n_sub_iter = 0


    namelist /sim_params/         delta_t, n_time_steps, beta, initial_distrib, initial_bfield, charge, mass, plasma_beta, particle_force, electrostatic

    namelist /output/             file_prefix, output_fields, output_particles

    namelist /grid_dims/          ng_x, x1_min, x1_max, jmean, degree_fem

    namelist /pic_params/         n_particles, sampling_case, splitting_case, spline_degree, splitting_type, boundary_fields, boundary_particles, smoothing, filtering, filter_iter, mode, with_control_variate, strong_ampere

    namelist /ctest/              ctest_case
    namelist /time_iterate/       n_sub_iter

    call sll_s_set_time_mark( start_init )

    ! Read parameters from file
    open(newunit = input_file, file=trim(filename), status='old', IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'init_pic_1d2v() failed to open file ', filename
       STOP
    end if

    read(input_file, sim_params)
    call sll_s_initial_distribution_new( trim(initial_distrib), [1,2], input_file, sim%init_distrib_params )
    read(input_file, output)
    read(input_file, grid_dims)
    read(input_file, pic_params)
    read(input_file, ctest )
    close (input_file)

    ! Set MPI parameters
    sim%world_size = sll_f_get_collective_size(sll_v_world_collective)
    sim%rank = sll_f_get_collective_rank(sll_v_world_collective)

    ! Copy the read parameters into the simulation parameters
    sim%delta_t = delta_t
    sim%n_time_steps = n_time_steps
    sim%beta = beta
    sim%charge = charge

    sim%plasma_betar = 1._f64/plasma_beta

    select case( particle_force)
    case( "attraction" )
       sim%force_sign = -1._f64
    case( "repulsion" )
       sim%force_sign = 1._f64
       if( sim%charge > 0._f64) sim%adiabatic_electrons = .true.
    case default
       sim%force_sign = 1._f64
       if( sim%charge > 0._f64) sim%adiabatic_electrons = .true.
    end select

    select case ( initial_bfield )
    case( "cos")
       sim%initial_bfield = sll_p_bfield_cos
    case( "sin" )
       sim%initial_bfield = sll_p_bfield_sin
    case( "constant" )
       sim%initial_bfield = sll_p_bfield_constant
    case( "cos_constant" )
       sim%initial_bfield = sll_p_bfield_cos_const
    case default
       print*, '#initial bfield must be either sin or cos or constant.'
    end select

    ! Output  
    sim%file_prefix = file_prefix
    sim%output_fields = output_fields
    sim%output_particles = output_particles

    sim%n_gcells = ng_x
    sim%domain = [x1_min, x1_max, x1_max - x1_min ]
    sim%delta_x = (x1_max - x1_min)/real(ng_x, f64)

    sim%n_particles = n_particles/sim%world_size
    sim%degree_smoother = spline_degree

    sim%strong_ampere = strong_ampere

    if ( sim%strong_ampere .eqv. .false. ) then
       sim%degree_fem = sim%degree_smoother ! Both degrees need to be synchronized for weak Ampere (since there actually is no smoother )
    else
       if ( degree_fem > -2 ) then
          sim%degree_fem = degree_fem
       else
          sim%degree_fem = sim%degree_smoother ! Default is same as zero form if not given
       end if
    end if

    ! filter
    select case( filtering )
    case( "binomial" )
       allocate( sll_t_binomial_filter :: sim%filter )
    case( "fft" )
       allocate( sll_t_fft_filter_1d :: sim%filter )
    case default
       allocate( sll_t_fft_filter_1d :: sim%filter )
    end select
    call sim%filter%init ( filter_iter, sim%n_gcells, mode )

    call sim%bfilter%init ( filter_iter, sim%n_gcells )

    call sim%sampler%init( trim(sampling_case), [1,2], sim%n_particles, sim%rank)
    sim%n_total_particles = sim%n_particles * sim%world_size

    if (with_control_variate .eqv. .true.) then
       sim%no_weights = 3
    else
       sim%no_weights = 1
    end if

    !boundary conditions
    select case(boundary_fields)
    case("clamped")
       sim%boundary = .true.
       if( sim%degree_smoother == 2 )then
          sim%boundary_fields=sll_p_boundary_clamped_square
       else if( sim%degree_smoother == 3 )then
          sim%boundary_fields=sll_p_boundary_clamped_cubic
       else
          sim%boundary_fields=sll_p_boundary_clamped
       end if
       sim%n_dofs0= ng_x+spline_degree
       sim%n_dofs1= ng_x+spline_degree-1
       call sll_s_spline_pp_init_1d( sim%spline0_pp, spline_degree, ng_x, sim%boundary_fields)
    case("clampeddiri")
       sim%boundary = .true.
       sim%boundary_fields=sll_p_boundary_clampeddiri
       sim%n_dofs0= ng_x+spline_degree
       sim%n_dofs1= ng_x+spline_degree-1
       call sll_s_spline_pp_init_1d( sim%spline0_pp, spline_degree, ng_x, sim%boundary_fields)
    case("periodic")
       sim%boundary = .true.
       sim%boundary_fields=sll_p_boundary_periodic
       sim%n_dofs0= ng_x+spline_degree
       sim%n_dofs1= ng_x+spline_degree-1
       call sll_s_spline_pp_init_1d( sim%spline0_pp, spline_degree, ng_x, sim%boundary_fields)
    case default
       sim%boundary = .false.
       sim%n_dofs0= ng_x
       sim%n_dofs1= ng_x
    end select

    select case(boundary_particles)
    case( "periodic" )
       sim%boundary_particles=sll_p_boundary_particles_periodic
    case( "singular" )
       sim%boundary_particles=sll_p_boundary_particles_singular
    case( "reflection" )
       sim%boundary_particles=sll_p_boundary_particles_reflection
    case( "absorption" )
       sim%boundary_particles=sll_p_boundary_particles_absorption
    end select

    select case(splitting_case)
    case("splitting_hs")
       sim%splitting_case = sll_p_splitting_hs
    case("splitting_boris")
       sim%splitting_case = sll_p_splitting_boris
    case("splitting_disgradE")
       sim%splitting_case = sll_p_splitting_disgradE
    case("splitting_cef")
       sim%splitting_case = sll_p_splitting_cef
    case("splitting_ecsim")
       sim%splitting_case = sll_p_splitting_ecsim
    case("splitting_ecsim2o")
       sim%splitting_case = sll_p_splitting_ecsim2o
    case("splitting_disgradEC")
       sim%splitting_case = sll_p_splitting_disgradEC
    case("splitting_disgradEC_sub")
       sim%splitting_case = sll_p_splitting_disgradEC_sub
    case("splitting_subcyc")
       sim%splitting_case = sll_p_splitting_subcyc
    case("splitting_zigsub")
       sim%splitting_case = sll_p_splitting_zigsub
    case("splitting_hs_trafo")
       sim%splitting_case = sll_p_splitting_hs_trafo
       sim%ct=.true.
       sim%domain = [0._f64, 1._f64, x1_max - x1_min ]
       sim%delta_x = 1._f64/real(ng_x, f64)
       call sim%map%init_from_file(filename)
    case("splitting_disgradE_trafo")
       sim%splitting_case = sll_p_splitting_disgradE_trafo
       sim%ct=.true.
       sim%domain = [0._f64, 1._f64, x1_max - x1_min ]
       sim%delta_x = 1._f64/real(ng_x, f64)
       call sim%map%init_from_file(filename)
    case("splitting_disgradEC_trafo")
       sim%splitting_case = sll_p_splitting_disgradEC_trafo
       sim%ct=.true.
       sim%domain = [0._f64, 1._f64, x1_max - x1_min ]
       sim%delta_x = 1._f64/real(ng_x, f64)
       call sim%map%init_from_file(filename)
    case ("splitting_momentum")
       sim%splitting_case = sll_p_splitting_momentum
    case default
       print*, '#splitting case ', splitting_case, ' not implemented.'
    end select

    select case(splitting_type)
    case("strang")
       sim%splitting_type=sll_p_strang_splitting
    case("strang_fourth")
       sim%splitting_type=sll_p_splitting_fourth
    case("lie")
       sim%splitting_type=sll_p_lie_splitting
    case("lie_back")
       sim%splitting_type=sll_p_lie_splitting_back
    case("fourth_10steps")
       sim%splitting_type=sll_p_splitting_fourth_10steps
    case("second_4steps")
       sim%splitting_type=sll_p_splitting_second_4steps
    case default
       sim%splitting_type=sll_p_strang_splitting
    end select

    ! Ctest
    select case( ctest_case)
    case ("hs")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_hs_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_hs_pic_vm_1d2v_cart_thdiag.dat"
    case ("hs_stronga")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_hs_stronga_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_hs_stronga_pic_vm_1d2v_cart_thdiag.dat"
    case("cef")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_cef_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_cef_pic_vm_1d2v_cart_thdiag.dat"
    case("ecsim")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_ecsim_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_ecsim_pic_vm_1d2v_cart_thdiag.dat"
    case("disgradE")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_disgradE_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_disgradE_pic_vm_1d2v_cart_thdiag.dat"
    case("disgradEC")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_disgradEC_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_disgradEC_pic_vm_1d2v_cart_thdiag.dat"
    case("disgradEC_stronga")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_disgradEC_stronga_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_disgradEC_stronga_pic_vm_1d2v_cart_thdiag.dat"
    case("disgradEC_sub")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_disgradEC_sub_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_disgradEC_sub_pic_vm_1d2v_cart_thdiag.dat"
    case("hs_smooth")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_hs_smooth_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_hs_smooth_pic_vm_1d2v_cart_thdiag.dat"
    case("disgradE_smooth")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_disgradE_smooth_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_disgradE_smooth_pic_vm_1d2v_cart_thdiag.dat"
    case("disgradEC_smooth")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_disgradEC_smooth_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_disgradEC_smooth_pic_vm_1d2v_cart_thdiag.dat"
    case("disgradEC_sub_smooth")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_disgradEC_sub_smooth_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_disgradEC_sub_smooth_pic_vm_1d2v_cart_thdiag.dat"
    case("disgradE_trafo")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_disgradE_trafo_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_disgradE_trafo_pic_vm_1d2v_cart_thdiag.dat"
    case("disgradEC_trafo")
       sim%make_ctest = .true.
       sim%ctest_ref_file_rho = "ctests/reffile_disgradEC_trafo_pic_vm_1d2v_cart_rho.dat"
       sim%ctest_ref_file_thdiag = "ctests/reffile_disgradEC_trafo_pic_vm_1d2v_cart_thdiag.dat"
    end select


    ! Initialize the particles   (mass set to 1.0 and charge set to -1.0)
    allocate( sim%particle_group )
    sim%particle_group%n_species = 1
    allocate( sll_t_particle_group_1d2v :: sim%particle_group%group(sim%particle_group%n_species) )
    select type ( qp => sim%particle_group%group(1) )
    type is (  sll_t_particle_group_1d2v )
       ! Note: This call produces a segmentation fault with the INTEL 17 compiler
       ! Therefore we manually initialize here
       ! TODO: Fix the problem with the init function
       call qp%init(sim%n_particles, &
            sim%n_total_particles, sim%charge, mass, sim%no_weights)
       !qp%n_particles = sim%n_particles
       !qp%n_total_particles = sim%n_total_particles
       !SLL_ALLOCATE( qp%particle_array(3+sim%no_weights,sim%n_particles), ierr)
       !allocate( qp%species, stat=ierr)
       !SLL_ASSERT(ierr==0)
       !call qp%species%init( sim%charge, 1.0_f64 )
       !qp%n_weights = sim%no_weights
    end select

    ! Initialize control variate
    allocate(sim%control_variate%cv(1) )
    call sim%control_variate%cv(1)%init(control_variate_equi, &
         distribution_params=sim%init_distrib_params)

    if (sim%rank == 0 ) then
       open(newunit=file_id, file=trim(sim%file_prefix)//'_parameters_used.dat')
       close(file_id)
    end if

    if( sim%boundary )then
       ! Initialize the field solver
       if (sim%ct) then
          allocate( sll_t_maxwell_clamped_1d_trafo :: sim%maxwell_solver )
          select type ( q=>sim%maxwell_solver )
          type is ( sll_t_maxwell_clamped_1d_trafo )
             call q%init_from_file( sim%domain(1:2), sim%n_gcells, sim%degree_smoother, sim%boundary_fields, sim%map,  trim(filename) )
          end select
       else
          allocate( sll_t_maxwell_clamped_1d_fem_sm :: sim%maxwell_solver )
          select type ( q=>sim%maxwell_solver )
          type is ( sll_t_maxwell_clamped_1d_fem_sm )
             call q%init_from_file( sim%domain(1:2), sim%n_gcells, &
                  sim%degree_smoother, sim%boundary_fields, trim(filename) )
          end select
       end if
       ! Initialize kernel smoother
       call sll_s_new_particle_mesh_coupling_spline_cl_1d(sim%kernel_smoother_1, &
            sim%domain(1:2), sim%n_gcells, &
            sim%n_particles, sim%degree_smoother-1, sll_p_galerkin, sim%boundary_fields) 
       call sll_s_new_particle_mesh_coupling_spline_cl_1d(sim%kernel_smoother_0, &
            sim%domain(1:2), sim%n_gcells, &
            sim%n_particles, sim%degree_smoother, sll_p_galerkin, sim%boundary_fields)
    else
       ! Initialize the field solver
       if (sim%ct) then
          allocate( sll_t_maxwell_1d_trafo :: sim%maxwell_solver )
          select type ( qtrafo=>sim%maxwell_solver )
          type is ( sll_t_maxwell_1d_trafo )
             call qtrafo%init_from_file( sim%domain(1:2), sim%n_gcells, sim%degree_smoother, sim%map,  trim(filename)  )
          end select
       else
          if ( sim%degree_fem > -1 ) then
             allocate( sll_t_maxwell_1d_fem :: sim%maxwell_solver )
             select type ( q=>sim%maxwell_solver )
             type is ( sll_t_maxwell_1d_fem )
                call q%init( sim%domain(1:2), sim%n_gcells, &
                     sim%degree_fem, delta_t*0.5_f64, strong_ampere = strong_ampere, &
                     force_sign=sim%force_sign, adiabatic_electrons=sim%adiabatic_electrons) ! Note: The time is defined for the matrices of the curl-curl solver. This is used for a half time step at the time (hence, the delta_t*0.5).
             end select
!!$             allocate( sll_t_maxwell_1d_fem_sm :: sim%maxwell_solver )
!!$             select type ( q=>sim%maxwell_solver )
!!$             type is ( sll_t_maxwell_1d_fem_sm )
!!$                call q%init_from_file( sim%domain(1:2), sim%n_gcells, &
!!$                     sim%degree_smoother, trim(filename) )
!!$             end select
          else
             allocate( sll_t_maxwell_1d_ps :: sim%maxwell_solver )
             select type ( q=>sim%maxwell_solver )
             type is ( sll_t_maxwell_1d_ps )
                call q%init( sim%domain(1:2), sim%n_gcells )
             end select
          end if
       end if
       ! Initialize kernel smoother
       if ( smoothing .eqv. .false. ) then
          if ( strong_ampere .eqv. .false. ) then
             call sll_s_new_particle_mesh_coupling_spline_1d(sim%kernel_smoother_1, &
                  sim%domain(1:2), sim%n_gcells, &
                  sim%n_particles, sim%degree_smoother-1, sll_p_galerkin) 
             call sll_s_new_particle_mesh_coupling_spline_1d(sim%kernel_smoother_0, &
                  sim%domain(1:2), sim%n_gcells, &
                  sim%n_particles, sim%degree_smoother, sll_p_galerkin)
          else
             if ( modulo(sim%degree_fem,2) == 0 ) then
                eval_grid_points = .false.
             else
                eval_grid_points = .true.
             end if

             call sll_s_new_particle_mesh_coupling_spline_strong_1d(sim%kernel_smoother_1, &
                  sim%domain(1:2), sim%n_gcells, &
                  sim%degree_smoother,  integ= .false., eval_grid_points=eval_grid_points ) 
             call sll_s_new_particle_mesh_coupling_spline_strong_1d(sim%kernel_smoother_0, &
                  sim%domain(1:2), sim%n_gcells, &
                  sim%degree_smoother+1, integ = .true., eval_grid_points=eval_grid_points )
          end if
       else
          call sll_s_new_particle_mesh_coupling_spline_smooth_1d(sim%kernel_smoother_1, &
               sim%domain(1:2), sim%n_gcells, &
               sim%n_particles, sim%degree_smoother-1, sll_p_galerkin) 
          call sll_s_new_particle_mesh_coupling_spline_smooth_1d(sim%kernel_smoother_0, &
               sim%domain(1:2), sim%n_gcells, &
               sim%n_particles, sim%degree_smoother, sll_p_galerkin) 
       end if
    end if

    ! Initialize the arrays for the spline coefficients of the fields
    SLL_ALLOCATE(sim%phi_dofs(sim%n_dofs0), ierr)
    SLL_ALLOCATE(sim%efield_dofs(sim%n_dofs0,2), ierr)
    SLL_ALLOCATE(sim%bfield_dofs(sim%n_dofs1), ierr)
    sim%phi_dofs = 0._f64
    sim%efield_dofs = 0._f64
    sim%bfield_dofs = 0._f64

    SLL_ALLOCATE(sim%rhob(1:sim%n_dofs0), ierr)
    SLL_ALLOCATE(sim%background_charge(1:sim%n_dofs0), ierr)
    sim%rhob = 0._f64
    sim%background_charge = 0.0_f64

    call add_background_charge(sim, sim%background_charge)

    ! Initialize the time-splitting propagator
    if (sim%splitting_case == sll_p_splitting_hs) then
       if (sim%no_weights == 1) then
          call sll_s_new_time_propagator_pic_vm_1d2v_hs(&
               sim%propagator, sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), sim%filter, boundary_particles=sim%boundary_particles, force_sign=sim%force_sign, jmean=jmean, betar = sim%plasma_betar(1:2), electrostatic=electrostatic, rhob = sim%rhob)
       else
          call sll_s_new_time_propagator_pic_vm_1d2v_hs(&
               sim%propagator, sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), sim%filter, sim%boundary_particles, sim%force_sign, jmean, sim%control_variate, sim%no_weights, betar = sim%plasma_betar(1:2), electrostatic=electrostatic, rhob = sim%rhob)
       end if
       sim%efield_dofs_n => sim%efield_dofs
    elseif( sim%splitting_case == sll_p_splitting_boris) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_boris :: sim%propagator )
       select type( qp=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_boris)
          call qp%init( sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3))
          sim%efield_dofs_n => qp%efield_dofs_mid
       end select
    elseif( sim%splitting_case == sll_p_splitting_disgradE) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_disgradE :: sim%propagator )
       select type( qpdisgradE=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_disgradE )
          if (sim%no_weights == 1) then
             call qpdisgradE%init_from_file( sim%maxwell_solver, &
                  sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(1), sim%domain(3), sim%filter, trim(filename), boundary_particles=sim%boundary_particles, force_sign=sim%force_sign, betar=sim%plasma_betar(1:2), electrostatic=electrostatic, jmean=jmean  )
             sim%efield_dofs_n => qpdisgradE%helper%efield_dofs
          else
             call qpdisgradE%init_from_file( sim%maxwell_solver, &
                  sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(1), sim%domain(3), sim%filter, trim(filename), sim%boundary_particles, sim%force_sign,  sim%control_variate, sim%no_weights, sim%plasma_betar(1:2), electrostatic, jmean=jmean  )
          end if
          sim%efield_dofs_n => qpdisgradE%helper%efield_dofs
       end select
    elseif  (sim%splitting_case == sll_p_splitting_cef) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_cef :: sim%propagator )
       select type( qpcef=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_cef )
          call qpcef%init( sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, sim%domain(1), sim%domain(3), &
               boundary_particles=sim%boundary_particles, electrostatic =electrostatic, rhob = sim%rhob )
          sim%efield_dofs_n => qpcef%efield_dofs
       end select
    elseif( sim%splitting_case == sll_p_splitting_ecsim) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_ecsim :: sim%propagator )
       select type( qpecsim=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_ecsim )
          call qpecsim%init_from_file(  sim%kernel_smoother_0, sim%kernel_smoother_1, &
               sim%particle_group, sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), trim(filename))
          sim%efield_dofs_n => qpecsim%efield_dofs
       end select
    elseif( sim%splitting_case == sll_p_splitting_ecsim2o) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_ecsim2o :: sim%propagator )
       select type( qpecsim2o=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_ecsim2o )
          call qpecsim2o%init(  sim%kernel_smoother_0, sim%kernel_smoother_1, &
               sim%particle_group, sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3))
          sim%efield_dofs_n => qpecsim2o%efield_dofs
       end select
    elseif  (sim%splitting_case == sll_p_splitting_disgradEC) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_disgradEC :: sim%propagator )
       select type( qpdisgradEC=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_disgradEC )
          call qpdisgradEC%init_from_file( sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), sim%filter, trim(filename), &
               boundary_particles=sim%boundary_particles, electrostatic =electrostatic)
          sim%efield_dofs_n => qpdisgradEC%helper%efield_dofs
       end select
    elseif  (sim%splitting_case == sll_p_splitting_disgradEC_sub) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub :: sim%propagator )
       select type( qpdgs=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub )
          call qpdgs%init_from_file( sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), sim%bfilter, trim(filename) )
          sim%efield_dofs_n => qpdgs%helper%efield_dofs
       end select
    elseif  (sim%splitting_case == sll_p_splitting_subcyc) then
       ! Read parameters from file
       open(newunit = input_file, file=trim(filename), status='old', IOStat=io_stat)
       if (io_stat /= 0) then
          print*, 'init_pic_1d2v() failed to open file ', filename
          STOP
       end if
       read(input_file, time_iterate)
       close(input_file)
       allocate( sll_t_time_propagator_pic_vm_1d2v_subcyc :: sim%propagator )
       select type( qpsc=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_subcyc )
          call qpsc%init(  sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), n_sub_iter, sim%file_prefix, &
               sim%bfilter, jmean, sim%control_variate, sim%no_weights)
          sim%efield_dofs_n => qpsc%efield_dofs
       end select
    elseif( sim%splitting_case == sll_p_splitting_hs_trafo) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_hs_trafo :: sim%propagator )
       select type( qphstrafo=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_hs_trafo )
          call qphstrafo%init( sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), sim%map, electrostatic=electrostatic )!,betar=sim%plasma_betar(1:2))
          sim%efield_dofs_n => qphstrafo%efield_dofs
       end select
    elseif( sim%splitting_case == sll_p_splitting_disgradE_trafo) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_disgradE_trafo :: sim%propagator )
       select type( qptrafo=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_disgradE_trafo )
          call qptrafo%init_from_file( sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), sim%map, trim(filename), sim%boundary_particles, force_sign=sim%force_sign, electrostatic=electrostatic, jmean=jmean  )!,betar=sim%plasma_betar(1:2))
          sim%efield_dofs_n => qptrafo%helper%efield_dofs
       end select
    elseif( sim%splitting_case == sll_p_splitting_disgradEC_trafo) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_disgradEC_trafo :: sim%propagator )
       select type( qptrafo=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_disgradEC_trafo )
          call qptrafo%init_from_file( sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), sim%map, trim(filename), sim%boundary_particles, force_sign=sim%force_sign, electrostatic=electrostatic, jmean=jmean  )!,betar=sim%plasma_betar(1:2))
          sim%efield_dofs_n => qptrafo%helper%efield_dofs
       end select
    elseif  (sim%splitting_case == sll_p_splitting_zigsub) then
       ! Read parameters from file
       open(newunit = input_file, file=trim(filename), status='old', IOStat=io_stat)
       if (io_stat /= 0) then
          print*, 'init_pic_1d2v() failed to open file ', filename
          STOP
       end if
       read(input_file, time_iterate)
       close(input_file)
       allocate( sll_t_time_propagator_pic_vm_1d2v_zigsub :: sim%propagator )
       select type( qpsc=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_zigsub )
          call qpsc%init(  sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), n_sub_iter, &
               sim%bfilter, jmean, sim%control_variate, sim%no_weights)
          sim%efield_dofs_n => qpsc%efield_dofs
       end select
    elseif (sim%splitting_case == sll_p_splitting_momentum) then
       call sll_s_new_time_propagator_pic_vm_1d2v_momentum(&
            sim%propagator, sim%maxwell_solver, &
            sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
            sim%efield_dofs, sim%bfield_dofs, &
            sim%domain(1), sim%domain(3), sim%bfilter, jmean, sim%control_variate, sim%no_weights)
       sim%efield_dofs_n => sim%efield_dofs
    end if

    call sll_s_set_time_mark( end_init )
    if (sim%rank == 0 ) then
       sim%time = sll_f_time_elapsed_between( start_init, end_init)
       write(*, "(A, F10.3)") "Init run time [s] = ", sim%time

       open(newunit=file_id, file=trim(sim%file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=ierr)
       write(file_id, *) 'delta t:', sim%delta_t
       write(file_id, *) 'n_time_steps:', sim%n_time_steps
       write(file_id, *) 'beta:', sim%beta
       write(file_id, *) 'charge:', sim%charge
       write(file_id, *) 'mass:', mass
       write(file_id, *) 'plasma betar:', sim%plasma_betar
       write(file_id, *) 'force sign simulation:', sim%force_sign
       write(file_id, *) 'electrostatic simulation:', electrostatic
       write(file_id, *) 'output filename:', sim%file_prefix
       write(file_id, *) 'output fields:', sim%output_fields
       write(file_id, *) 'output particles:', sim%output_particles
       write(file_id, *) 'n_cells:', sim%n_gcells
       write(file_id, *) 'domain:', sim%domain
       write(file_id, *) 'delta x:', sim%delta_x
       write(file_id, *) 'jmean:', jmean
       write(file_id, *) 'n_particles:', sim%n_total_particles
       write(file_id, *) 'spline degree:', sim%degree_smoother
       write(file_id, *) 'no_weights:', sim%no_weights
       close(file_id)
    end if

  end subroutine init_pic_vm_1d2v


  !------------------------------------------------------------------------------!
  !> Run simulation
  subroutine run_pic_vm_1d2v (sim)
    class(sll_t_sim_pic_vm_1d2v_cart), intent(inout) :: sim !< Singlespecies simulation
    ! Local variables
    sll_int32 :: j, ierr, i_part
    sll_real64, allocatable :: rho(:), rho_local(:), scratch(:)
    sll_int32 :: th_diag_id, file_id
    character(len=4) :: crank
    character(len=6) :: step
    character(len=256) :: diag_file_name
    sll_real64 :: wi(1)
    sll_real64 :: xi(3)
    type(sll_t_time_mark) :: start_loop, end_loop
    sll_int32 :: degreeb


    ! Initialize file for diagnostics
    if (sim%rank == 0) then
       diag_file_name = trim(sim%file_prefix)//"_thdiag.dat"

       call sll_s_ascii_file_create(trim(diag_file_name), th_diag_id, ierr)
    end if

    if(sim%ct) then
       call sim%sampler%sample( sim%particle_group%group(1), sim%init_distrib_params, sim%domain(1:1), sim%domain(3:3), sim%map )
    else
       if (sim%no_weights == 1) then
          call sim%sampler%sample( sim%particle_group%group(1), sim%init_distrib_params, sim%domain(1:1), sim%domain(3:3) )
       else
          call sim%sampler%sample_cv( sim%particle_group%group(1), sim%init_distrib_params, sim%domain(1:1), sim%domain(3:3), sim%control_variate%cv(1) )
       end if
    end if


    ! Print particle array to file
    if ( sim%output_particles ) then
       call sll_s_int2string( sim%rank, crank )
       call sim%particle_group%group(1)%print(trim(sim%file_prefix)//'_particles_start_'//crank//'.dat')
    end if

    ! Set the initial fields
    SLL_ALLOCATE(rho_local(sim%n_dofs0), ierr)
    SLL_ALLOCATE(rho(sim%n_dofs0), ierr)
    SLL_ALLOCATE(scratch(sim%n_dofs0), ierr)
    scratch = 0._f64
    ! Efield 1 by Poisson
    call solve_poisson( sim, rho_local, rho )
    ! Efield 2 to zero
    sim%efield_dofs(:,2) = 0.0_f64
    ! Bfield = beta*cos(kx): Use b = M{-1}(N_i,beta*cos(kx))

    if ( sim%strong_ampere .eqv. .true. ) then
       degreeb = sim%degree_fem
    else
       degreeb = sim%degree_fem-1
    end if

    select case( sim%initial_bfield )
    case (  sll_p_bfield_cos )
       call sim%maxwell_solver%L2projection( beta_cos_k, degreeb, &
            sim%bfield_dofs)
    case ( sll_p_bfield_sin )
       call sim%maxwell_solver%L2projection( beta_sin_k, degreeb, &
            sim%bfield_dofs)
    case ( sll_p_bfield_constant )
       call sim%maxwell_solver%L2projection( beta_constant, degreeb, &
            sim%bfield_dofs)
    case ( sll_p_bfield_cos_const )
       call sim%maxwell_solver%L2projection( beta_cos_const, sim%degree_smoother-1, &
            sim%bfield_dofs)
    end select
    ! End field initialization
    call sim%propagator%reinit_fields()


    ! In case we use the Boris propagator, we need to initialize the staggering used in the scheme.
    select type( qp=>sim%propagator )
    type is ( sll_t_time_propagator_pic_vm_1d2v_boris)
       call qp%staggering( sim%delta_t )
    end select



    ! Diagnostis
    call sll_s_time_history_diagnostics_pic_vm_1d2v( &
         sim, 0.0_f64, th_diag_id, rho, scratch)
    if (sim%rank == 0 ) then
       call sll_s_set_time_mark(start_loop )
    end if

    ! Time loop
    select case (sim%splitting_type)
    case(sll_p_strang_splitting)
       do j=1, sim%n_time_steps
          call sim%propagator%strang_splitting(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_1d2v( &
               sim,  sim%delta_t*real(j,f64), th_diag_id, rho, scratch)

          if( sim%output_particles ) then
             if( modulo(j, 100) == 0 ) then
                call sll_s_int2string( sim%rank, crank )
                call sll_s_int2string( j, step )
                call sim%particle_group%group(1)%print(trim(sim%file_prefix)//step//'_particles_'//crank//'.dat')
             end if
          end if

          if (sim%rank == 0 ) then
             if ( sim%output_fields ) then
                if( modulo(j, 100) == 0 ) then
                   call sll_s_int2string( j, step )
                   open(newunit=file_id, file=trim(sim%file_prefix)//step//'_efield.dat')
                   write(file_id, *) sim%efield_dofs
                   close(file_id)
                   open(newunit=file_id, file=trim(sim%file_prefix)//step//'_bfield.dat')
                   write(file_id, *) sim%bfield_dofs
                   close(file_id)
                end if
             end if
          end if

       end do
    case(sll_p_splitting_fourth)
       do j=1, sim%n_time_steps
          call sim%propagator%splitting_fourth(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_1d2v( &
               sim,  sim%delta_t*real(j,f64), th_diag_id, rho, scratch)
       end do
    case(sll_p_lie_splitting)
       do j=1, sim%n_time_steps
          call sim%propagator%lie_splitting(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_1d2v( &
               sim,  sim%delta_t*real(j,f64), th_diag_id, rho, scratch)
       end do
    case(sll_p_lie_splitting_back)
       do j=1, sim%n_time_steps
          call sim%propagator%lie_splitting_back(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_1d2v( &
               sim,  sim%delta_t*real(j,f64), th_diag_id, rho, scratch)
       end do
    case(sll_p_splitting_fourth_10steps)
       do j=1, sim%n_time_steps
          call sim%propagator%splitting_fourth_10steps(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_1d2v( &
               sim,  sim%delta_t*real(j,f64), th_diag_id, rho, scratch)
       end do
    case(sll_p_splitting_second_4steps)
       do j=1, sim%n_time_steps
          !print*, 'TIME STEP', j
          call sim%propagator%splitting_second_4steps(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_1d2v( &
               sim,  sim%delta_t*real(j,f64), th_diag_id, rho, scratch)
       end do
    case default
       print*, 'this splitting type is not implemented'
    end select



    if (sim%rank == 0 ) then
       call sll_s_set_time_mark( end_loop )
       write(*, "(A, F10.3)") "Main loop run time [s] = ", sll_f_time_elapsed_between( start_loop, end_loop)
       close(th_diag_id)
    end if

    ! Print particle array to file
    if (sim%output_particles ) then
       call sll_s_int2string( sim%rank, crank )
       call sll_s_int2string( sim%n_time_steps, step )
       call sim%particle_group%group(1)%print(trim(sim%file_prefix)//step//'_particles_'//crank//'.dat')
       if(sim%output_fields)then
          call sll_s_int2string( sim%n_time_steps, step )
          open(newunit=file_id, file=trim(sim%file_prefix)//step//'_efield.dat')
          write(file_id, *) sim%efield_dofs
          close(file_id)
          open(newunit=file_id, file=trim(sim%file_prefix)//step//'_bfield.dat')
          write(file_id, *) sim%bfield_dofs
          close(file_id)
       end if
    end if

    ! Part for ctest
    if ( sim%make_ctest .eqv. .true. ) then
       ! Compute final rho
       rho_local = 0.0_f64
       do i_part = 1, sim%particle_group%group(1)%n_particles
          xi = sim%particle_group%group(1)%get_x(i_part)
          wi(1) = sim%particle_group%group(1)%get_charge( i_part)
          call sim%kernel_smoother_0%add_charge(xi(1), wi(1), rho_local)
       end do
       ! MPI to sum up contributions from each processor
       rho = 0.0_f64
       call sll_o_collective_allreduce( sll_v_world_collective, &
            rho_local, &
            sim%n_dofs0, MPI_SUM, rho)
       if (sim%rank == 0) then
          !print*, rho
          call ctest( rho, rho_local, trim(sim%ctest_ref_file_rho), sim%ctest_passed )
          call sll_s_check_diagnostics(trim(sim%ctest_ref_file_thdiag),'ctest_thdiag.dat', 3E-13_f64, sim%ctest_passed)
       end if
    end if


  contains
    function beta_cos_k(x)
      sll_real64             :: beta_cos_k
      sll_real64, intent(in) :: x

      beta_cos_k = sim%beta * cos(sll_p_twopi*x/sim%domain(3))
    end function beta_cos_k

    function beta_sin_k(x)
      sll_real64             :: beta_sin_k
      sll_real64, intent(in) :: x

      beta_sin_k = sim%beta * sin(sll_p_twopi*x/sim%domain(3))
    end function beta_sin_k

    function beta_constant(x)
      sll_real64             :: beta_constant
      sll_real64, intent(in) :: x

      beta_constant = sim%beta
    end function beta_constant

    function beta_cos_const(x)
      sll_real64             :: beta_cos_const
      sll_real64, intent(in) :: x

      beta_cos_const = sim%beta + 0.001_f64 * cos(sll_p_twopi*x/sim%domain(3))

    end function beta_cos_const


  end subroutine run_pic_vm_1d2v


  !------------------------------------------------------------------------------!
  !> Finalize simulation
  subroutine delete_pic_vm_1d2v (sim)
    class(sll_t_sim_pic_vm_1d2v_cart), intent(inout) :: sim !< Singlespecies simulation

    SLL_ASSERT(storage_size(sim)>0)

    call sim%propagator%free()
    deallocate(sim%propagator)
    call sim%particle_group%group(1)%free()
    deallocate (sim%particle_group)
    call sim%maxwell_solver%free()
    deallocate(sim%maxwell_solver)
    deallocate(sim%kernel_smoother_0)
    deallocate(sim%kernel_smoother_1)
    call sim%control_variate%cv(1)%free()

    deallocate(sim%efield_dofs)
    deallocate(sim%bfield_dofs)

    call sim%init_distrib_params%free()
    deallocate(sim%init_distrib_params)
    call sim%sampler%free()

    if( sim%boundary ) then
       call sll_s_spline_pp_free_1d(sim%spline0_pp )
    end if

  end subroutine delete_pic_vm_1d2v


  !------------------------------------------------------------------------------!
  ! local subroutine to handle ctest
  subroutine ctest(rho_simulated, rho_ref, ctest_ref_file, passed)
    sll_real64, intent(in   ) :: rho_simulated(:) !< simulated charge
    sll_real64, intent(inout) :: rho_ref(:) !< reference charge
    character(*), intent(in   ) :: ctest_ref_file !< ctest reference file
    logical,    intent(  out) :: passed !< true if ctest checks out      

    ! For testing
    character(len=256) :: reffile
    sll_real64 :: error

    call sll_s_concatenate_filename_and_path( trim(ctest_ref_file), __FILE__,&
         reffile)
    call sll_s_read_data_real_array( reffile, rho_ref)

    rho_ref = rho_ref -  rho_simulated
    error = maxval(rho_ref)
    print*, 'Maximum error in rho is', error, '.'
    if (abs(error)> 1d-14) then
       passed = .FALSE.
    else
       passed = .TRUE.
    end if

  end subroutine ctest


  !> As a control variate, we use the equilibrium (v part of the initial distribution)
  function control_variate_equi( this, xi, vi, time) result(sll_f_control_variate)
    class(sll_t_control_variate) :: this !< control variate
    sll_real64, optional,  intent( in ) :: xi(:) !< particle position
    sll_real64, optional,  intent( in ) :: vi(:) !< particle velocity
    sll_real64, optional,  intent( in ) :: time  !< current time
    sll_real64               :: sll_f_control_variate


    sll_f_control_variate = &
         this%control_variate_distribution_params%eval_v_density( vi(1:2) )


  end function control_variate_equi


  !> Add background charge
  subroutine add_background_charge(sim, background_charge)
    class(sll_t_sim_pic_vm_1d2v_cart), intent(inout) :: sim !< Singlespecies simulation
    sll_real64,                        intent(  out) :: background_charge(:) !< Coefficient vector of the charge distribution
    !local variables
    sll_int32 :: i, quad, cell, index1d, q
    sll_real64, allocatable :: quad_xw(:,:), spline_val(:,:)
    sll_real64 :: jacobian, c

    background_charge=0._f64
    q = sim%degree_smoother+1
    allocate(quad_xw(2, q))
    allocate(spline_val( sim%degree_smoother+1, q ))
    quad_xw = sll_f_gauss_legendre_points_and_weights( q , 0.0_f64, 1.0_f64 )

    do i = 1, q
       call sll_s_uniform_bsplines_eval_basis( sim%degree_smoother, quad_xw(1,i), spline_val(:,i) )
    end do

    if( sim%ct) then
       if( sim%boundary )then
          do i = 1, sim%degree_smoother-1
             ! loop over Gauss points
             do quad = 1, q
                c =  (quad_xw(1,quad)+ real(i-1,f64))/real(sim%n_gcells, f64)
                jacobian = sim%map%jacobian( [c, 0._f64, 0._f64] )
                do cell = 1, sim%degree_smoother + 1
                   index1d = i+cell-1
                   background_charge(index1d) = background_charge(index1d) + &
                        sll_f_spline_pp_horner_1d( sim%degree_smoother, sim%spline0_pp%poly_coeffs_boundary_left(:,:,i), quad_xw(1,quad), cell) * quad_xw(2, quad) * abs(jacobian)
                end do
             enddo
          enddo
          do i = sim%degree_smoother, sim%n_gcells-sim%degree_smoother+1
             ! loop over Gauss points
             do quad = 1, q
                c =  (quad_xw(1,quad)+ real(i-1,f64))/real(sim%n_gcells, f64)
                jacobian = sim%map%jacobian( [c, 0._f64, 0._f64] )
                do cell = 1, sim%degree_smoother + 1
                   index1d = i+cell-1
                   background_charge(index1d) = background_charge(index1d) + &
                        spline_val(cell, quad) * quad_xw(2, quad) * abs(jacobian)
                end do
             enddo
          enddo
          do i = sim%n_gcells-sim%degree_smoother+2, sim%n_gcells 
             ! loop over Gauss points
             do quad = 1, q
                c =  (quad_xw(1,quad)+ real(i-1,f64))/real(sim%n_gcells, f64)
                jacobian = sim%map%jacobian( [c, 0._f64, 0._f64] )
                do cell = 1, sim%degree_smoother + 1
                   index1d = i+cell-1
                   background_charge(index1d) = background_charge(index1d) + &
                        sll_f_spline_pp_horner_1d( sim%degree_smoother, sim%spline0_pp%poly_coeffs_boundary_right(:,:,i-sim%n_gcells+sim%degree_smoother-1), quad_xw(1,quad), cell)* quad_xw(2, quad) * abs(jacobian)
                end do
             enddo
          enddo
          background_charge = background_charge/real(sim%n_gcells, f64)
       else
          do i = 1, sim%n_gcells
             do quad = 1, q
                c =  (quad_xw(1,quad)+ real(i-1,f64))/real(sim%n_gcells, f64)
                jacobian = sim%map%jacobian( [c, 0._f64, 0._f64] )

                do cell = 1, sim%degree_smoother + 1
                   index1d = modulo(i-sim%degree_smoother+cell-2,sim%n_gcells) + 1
                   background_charge(index1d) = background_charge(index1d) + &
                        spline_val(cell, quad) * quad_xw(2, quad) * abs(jacobian)
                end do
             end do
          end do
          background_charge = background_charge/real(sim%n_gcells, f64)
       end if
    else
       if( sim%boundary )then
          do i = 1, sim%degree_smoother-1
             ! loop over Gauss points
             do quad = 1, q
                do cell = 1, sim%degree_smoother + 1
                   index1d = i+cell-1
                   background_charge(index1d) = background_charge(index1d) + &
                        sll_f_spline_pp_horner_1d( sim%degree_smoother, sim%spline0_pp%poly_coeffs_boundary_left(:,:,i), quad_xw(1,quad), cell) * quad_xw(2, quad) 
                end do
             enddo
          enddo
          do i = sim%degree_smoother, sim%n_gcells-sim%degree_smoother+1
             ! loop over Gauss points
             do quad = 1, q
                do cell = 1, sim%degree_smoother + 1
                   index1d = i+cell-1
                   background_charge(index1d) = background_charge(index1d) + &
                        spline_val(cell, quad) * quad_xw(2, quad) 
                end do
             enddo
          enddo
          do i = sim%n_gcells-sim%degree_smoother+2, sim%n_gcells 
             ! loop over Gauss points
             do quad = 1, q
                do cell = 1, sim%degree_smoother + 1
                   index1d = i+cell-1
                   background_charge(index1d) = background_charge(index1d) + &
                        sll_f_spline_pp_horner_1d( sim%degree_smoother, sim%spline0_pp%poly_coeffs_boundary_right(:,:,i-sim%n_gcells+sim%degree_smoother-1), quad_xw(1,quad), cell)* quad_xw(2, quad) 
                end do
             enddo
          enddo
          background_charge = sim%delta_x * background_charge
       else
          background_charge = sim%delta_x
       end if
    end if

  end subroutine add_background_charge


  !------------------------------------------------------------------------------!
  !Diagnostic functions and other helper functions
  !> Diagnostics for PIC Vlasov-Maxwell 1d2v 
  !> @todo (should be part of the library)
  subroutine sll_s_time_history_diagnostics_pic_vm_1d2v(&
       sim,&
       time, &
       file_id, &
       rho, scratch)
    class(sll_t_sim_pic_vm_1d2v_cart), intent(inout) :: sim !< Singlespecies simulation
    sll_real64,                        intent(in   ) :: time !< time
    sll_int32,                         intent(in   ) :: file_id !< file ide
    sll_real64,                        intent(  out) :: rho(:) !< scratch data
    sll_real64,                        intent(  out) :: scratch(:) !< scratch data
    ! local variables
    sll_real64 :: diagnostics_local(4)
    sll_real64 :: diagnostics(4)
    sll_real64 :: potential_energy(3)
    sll_int32  :: i_part
    sll_int32  :: degree
    sll_real64 :: vi(3),  xi(3)
    sll_real64 :: wi(1)
    sll_real64 :: transfer(1), vvb(1), poynting
    sll_real64 :: error, phi, ecb(2)

    degree = sim%degree_fem

    diagnostics_local = 0.0_f64
    do i_part=1,sim%particle_group%group(1)%n_particles
       vi = sim%particle_group%group(1)%get_v(i_part)
       xi = sim%particle_group%group(1)%get_x(i_part)
       wi = sim%particle_group%group(1)%get_mass(i_part)

       ! Kinetic energy
       diagnostics_local(1) = diagnostics_local(1) + &
            0.5_f64*(vi(1)**2)*wi(1)
       diagnostics_local(2) = diagnostics_local(2) + &
            0.5_f64*(vi(2)**2)*wi(1)
       ! Momentum 1
       diagnostics_local(3) = diagnostics_local(3) + &
            vi(1)*wi(1)
       ! Momentum 2
       diagnostics_local(4) = diagnostics_local(4) + &
            vi(2)*wi(1)

    end do
    diagnostics = 0.0_f64
    call sll_s_collective_reduce_real64(sll_v_world_collective, diagnostics_local, 4,&
         MPI_SUM, 0, diagnostics)

    ! Add ExB part
    if( sim%boundary )then
       if( sim%ct ) then
          call compute_e_cross_b_curvilinear( sim%kernel_smoother_0, sim%kernel_smoother_1, degree, sim%map, sim%efield_dofs(1:sim%n_dofs1,1), sim%efield_dofs(:,2), sim%bfield_dofs, ecb)
       else
          call compute_e_cross_b_gauss( sim%kernel_smoother_0, sim%kernel_smoother_1, degree, sim%efield_dofs(1:sim%n_dofs1,1), sim%efield_dofs(:,2), sim%bfield_dofs, ecb)
       end if
       diagnostics(3:4) = diagnostics(3:4) + ecb 
    else
       if ( sim%strong_ampere .eqv. .false. ) then
          diagnostics(3) = diagnostics(3) + sim%maxwell_solver%inner_product( sim%efield_dofs(:,2), sim%bfield_dofs, degree, degree-1 )

          diagnostics(4) = diagnostics(4) - sim%maxwell_solver%inner_product( sim%efield_dofs(1:sim%n_dofs1,1), sim%bfield_dofs, degree-1 )
       else
          if ( degree == -1) then             
             diagnostics(3) = diagnostics(3) + sim%maxwell_solver%inner_product( sim%efield_dofs(:,2), sim%bfield_dofs, degree-1, degree )
             diagnostics(4) = diagnostics(4) - sim%maxwell_solver%inner_product( sim%efield_dofs(1:sim%n_dofs1,1), sim%bfield_dofs, degree, degree )
          else
             scratch(1:sim%n_gcells-1) = 0.5_f64 * (sim%efield_dofs(1:sim%n_gcells-1,2)+sim%efield_dofs(2:sim%n_gcells,2) )
             scratch(sim%n_gcells) = 0.5_f64 * (sim%efield_dofs(1,2)+sim%efield_dofs(sim%n_gcells,2) )
             diagnostics(3) = diagnostics(3) + sim%maxwell_solver%inner_product( sim%bfield_dofs, scratch,  degree )
             diagnostics(4) = diagnostics(4) - sim%maxwell_solver%inner_product( sim%efield_dofs(1:sim%n_dofs1,1), sim%bfield_dofs, degree )
          end if
       end if
    end if

    call sll_s_pic_diagnostics_transfer( sim%particle_group%group(1), sim%kernel_smoother_0, sim%kernel_smoother_1, sim%efield_dofs(1:sim%n_dofs1,1), sim%efield_dofs(:,2), transfer )

    call sll_s_pic_diagnostics_vvb( sim%particle_group%group(1), sim%kernel_smoother_1, &
         sim%bfield_dofs, vvb )
    call sll_s_pic_diagnostics_poynting( sim%maxwell_solver, degree, sim%efield_dofs(:,2), sim%bfield_dofs, &
         scratch, poynting )

    ! Check error in Gauss law
    call check_gauss_law( sim, rho, scratch, error )


    if (sim%rank == 0) then
       if (sim%strong_ampere .eqv. .false. ) then
          potential_energy(1) = sim%maxwell_solver%inner_product( sim%efield_dofs(1:sim%n_dofs1,1),  sim%efield_dofs_n(1:sim%n_dofs1,1), degree-1 )/sim%plasma_betar(2)
          potential_energy(2) = sim%maxwell_solver%inner_product( sim%efield_dofs(:,2),  sim%efield_dofs_n(:,2), degree )/sim%plasma_betar(2)
       else
          potential_energy(1) = sim%maxwell_solver%inner_product( sim%efield_dofs(:,1),  sim%efield_dofs_n(:,1), degree )/sim%plasma_betar(2)
          potential_energy(2) = sim%maxwell_solver%inner_product( sim%efield_dofs(1:sim%n_dofs1,2),  sim%efield_dofs_n(1:sim%n_dofs1,2), degree-1 )/sim%plasma_betar(2)
       end if
       if ( sim%strong_ampere .eqv. .false. ) then
          potential_energy(3) = sim%maxwell_solver%L2norm_squared( sim%bfield_dofs, degree-1 )*sim%plasma_betar(3)
       else
          potential_energy(3) = sim%maxwell_solver%L2norm_squared( sim%bfield_dofs, degree )*sim%plasma_betar(3)
       end if
       if(sim%adiabatic_electrons) then
          phi = sim%maxwell_solver%L2norm_squared( sim%phi_dofs, degree )
          write(file_id,'(f12.5,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16)' ) &
               time,  0.5_f64*potential_energy, 0.5_f64*phi, diagnostics(1), diagnostics(2), &
               sum(diagnostics(1:2)) + 0.5_f64*phi, diagnostics(3:4), -transfer+vvb+poynting, &
               error

       else
          write(file_id,'(f12.5,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16)' ) &
               time, 0.5_f64*potential_energy, diagnostics(1:2), &
               sum(diagnostics(1:2)) + sim%force_sign * 0.5_f64*sum(potential_energy), diagnostics(3:4), -transfer+vvb+poynting, &
               error
       end if
    end if

  end subroutine sll_s_time_history_diagnostics_pic_vm_1d2v


  !> Compute ExB with Gauss quadrature
  subroutine compute_e_cross_b_gauss( kernel_smoother_0, kernel_smoother_1, degree, efield_dofs1, efield_dofs2, bfield_dofs, ecb)
    class(sll_c_particle_mesh_coupling_1d) :: kernel_smoother_0  !< Kernel smoother (order p+1)
    class(sll_c_particle_mesh_coupling_1d) :: kernel_smoother_1  !< Kernel smoother (order p)   
    sll_int32,  intent( in    )            :: degree     !< maximal spline deg
    sll_real64, intent( in    )            :: efield_dofs1(:) !< coefficients of efield
    sll_real64, intent( in    )            :: efield_dofs2(:) !< coefficients of efield
    sll_real64, intent( in    )            :: bfield_dofs(:)   !< DoFs describing the magnetic field
    sll_real64, intent(   out )            :: ecb(2) !< E cross B 
    !local variables
    sll_int32  :: row, quad, q
    sll_real64, allocatable :: xw_gauss(:,:)
    sll_real64 :: xi(1)
    sll_real64 :: efield(2), bfield

    efield = 0.0_f64
    bfield = 0.0_f64
    ecb = 0.0_f64

    q = 2*degree
    allocate(xw_gauss(2, q))
    xw_gauss = sll_f_gauss_legendre_points_and_weights( q , 0.0_f64, 1.0_f64 )
    do row = 1, kernel_smoother_0%n_cells
       do quad = 1, q
          xi = kernel_smoother_0%delta_x*(xw_gauss(1,quad) + real(row - 1,f64))
          call kernel_smoother_1%evaluate &
               (xi(1), efield_dofs1, efield(1))
          call kernel_smoother_0%evaluate &
               (xi(1), efield_dofs2, efield(2))
          call kernel_smoother_1%evaluate &
               (xi(1), bfield_dofs, bfield)

          ecb(1) = ecb(1) + bfield*efield(2) * xw_gauss(2,quad)*kernel_smoother_0%delta_x
          ecb(2) = ecb(2) - bfield*efield(1) * xw_gauss(2,quad)*kernel_smoother_0%delta_x
       end do
    end do

  end subroutine compute_e_cross_b_gauss


  !> Compute ExB with coordinate transformation
  subroutine compute_e_cross_b_curvilinear( kernel_smoother_0, kernel_smoother_1, degree, map, efield_dofs1, efield_dofs2, bfield_dofs, ecb)
    class(sll_c_particle_mesh_coupling_1d)     :: kernel_smoother_0  !< Kernel smoother (order p+1)
    class(sll_c_particle_mesh_coupling_1d)     :: kernel_smoother_1  !< Kernel smoother (order p)   
    sll_int32,  intent( in    )                :: degree     !< maximal spline deg
    type(sll_t_mapping_3d), intent( inout    ) :: map        !< coordinate transformation
    sll_real64, intent( in    )                :: efield_dofs1(:) !< coefficients of efield
    sll_real64, intent( in    )                :: efield_dofs2(:) !< coefficients of efield
    sll_real64, intent( in    )                :: bfield_dofs(:)   !< DoFs describing the magnetic field
    sll_real64, intent(   out )                :: ecb(2) !< E cross B 
    !local variables
    sll_int32  :: row, quad, q
    sll_real64, allocatable :: xw_gauss(:,:)
    sll_real64 :: xi(3)
    sll_real64 :: efield(2), bfield

    efield = 0.0_f64
    bfield = 0.0_f64
    ecb = 0.0_f64
    xi = 0._f64

    q = 2*degree
    allocate(xw_gauss(2, q))
    xw_gauss = sll_f_gauss_legendre_points_and_weights( q , 0.0_f64, 1.0_f64 )
    do row = 1, kernel_smoother_0%n_cells
       do quad = 1, q
          xi(1) = kernel_smoother_0%delta_x*(xw_gauss(1,quad) + real(row - 1,f64))
          call kernel_smoother_1%evaluate &
               (xi(1), efield_dofs1, efield(1))
          call kernel_smoother_0%evaluate &
               (xi(1), efield_dofs2, efield(2))
          call kernel_smoother_1%evaluate &
               (xi(1), bfield_dofs, bfield)

          ecb(1) = ecb(1) + bfield*efield(2) * xw_gauss(2,quad)*kernel_smoother_0%delta_x
          ecb(2) = ecb(2) - bfield*efield(1) * xw_gauss(2,quad)*kernel_smoother_0%delta_x/map%jacobian(xi )
       end do
    end do

  end subroutine compute_e_cross_b_curvilinear


  !> Accumulate rho and solve Poisson
  subroutine solve_poisson( sim, rho_local, rho )
    class(sll_t_sim_pic_vm_1d2v_cart), intent(inout) :: sim !< Singlespecies simulation
    sll_real64,                        intent(  out) :: rho_local(:) !< charge local to mpi processors
    sll_real64,                        intent(  out) :: rho(:) !< charge
    !local variables
    sll_int32 :: i_part
    sll_real64 :: xi(3), wi(1)
    sll_real64 :: rho_cv(sim%n_dofs0)

    rho_local = 0.0_f64
    do i_part = 1, sim%particle_group%group(1)%n_particles
       xi = sim%particle_group%group(1)%get_x(i_part)
       ! Get charge for accumulation of rho
       wi(1) = sim%particle_group%group(1)%get_charge( i_part, sim%no_weights)
       call sim%kernel_smoother_0%add_charge(xi(1), wi(1), rho_local)
    end do
    ! MPI to sum up contributions from each processor
    rho = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, rho_local, &
         sim%n_dofs0, MPI_SUM, rho)

    call sim%filter%apply_inplace( rho )

    if ( sim%no_weights == 3 ) then
       call sim%maxwell_solver%compute_rhs_from_function( one, sim%degree_smoother, rho_cv)
       rho = rho + sim%charge*rho_cv
    end if

    rho_local = 0._f64
    if ( sim%strong_ampere) then
       rho = rho - sim%charge
       ! Add this if there should not be an offset in rho
       rho = rho - sum(rho)/ real(sim%n_dofs0, f64)
    else
       ! Add neutralizing background distribution
       rho = sim%force_sign*(rho - sim%charge * sim%background_charge)
    end if

    ! Solve Poisson problem
    if( sim%adiabatic_electrons) then
       call sim%maxwell_solver%compute_phi_from_rho ( rho, sim%phi_dofs, sim%efield_dofs(1:sim%n_dofs1,1) )
    else
       rho = rho * sim%plasma_betar(2)
       call sim%maxwell_solver%compute_E_from_rho( rho, sim%efield_dofs(1:sim%n_dofs1,1) )
    end if

  contains

    function one( x)
      sll_real64 :: one
      sll_real64, intent(in) :: x

      one = 1.0_f64

    end function one

    function sin_l(eta)
      sll_real64             :: sin_l
      sll_real64, intent(in) :: eta
      !local variables
      sll_real64             :: x1
      sll_real64             :: k=0.8_f64
      sll_real64             :: alpha=0.01_f64

      if( sim%ct )then
         x1=sim%map%get_x1( [eta, 0._f64, 0._f64] )
      else
         x1 = sim%domain(3)*eta
      end if
      sin_l = alpha/k*sin(k*x1)
    end function sin_l


  end subroutine solve_poisson


  !> check Gauss' law
  subroutine check_gauss_law( sim, rho, rho_gauss, error )
    class(sll_t_sim_pic_vm_1d2v_cart), intent(inout) :: sim !< Singlespecies simulation
    sll_real64,                        intent(  out) :: rho(:) !< charge
    sll_real64,                        intent(  out) :: rho_gauss(:) !< charge local to mpi processes
    sll_real64,                        intent(  out) :: error !< error in Gauss' law
    !local variables
    sll_int32 :: i_part
    sll_real64 :: xi(3), wi(1)
    sll_real64 :: rho_cv(sim%n_dofs0)

    rho_gauss = 0.0_f64
    do i_part = 1, sim%particle_group%group(1)%n_particles
       xi = sim%particle_group%group(1)%get_x(i_part)
       if (xi(1) >= sim%domain(1) .and. xi(1) <= sim%domain(2)) then
          wi(1) = sim%particle_group%group(1)%get_charge( i_part, sim%no_weights)
          call sim%kernel_smoother_0%add_charge(xi(1), wi(1), rho_gauss)
       end if
    end do
    ! MPI to sum up contributions from each processor
    rho = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, &
         rho_gauss, sim%n_dofs0, MPI_SUM, rho)
    rho_gauss = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, &
         sim%rhob, sim%n_dofs0, MPI_SUM, rho_gauss)
    rho = rho + rho_gauss

    call sim%filter%apply_inplace( rho )

    if ( sim%no_weights == 3 ) then
       call sim%maxwell_solver%compute_rhs_from_function( one, sim%degree_smoother, rho_cv)
       rho = rho+sim%charge*rho_cv
    end if

    ! Add neutralizing background distribution
    if ( sim%strong_ampere) then
       rho = rho - sim%charge
    else
       rho = sim%force_sign*(rho - sim%charge * sim%background_charge)
    end if

    rho_gauss=0._f64
    if( sim%adiabatic_electrons) then
       call sim%maxwell_solver%multiply_mass( sim%phi_dofs, rho_gauss, sim%degree_smoother )
    else
       if( .not. sim%boundary ) then
          ! subtract the constant error in rho
          error = 0._f64
          do i_part = 1, sim%n_dofs0
             error = error + rho(i_part)
          end do
          error = error/real(sim%n_dofs0, f64)
          rho = rho -error
       end if
       call sim%maxwell_solver%compute_rho_from_e( sim%efield_dofs_n(:,1), rho_gauss )
    end if
    error = maxval(abs(rho-rho_gauss))

  contains

    function one( x)
      sll_real64 :: one
      sll_real64, intent(in) :: x

      one = 1.0_f64

    end function one

  end subroutine check_gauss_law


  !> Compute \sum(particles)w_p( v_1,p e_1(x_p) + v_2,p e_2(x_p))
  subroutine sll_s_pic_diagnostics_transfer ( particle_group, kernel_smoother_0, kernel_smoother_1, efield_dofs1, efield_dofs2, transfer)
    class(sll_c_particle_group_base), intent( in   )  :: particle_group !< particle group object
    class(sll_c_particle_mesh_coupling_1d) :: kernel_smoother_0  !< Kernel smoother (order p+1)
    class(sll_c_particle_mesh_coupling_1d) :: kernel_smoother_1  !< Kernel smoother (order p)   
    sll_real64, intent( in    ) :: efield_dofs1(:) !< coefficients of efield
    sll_real64, intent( in    ) :: efield_dofs2(:) !< coefficients of efield
    sll_real64, intent(   out ) :: transfer(1) !< result

    sll_int32 :: i_part
    sll_real64 :: xi(3), vi(3), wi, efield(2), transfer_local(1)

    transfer_local = 0.0_f64
    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x( i_part )
       wi = particle_group%get_charge( i_part )
       vi = particle_group%get_v( i_part )

       call kernel_smoother_1%evaluate &
            (xi(1), efield_dofs1, efield(1))
       call kernel_smoother_0%evaluate &
            (xi(1), efield_dofs2, efield(2))

       transfer_local(1) = transfer_local(1) + (vi(1) * efield(1) + vi(2) * efield(2))*wi

    end do

    call sll_o_collective_allreduce( sll_v_world_collective, transfer_local, 1, MPI_SUM, transfer )

  end subroutine sll_s_pic_diagnostics_transfer


  !> Compute \sum(particles) w_p v_1,p b(x_p) v_2,p
  subroutine sll_s_pic_diagnostics_vvb ( particle_group, kernel_smoother_1, bfield_dofs, vvb )
    class(sll_c_particle_group_base), intent( in   )  :: particle_group   !< particle group object
    class(sll_c_particle_mesh_coupling_1d), intent( inout ) :: kernel_smoother_1  !< Kernel smoother (order p)  
    sll_real64,           intent( in    ) :: bfield_dofs(:) !< coefficients of bfield
    sll_real64,           intent(   out ) :: vvb(1) !< result

    sll_int32 :: i_part
    sll_real64 :: xi(3), vi(3), wi, bfield, vvb_local(1)

    vvb_local = 0.0_f64
    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x( i_part )
       wi = particle_group%get_charge( i_part )
       vi = particle_group%get_v( i_part )

       call kernel_smoother_1%evaluate &
            ( xi(1), bfield_dofs, bfield )

       vvb_local = vvb_local + wi * vi(1) * vi(2) * bfield

    end do

    call sll_o_collective_allreduce( sll_v_world_collective, vvb_local, 1, MPI_SUM, vvb )

  end subroutine sll_s_pic_diagnostics_vvb


  !> Compute e^T M_0^{-1}  R^T b
  subroutine sll_s_pic_diagnostics_poynting ( maxwell_solver, degree, efield_dofs, bfield_dofs, scratch, poynting )
    class(sll_c_maxwell_1d_base) :: maxwell_solver !< maxwell solver object
    sll_int32, intent( in    ) :: degree !< degree of finite element
    sll_real64, intent( in    ) :: efield_dofs(:) !< coefficients of efield
    sll_real64, intent( in    ) :: bfield_dofs(:) !< coefficients of bfield
    sll_real64, intent(   out ) :: scratch(:) !< scratch data 
    sll_real64, intent(   out ) :: poynting !< value of  e^T M_0^{-1}  R^T b

    scratch = 0.0_f64
    ! Multiply B by M_0^{-1}  R^T
    call maxwell_solver%compute_e_from_b ( 1.0_f64, bfield_dofs, scratch )

    poynting =  maxwell_solver%inner_product( efield_dofs, scratch, degree, degree )

  end subroutine sll_s_pic_diagnostics_poynting


  !> Check diagnostics
  subroutine sll_s_check_diagnostics(reffile, simfile, tol_error, passed)
    character(*), intent(in) :: reffile !< Name of reference file (stored in same folder as source file)
    character(*), intent(in) :: simfile !< Name of file with simulation result
    sll_real64, intent(in)   :: tol_error !< tolerance 
    logical, intent(out)     :: passed !< true if diagnostics checks out

    sll_real64 :: error
    sll_real64 :: data_sim(3,11)
    sll_real64 :: data_ref(3,11)
    sll_int32  :: file_id
    sll_int32   :: io_stat
    character(len=256) :: reffile_full   
    ! Read simulation result
    open(newunit=file_id, file=simfile, status='old', action='read', IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'sll_s_check_diagnostics() failed to open file ', simfile
       STOP
    end if
    read(unit=file_id,fmt=*) data_sim
    close(file_id)

    ! Read reference
    call sll_s_concatenate_filename_and_path( trim(reffile), __FILE__,&
         reffile_full)
    open(newunit=file_id, file=reffile_full, status='old', action='read', IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'sll_s_check_diagnostics() failed to open file ', reffile
       STOP
    end if
    read(unit=file_id,fmt=*) data_ref
    close(file_id)


    ! Compare
    data_sim = data_sim - data_ref
    error = maxval(abs(data_sim))

    print*, 'Max error in time history diagnostics: ', error
    if ( passed .eqv. .true. ) then
       if (error < tol_error) then
          passed = .true.
       else
          passed = .false.
       end if
    end if

  end subroutine sll_s_check_diagnostics


end module sll_m_sim_pic_vm_1d2v_cart
