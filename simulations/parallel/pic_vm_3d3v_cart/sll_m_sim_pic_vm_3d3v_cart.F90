! Simulation of 3d3v Vlasov-Maxwell with simple PIC method, FEM with splines

! author: Katharina Kormann, Benedikt Perse, IPP

module sll_m_sim_pic_vm_3d3v_cart
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_ascii_io, only: &
       sll_s_ascii_file_create

  use sll_m_collective, only: &
       sll_o_collective_allreduce, &
       sll_s_collective_reduce_real64, &
       sll_f_get_collective_rank, &
       sll_f_get_collective_size, &
       sll_v_world_collective

  use sll_m_constants, only: &
       sll_p_pi

  use sll_m_control_variate, only: &
       sll_t_control_variate, &
       sll_t_control_variates

  use sll_m_filter_base_3d, only: &
       sll_c_filter_base_3d

  use sll_m_fft_filter_3d, only: &
       sll_t_fft_filter_3d

  use sll_m_gauss_legendre_integration, only: &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_initial_distribution, only : &
       sll_c_distribution_params, &
       sll_s_initial_distribution_new, &
       sll_t_params_cos_gaussian

  use sll_m_io_utilities, only : &
       sll_s_read_data_real_array, &
       sll_s_concatenate_filename_and_path

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_m_maxwell_3d_fem, only: &
       sll_t_maxwell_3d_fem

  use sll_m_maxwell_3d_fem_fft, only: &
       sll_t_maxwell_3d_fem_fft

  use sll_m_maxwell_3d_trafo_parallel, only:&
       sll_t_maxwell_3d_trafo_parallel

  use sll_m_maxwell_clamped_3d_fem, only: &
       sll_t_maxwell_clamped_3d_fem

  use sll_m_maxwell_clamped_3d_trafo_parallel, only:&
       sll_t_maxwell_clamped_3d_trafo_parallel

  use sll_mpi, only: &
       mpi_sum

  use sll_m_particle_mesh_coupling_base_3d, only: &
       sll_c_particle_mesh_coupling_3d

  use sll_m_particle_mesh_coupling_spline_3d_feec, only: &
       sll_t_particle_mesh_coupling_spline_3d_feec

  use sll_m_particle_mesh_coupling_spline_cl_3d_feec, only: &
       sll_t_particle_mesh_coupling_spline_cl_3d_feec

  use sll_m_particle_group_3d3v, only: &
       sll_t_particle_group_3d3v, &
       sll_s_new_particle_group_3d3v_ptr

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base, &
       sll_t_particle_array

  use sll_m_particle_sampling, only: &
       sll_t_particle_sampling

  use sll_m_profile_functions, only: &
       sll_t_profile_functions

  use sll_m_sim_base, only: &
       sll_c_simulation_base_class

  use sll_m_splines_pp, only: &
       sll_f_spline_pp_horner_1d, &
       sll_p_boundary_periodic, &
       sll_p_boundary_clamped, &
       sll_p_boundary_clamped_square, &
       sll_p_boundary_clamped_cubic

  use sll_m_time_propagator_base, only: &
       sll_c_time_propagator_base

  use sll_m_time_propagator_pic_vm_3d3v_cef, only: &
       sll_t_time_propagator_pic_vm_3d3v_cef

  use sll_m_time_propagator_pic_vm_3d3v_cef_trafo, only: &
       sll_t_time_propagator_pic_vm_3d3v_cef_trafo

  use sll_m_time_propagator_pic_vm_3d3v_cl_helper, only: &
       sll_p_boundary_particles_periodic, &
       sll_p_boundary_particles_singular, &
       sll_p_boundary_particles_reflection, &
       sll_p_boundary_particles_absorption

  use sll_m_time_propagator_pic_vm_3d3v_disgradE, only: &
       sll_t_time_propagator_pic_vm_3d3v_disgradE

  use sll_m_time_propagator_pic_vm_3d3v_disgradEC, only: &
       sll_t_time_propagator_pic_vm_3d3v_disgradEC

  use sll_m_time_propagator_pic_vm_3d3v_disgradEC_trafo, only: &
       sll_t_time_propagator_pic_vm_3d3v_disgradEC_trafo

  use sll_m_time_propagator_pic_vm_3d3v_disgradE_trafo, only: &
       sll_t_time_propagator_pic_vm_3d3v_disgradE_trafo

  use sll_m_time_propagator_pic_vm_3d3v_disgradE_trunc, only: &
       sll_t_time_propagator_pic_vm_3d3v_disgradE_trunc

  use sll_m_time_propagator_pic_vm_3d3v_hs, only: &
       sll_t_time_propagator_pic_vm_3d3v_hs

  use sll_m_time_propagator_pic_vm_3d3v_hs_trafo, only: &
       sll_t_time_propagator_pic_vm_3d3v_hs_trafo

  use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_t_time_mark

  use sll_m_utilities, only : &
       sll_s_int2string


  implicit none

  public :: &
       sll_t_sim_pic_vm_3d3v_cart

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: sll_p_splitting_hs = 0
  sll_int32, parameter :: sll_p_splitting_hs_trafo = 1
  sll_int32, parameter :: sll_p_splitting_boris = 2
  sll_int32, parameter :: sll_p_splitting_disgradE = 3
  sll_int32, parameter :: sll_p_splitting_disgradE_trafo = 4
  sll_int32, parameter :: sll_p_splitting_nldisgradE = 5
  sll_int32, parameter :: sll_p_splitting_cef = 6
  sll_int32, parameter :: sll_p_splitting_cef_trafo = 7
  sll_int32, parameter :: sll_p_splitting_disgradE_trunc = 8
  sll_int32, parameter :: sll_p_splitting_disgradEC = 10
  sll_int32, parameter :: sll_p_splitting_disgradEC_trafo = 11


  sll_int32, parameter :: sll_p_onegaussian = 0
  sll_int32, parameter :: sll_p_twogaussian = 1

  sll_int32, parameter :: sll_p_bfield_cossum = 0
  sll_int32, parameter :: sll_p_bfield_sumcos = 1
  sll_int32, parameter :: sll_p_bfield_prodcos = 2
  sll_int32, parameter :: sll_p_bfield_polynom = 3
  sll_int32, parameter :: sll_p_bfield_constant = 4

  sll_int32, parameter :: sll_p_strang_splitting = 0
  sll_int32, parameter :: sll_p_splitting_fourth = 1
  sll_int32, parameter :: sll_p_lie_splitting = 2
  sll_int32, parameter :: sll_p_lie_splitting_back = 3
  sll_int32, parameter :: sll_p_splitting_fourth_10steps = 4
  sll_int32, parameter :: sll_p_splitting_second_4steps = 5


  type, extends(sll_c_simulation_base_class) :: sll_t_sim_pic_vm_3d3v_cart

     ! Abstract particle group
     class(sll_t_particle_array), allocatable :: particle_group !< Particle group

     ! Control variate
     type(sll_t_control_variates), pointer :: control_variate !< control variate
     sll_int32 :: no_weights !< weight number

     sll_real64, allocatable :: phi_dofs(:) !< DoFs describing the scalar potential
     sll_real64, allocatable :: efield_dofs(:) !< DoFs describing the three components of the electric field
     sll_real64, allocatable :: bfield_dofs(:) !< DoFs describing the three components of the magnetic field
     sll_real64, allocatable :: rhob(:) !< charge at the boundary

     sll_real64 :: time !< time

     ! Abstract Maxwell solver
     class(sll_c_maxwell_3d_base), allocatable :: maxwell_solver !< Maxwell solver

     ! Abstract kernel smoothers
     class( sll_c_particle_mesh_coupling_3d), allocatable :: particle_mesh_coupling !< Particle  mesh coupling

     ! Cartesian mesh
     sll_real64 :: delta_x(3) !< Grid spacing
     sll_real64 :: volume !< cell volume

     ! Specific operator splitting
     class(sll_c_time_propagator_base), allocatable :: propagator !< time propagator object 
     sll_int32 :: splitting_case !< time splitting case
     sll_int32 :: splitting_type !< time splitting type

     ! Physical parameters
     class(sll_c_distribution_params), allocatable :: init_distrib_params !< distribution parameter
     sll_real64 :: bkx(3) !< modenumbers for initial magnetic field
     sll_real64 :: balpha(3) !< magnitude of initial magnetic field
     sll_real64 :: bphase_shift(3) !< phase shift in initial magnetic field(e.g. cos or sinus)
     sll_real64 :: bconstant(3) !< constant part of initial magnetic field
     sll_real64 :: domain(3,3) !< domain length
     type(sll_t_particle_sampling) :: sampler !< particle sampler
     sll_real64 :: plasma_betar(3) !< reciprocal of plasma beta
     sll_real64 :: force_sign  !< sign of particle force
     logical :: external_field = .false.  !< true for run with external electric field
     logical :: adiabatic_electrons = .false. !< true for run with adiabatic electrons
     logical :: fft = .false. !< true for fft filter

     type(sll_t_profile_functions) :: profile !< temperature profile
     class(sll_c_filter_base_3d), allocatable :: filter !< filter
     ! Simulation parameters
     sll_real64 :: delta_t !< time step
     sll_int32  :: n_time_steps !< amount of time steps
     sll_int32  :: n_particles !< amount of particles per mpi process
     sll_int32  :: n_total_particles !< total number of particles
     sll_int32  :: degree_smoother(3) !< spline degrees
     sll_int32  :: n_gcells(3) !< numer of gridcells for each direction
     sll_int32  :: n_totalcells !< product of gridcells
     sll_int32  :: n_totaldofs0 !< total number of Dofs for 0-form
     sll_int32  :: n_totaldofs1  !< total number of Dofs for 1-form
     logical    :: boundary = .false. !< true for non periodic field boundary
     sll_int32  :: boundary_fields(3) = 100 !< field boundary conditions
     sll_int32  :: boundary_particles = 100 !< particle boundary conditions

     ! Parameters for MPI
     sll_int32  :: rank !< mpi rank
     sll_int32  :: world_size !< mpi world size

     ! Case definitions
     sll_int32  :: initial_bfield !< case for intial magnetic field
     sll_real64 :: charge !< charge of particle species

     ! Output
     character(len=256)   :: file_prefix !< name of diagnostic file
     logical              :: output_fields !< logical for print out fields
     logical              :: output_particles !< logical for print out particles

     ! For ctest
     logical    :: ctest_passed = .false. !< logical for ctest
     logical    :: make_ctest = .false. !< logical for ctest
     character(len=256)   :: ctest_ref_file !< name of ctest reffile

     ! For restart
     logical    :: restart = .false. !< true if restart run
     character(len=256) :: restart_file !< name of restart data
     sll_int32 :: restart_steps = 0  !< restart steps

     !coordinate transformation
     type(sll_t_mapping_3d) :: map !< coordinate transformation
     logical                :: ct = .false. !< true if coordinate transformation is used

   contains
     procedure :: init_from_file => init_pic_vm_3d3v !< Initialization
     procedure :: run => run_pic_vm_3d3v !< Simulation
     procedure :: delete => delete_pic_vm_3d3v !< Finalization

  end type sll_t_sim_pic_vm_3d3v_cart


contains


  !------------------------------------------------------------------------------!
  ! Read in the simulation parameters from input file
  subroutine init_pic_vm_3d3v (sim, filename)
    class(sll_t_sim_pic_vm_3d3v_cart), intent(inout) :: sim !< Singlespecies simulation
    character(len=*),                  intent(in)    :: filename !< filename
    !local variables
    sll_int32   :: io_stat
    sll_int32   :: input_file, file_id
    sll_int32   :: ierr
    type(sll_t_time_mark) :: start_init, end_init
    sll_real64         :: delta_t
    sll_int32          :: n_time_steps
    character(len=256) :: initial_distrib
    character(len=256) :: initial_bfield
    sll_real64         :: charge = -1._f64
    sll_real64         :: mass = 1._f64
    sll_real64         :: plasma_beta(3) = 1._f64
    character(len=256) :: particle_force
    logical            :: electrostatic = .false.
    logical            :: restart = .false.
    character(len=256) :: restart_file
    sll_int32          :: restart_steps = 0
    sll_real64         :: bkx(3) = 0._f64
    sll_real64         :: balpha(3) = 0._f64
    sll_real64         :: bphase_shift(3) = 0._f64
    sll_real64         :: bconstant(3) = 0._f64
    character(len=256) :: file_prefix
    logical            :: output_fields = .false.
    logical            :: output_particles = .false.
    sll_int32          :: ng_x(3)
    sll_real64         :: x_min(3), x_max(3)
    logical            :: jmean = .false.
    sll_int32          :: n_particles
    character(len=256) :: sampling_case
    logical            :: delta_perturb = .false.
    sll_real64         :: delta_eps(6) 
    character(len=256) :: splitting_case
    sll_int32          :: spline_degree(3)
    character(len=256) :: splitting_type
    character(len=256) :: boundary_fields = "none"
    character(len=256) :: boundary_particles = "none"
    logical            :: with_control_variate = .false.
    logical            :: lindf = .false.
    character(len=256) :: filtering
    sll_int32          :: filter_iter = 0
    sll_int32          :: mode(3) = 0
    logical            :: fft_diag = .false.
    character(len=256) :: ctest_case ="No"
    sll_int32          :: n_sub_iter


    namelist /sim_params/         delta_t, n_time_steps, initial_distrib, initial_bfield, charge, mass, plasma_beta, particle_force, electrostatic, restart, restart_file, restart_steps
    namelist /bfield/             bkx, balpha, bphase_shift, bconstant

    namelist /output/             file_prefix, output_fields, output_particles

    namelist /grid_dims/          ng_x, x_min, x_max, jmean

    namelist /pic_params/         n_particles, sampling_case, delta_perturb, delta_eps, splitting_case, spline_degree, splitting_type, boundary_fields, boundary_particles, with_control_variate, lindf, filtering, filter_iter, mode, fft_diag
    namelist /ctest/              ctest_case
    namelist /time_iterate/       n_sub_iter


    call sll_s_set_time_mark( start_init )

    ! Read parameters from file
    open(newunit = input_file, file=trim(filename), IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'init_pic_3d3v() failed to open file ', filename
       STOP
    end if

    read(input_file, sim_params)
    if( restart ) then
       sim%restart = .true.
       sim%restart_file = restart_file
       sim%restart_steps = restart_steps
    end if
    call sim%profile%init(input_file)
    close(input_file)
    open(newunit = input_file, file=trim(filename), IOStat=io_stat)
    call sll_s_initial_distribution_new( trim(initial_distrib), [3,3], input_file, sim%init_distrib_params, sim%profile )
    read(input_file, bfield)
    read(input_file, output)
    read(input_file, grid_dims)
    read(input_file, pic_params)
    read(input_file, ctest)
    close(input_file)

    ! Set MPI paramseters
    sim%world_size = sll_f_get_collective_size(sll_v_world_collective)
    sim%rank = sll_f_get_collective_rank(sll_v_world_collective)

    ! Copy the read parameters into the simulation parameters
    sim%delta_t = delta_t
    sim%n_time_steps = n_time_steps
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
    case( "cossum")
       sim%bkx = bkx
       sim%balpha = balpha
       sim%bphase_shift = bphase_shift
       sim%bconstant = bconstant
       sim%initial_bfield = sll_p_bfield_cossum
    case( "sumcos")
       sim%bkx = bkx
       sim%balpha = balpha
       sim%bphase_shift = bphase_shift
       sim%bconstant = bconstant
       sim%initial_bfield = sll_p_bfield_sumcos
    case( "prodcos")
       sim%bkx = bkx
       sim%balpha = balpha
       sim%bphase_shift = bphase_shift
       sim%bconstant = bconstant
       sim%initial_bfield = sll_p_bfield_prodcos
    case( "polynom")
       sim%bkx = bkx
       sim%balpha = balpha
       sim%bconstant = bconstant
       sim%initial_bfield = sll_p_bfield_polynom
    case( "constant")
       sim%bconstant = bconstant
       sim%initial_bfield = sll_p_bfield_constant
    case default
       print*, '#initial bfield must be cossum, sumcos or constant as alternatives are not yet implemented'
    end select

    ! Output  
    sim%file_prefix = file_prefix
    sim%output_fields = output_fields
    sim%output_particles = output_particles

    sim%n_gcells = ng_x
    sim%n_totalcells = product(ng_x)
    sim%domain(:,1) = x_min
    sim%domain(:,2) = x_max
    sim%domain(:,3) = x_max - x_min
    sim%delta_x = (x_max - x_min)/real(ng_x, f64)
    sim%volume =  product(sim%delta_x)

    sim%n_particles = n_particles/sim%world_size
    sim%degree_smoother = spline_degree

    call sim%sampler%init( trim(sampling_case), [3,3], sim%n_particles, sim%rank, delta_perturb, delta_eps )
    sim%n_total_particles = sim%n_particles * sim%world_size

    !boundary conditions
    select case(trim(boundary_fields))
    case( "clamped" )
       sim%boundary = .true.
       sim%boundary_fields=sll_p_boundary_periodic
       if( sim%degree_smoother(1) == 2 )then
          sim%boundary_fields(1)=sll_p_boundary_clamped_square
       else if( sim%degree_smoother(1) == 3 )then
          sim%boundary_fields(1)=sll_p_boundary_clamped_cubic
       else
          sim%boundary_fields(1)=sll_p_boundary_clamped
       end if
       sim%n_totaldofs0 = (ng_x(1)+spline_degree(1))*ng_x(2)*ng_x(3)
       sim%n_totaldofs1 = (ng_x(1)+spline_degree(1)-1)*ng_x(2)*ng_x(3)
    case("periodic")
       sim%boundary = .true.
       sim%boundary_fields=sll_p_boundary_periodic
       sim%n_totaldofs0 = (ng_x(1)+spline_degree(1))*ng_x(2)*ng_x(3)
       sim%n_totaldofs1 = (ng_x(1)+spline_degree(1)-1)*ng_x(2)*ng_x(3)
    case default
       sim%boundary = .false.
       sim%n_totaldofs0 = sim%n_totalcells
       sim%n_totaldofs1 = sim%n_totalcells
    end select

    select case(trim(boundary_particles))
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
    case("splitting_nldisgradE")
       sim%splitting_case = sll_p_splitting_nldisgradE
    case("splitting_cef")
       sim%splitting_case = sll_p_splitting_cef
    case("splitting_disgradE_trunc")
       sim%splitting_case = sll_p_splitting_disgradE_trunc
    case("splitting_disgradEC")
       sim%splitting_case = sll_p_splitting_disgradEC
    case("splitting_disgradE_trafo")
       sim%splitting_case = sll_p_splitting_disgradE_trafo
       sim%ct=.true.
       !for curvilinear set domain in logical coordinates due to the init of mesh coupling
       sim%domain(:,1) = 0._f64
       sim%domain(:,2) = 1._f64
       sim%domain(:,3) = 1._f64
       sim%delta_x = 1._f64/real(ng_x, f64)
       !sim%volume =  product(sim%delta_x)
       call sim%map%init_from_file(filename)
    case("splitting_hs_trafo")
       sim%splitting_case = sll_p_splitting_hs_trafo
       sim%ct=.true.
       !for curvilinear set domain in logical coordinates due to the init of mesh coupling
       sim%domain(:,1) = 0._f64
       sim%domain(:,2) = 1._f64
       sim%domain(:,3) = 1._f64
       sim%delta_x = 1._f64/real(ng_x, f64)
       call sim%map%init_from_file(filename)
    case("splitting_cef_trafo")
       sim%splitting_case = sll_p_splitting_cef_trafo
       sim%ct=.true.
       !for curvilinear set domain in logical coordinates due to the init of mesh coupling
       sim%domain(:,1) = 0._f64
       sim%domain(:,2) = 1._f64
       sim%domain(:,3) = 1._f64
       sim%delta_x = 1._f64/real(ng_x, f64)
       call sim%map%init_from_file(filename)
    case("splitting_disgradEC_trafo")
       sim%splitting_case = sll_p_splitting_disgradEC_trafo
       sim%ct=.true.
       !for curvilinear set domain in logical coordinates due to the init of mesh coupling
       sim%domain(:,1) = 0._f64
       sim%domain(:,2) = 1._f64
       sim%domain(:,3) = 1._f64
       sim%delta_x = 1._f64/real(ng_x, f64)
       call sim%map%init_from_file(filename)
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

    ! Control variate
    if (with_control_variate .eqv. .true.) then
       sim%no_weights = 3
       allocate(sim%control_variate)
       allocate(sim%control_variate%cv(1))
       if( sim%adiabatic_electrons ) then
          call sim%control_variate%cv(1)%init( control_variate_xi, &
               distribution_params=sim%init_distrib_params )
       else
          call sim%control_variate%cv(1)%init( control_variate_equi, &
               distribution_params=sim%init_distrib_params )
       end if
    else
       sim%no_weights = 1
    end if

    ! filter
    select case( filtering )
    case( "fft" )
       allocate( sll_t_fft_filter_3d :: sim%filter )
    case default
       allocate( sll_t_fft_filter_3d :: sim%filter )
    end select
    call sim%filter%init( filter_iter, sim%n_gcells, mode )

    sim%fft = fft_diag

    !ctest
    select case( ctest_case)
    case ("hs")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_hs.dat"
    case("hs_trafo")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_hs_trafo.dat"
    case("cl_hs_trafo")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_cl_hs_trafo.dat"
    case ("cef")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_cef.dat"
    case ("cl_cef")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_cl_cef.dat"
    case("cef_trafo")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_cef_trafo.dat"
    case("cl_cef_trafo")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_cl_cef_trafo.dat"
    case("disgradE")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_disgradE.dat"
    case("cl_disgradE")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_cl_disgradE.dat"
    case("disgradE_trafo")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_disgradE_trafo.dat"
    case("cl_disgradE_trafo")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_cl_disgradE_trafo.dat"
    case("disgradE_trunc")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_disgradE_trunc.dat"
    case("disgradEC")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_disgradEC.dat"
    case("disgradEC_trafo")
       sim%make_ctest = .true.
       sim%ctest_ref_file = "reffile_pic_vm_3d3v_cart_disgradEC_trafo.dat"
    end select

    ! Initialize the particles    (mass set to 1.0 and charge set to -1.0)
    allocate( sim%particle_group )
    sim%particle_group%n_species = 1
    allocate( sll_t_particle_group_3d3v :: sim%particle_group%group(sim%particle_group%n_species) )
    select type ( qp => sim%particle_group%group(1) )
    type is (  sll_t_particle_group_3d3v )
       ! Note: This call produces a segmentation fault with the INTEL 17 compiler
       ! Therefore we manually initialize here
       ! TODO: Fix the problem with the init function

       !call qp%init(sim%n_particles, &
       !     sim%n_total_particles, -1.0_f64, 1.0_f64, sim%no_weights )
       qp%n_particles = sim%n_particles
       qp%n_total_particles = sim%n_total_particles

       SLL_ALLOCATE(qp%particle_array(6+sim%no_weights, qp%n_particles), ierr)
       qp%particle_array = 0._f64
       allocate(qp%species, stat=ierr)
       SLL_ASSERT( ierr == 0)
       call qp%species%init( sim%charge, mass )
       qp%n_weights = sim%no_weights
    end select
    !call sll_s_new_particle_group_3d3v_ptr(sim%particle_group, sim%n_particles, &
    !     sim%n_total_particles, -1.0_f64, 1.0_f64, sim%no_weights)


    if (sim%rank == 0 ) then
       open(newunit=file_id, file=trim(sim%file_prefix)//'_parameters_used.dat')
       close(file_id)
    end if
    if( sim%boundary ) then
       ! Initialize the field solver
       if (sim%ct) then
          if(sim%map%singular) then
             print*, 'Error: singular mapping not yet implemented'
          else
             allocate( sll_t_maxwell_clamped_3d_trafo_parallel :: sim%maxwell_solver )
             select type ( q=>sim%maxwell_solver )
             type is ( sll_t_maxwell_clamped_3d_trafo_parallel )
                call q%init_from_file( sim%domain(:,1:2), sim%n_gcells, sim%degree_smoother,sim%boundary_fields, sim%map, trim(filename), sim%adiabatic_electrons, sim%profile )
             end select
          end if
       else
          allocate( sll_t_maxwell_clamped_3d_fem :: sim%maxwell_solver )
          select type ( q=>sim%maxwell_solver )
          type is ( sll_t_maxwell_clamped_3d_fem )
             call q%init_from_file( sim%domain(:,1:2), sim%n_gcells, sim%degree_smoother,sim%boundary_fields, trim(filename), sim%adiabatic_electrons, sim%profile )
          end select
       end if
       ! Initialize kernel smoother
       allocate( sll_t_particle_mesh_coupling_spline_cl_3d_feec :: sim%particle_mesh_coupling)
       select type ( mc=>sim%particle_mesh_coupling )
       type is ( sll_t_particle_mesh_coupling_spline_cl_3d_feec )
          call mc%init( sim%n_gcells, sim%domain(:,1:2), &
               sim%degree_smoother, sim%boundary_fields, sim%n_particles)
       end select
    else
       ! Initialize the field solver
       if (sim%ct) then
          allocate( sll_t_maxwell_3d_trafo_parallel :: sim%maxwell_solver )
          select type ( q=>sim%maxwell_solver )
          type is ( sll_t_maxwell_3d_trafo_parallel )
             call q%init_from_file( sim%domain(:,1:2), sim%n_gcells, sim%degree_smoother,sim%map, trim(filename), sim%adiabatic_electrons, sim%profile )
          end select
       else
          if( sim%adiabatic_electrons ) then
             allocate( sll_t_maxwell_3d_fem :: sim%maxwell_solver )
             select type ( q=>sim%maxwell_solver )
             type is ( sll_t_maxwell_3d_fem )
                call q%init_from_file( sim%domain(:,1:2), sim%n_gcells, sim%degree_smoother, trim(filename), sim%adiabatic_electrons, sim%profile )
             end select
          else
             allocate( sll_t_maxwell_3d_fem_fft :: sim%maxwell_solver )
             select type ( q=>sim%maxwell_solver )
             type is ( sll_t_maxwell_3d_fem_fft )
                call q%init( sim%domain(:,1:2), sim%n_gcells, sim%degree_smoother )
             end select
          end if
       end if
       ! Initialize kernel smoother
       allocate( sll_t_particle_mesh_coupling_spline_3d_feec :: sim%particle_mesh_coupling)
       select type ( mc=>sim%particle_mesh_coupling )
       type is ( sll_t_particle_mesh_coupling_spline_3d_feec )
          call mc%init( sim%n_gcells, sim%domain(:,1:2), &
               sim%degree_smoother, sim%n_particles )
       end select
    end if

    ! Initialize the arrays for the spline coefficients of the fields
    SLL_ALLOCATE(sim%phi_dofs(1:sim%n_totaldofs0), ierr)
    SLL_ALLOCATE(sim%efield_dofs(1:sim%n_totaldofs1+sim%n_totaldofs0*2), ierr)
    SLL_ALLOCATE(sim%bfield_dofs(1:sim%n_totaldofs0+sim%n_totaldofs1*2), ierr)
    sim%phi_dofs = 0._f64
    sim%efield_dofs = 0._f64
    sim%bfield_dofs = 0._f64

    SLL_ALLOCATE(sim%rhob(1:sim%n_totaldofs0), ierr)
    sim%rhob = 0._f64

    ! Initialize the time-splitting propagator
    if (sim%splitting_case == sll_p_splitting_hs) then
       allocate( sll_t_time_propagator_pic_vm_3d3v_hs :: sim%propagator )
       select type( qp=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_3d3v_hs )
          if (sim%no_weights == 1) then
             call qp%init( sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%filter, &
                  sim%boundary_particles, betar=sim%plasma_betar(1:2), &
                  force_sign=sim%force_sign, electrostatic=electrostatic, &
                  rhob = sim%rhob, jmean = jmean)
          else
             call qp%init( sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%filter, &
                  sim%boundary_particles, betar=sim%plasma_betar(1:2), &
                  force_sign=sim%force_sign, electrostatic=electrostatic, &
                  rhob = sim%rhob,  control_variate=sim%control_variate, &
                  jmean = jmean, lindf = lindf)
          end if
       end select
    elseif (sim%splitting_case == sll_p_splitting_disgradE) then
       allocate( sll_t_time_propagator_pic_vm_3d3v_disgradE :: sim%propagator )
       select type( qpdisgradE=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_3d3v_disgradE )
          if (sim%no_weights == 1) then
             call qpdisgradE%init_from_file(sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%filter, &
                  trim(filename), sim%boundary_particles, betar=sim%plasma_betar(1:2), &
                  force_sign=sim%force_sign, electrostatic=electrostatic, rhob = sim%rhob, jmean = jmean )
          else
             call qpdisgradE%init_from_file(sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%filter, &
                  trim(filename), sim%boundary_particles, betar=sim%plasma_betar(1:2), &
                  force_sign=sim%force_sign, electrostatic=electrostatic, rhob = sim%rhob, &
                  control_variate=sim%control_variate, lindf = lindf, jmean = jmean )
          end if
       end select
    elseif (sim%splitting_case == sll_p_splitting_disgradEC) then
       if(sim%boundary ) then
          print*, 'Error: not implemented'
       else
          allocate( sll_t_time_propagator_pic_vm_3d3v_disgradEC :: sim%propagator )
          select type( qpdg=>sim%propagator )
          type is ( sll_t_time_propagator_pic_vm_3d3v_disgradEC )
             if (sim%no_weights == 1) then
                call qpdg%init_from_file(sim%maxwell_solver, &
                     sim%particle_mesh_coupling, sim%particle_group, &
                     sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                     sim%domain(:,1), sim%domain(:,3), sim%filter, &
                     trim(filename), betar=sim%plasma_betar(1:2), force_sign=sim%force_sign, electrostatic=electrostatic, jmean = jmean)
             else
                call qpdg%init_from_file(sim%maxwell_solver, &
                     sim%particle_mesh_coupling, sim%particle_group, &
                     sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                     sim%domain(:,1), sim%domain(:,3), sim%filter, &
                     trim(filename), betar=sim%plasma_betar(1:2), force_sign=sim%force_sign, electrostatic=electrostatic, jmean = jmean)
             end if
          end select
       end if
    elseif (sim%splitting_case == sll_p_splitting_disgradE_trafo) then
       allocate( sll_t_time_propagator_pic_vm_3d3v_disgradE_trafo :: sim%propagator )
       select type( qdisgradEtrafo=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_3d3v_disgradE_trafo )
          if (sim%no_weights == 1) then
             call qdisgradEtrafo%init_from_file(sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%map, &
                  trim(filename), sim%boundary_particles, betar=sim%plasma_betar(1:2), &
                  force_sign=sim%force_sign, electrostatic=electrostatic, &
                  rhob = sim%rhob, jmean = jmean)
          else
             call qdisgradEtrafo%init_from_file(sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%map, &
                  trim(filename), sim%boundary_particles, betar=sim%plasma_betar(1:2), &
                  force_sign=sim%force_sign, electrostatic=electrostatic, &
                  control_variate=sim%control_variate, rhob = sim%rhob, lindf= lindf, jmean = jmean)
          end if
       end select
    elseif (sim%splitting_case == sll_p_splitting_hs_trafo) then
       allocate( sll_t_time_propagator_pic_vm_3d3v_hs_trafo :: sim%propagator )
       select type( qhstrafo=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_3d3v_hs_trafo )
          if (sim%no_weights == 1) then
             call qhstrafo%init_from_file(sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%map, &
                  trim(filename), sim%boundary_particles, betar=sim%plasma_betar(1:2), electrostatic=electrostatic, rhob = sim%rhob)
          else
             call qhstrafo%init_from_file(sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%map, &
                  trim(filename), sim%boundary_particles, betar=sim%plasma_betar(1:2), electrostatic=electrostatic, rhob = sim%rhob, &
                  control_variate=sim%control_variate, lindf = lindf)
          end if
       end select
    elseif (sim%splitting_case == sll_p_splitting_cef_trafo) then
       allocate( sll_t_time_propagator_pic_vm_3d3v_cef_trafo :: sim%propagator )
       select type( qceftrafo=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_3d3v_cef_trafo )
          if (sim%no_weights == 1) then
             call qceftrafo%init_from_file(sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%map, &
                  trim(filename), boundary_particles=sim%boundary_particles, betar=sim%plasma_betar(1:2), electrostatic=electrostatic, rhob = sim%rhob)
          else
             call qceftrafo%init_from_file(sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%map, &
                  trim(filename), boundary_particles=sim%boundary_particles, betar=sim%plasma_betar(1:2), electrostatic=electrostatic, rhob = sim%rhob, control_variate=sim%control_variate, lindf = lindf)
          end if
       end select
    elseif (sim%splitting_case == sll_p_splitting_disgradEC_trafo) then
       if(sim%boundary ) then
          print*, 'Error: not implemented'
       else
          allocate( sll_t_time_propagator_pic_vm_3d3v_disgradEC_trafo :: sim%propagator )
          select type( qdisgtrafo=>sim%propagator )
          type is ( sll_t_time_propagator_pic_vm_3d3v_disgradEC_trafo )
             call qdisgtrafo%init_from_file(sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%map, &
                  trim(filename), betar=sim%plasma_betar(1:2), electrostatic=electrostatic )
          end select
       end if
    elseif (sim%splitting_case == sll_p_splitting_cef) then
       allocate( sll_t_time_propagator_pic_vm_3d3v_cef :: sim%propagator )
       select type( qpcef=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_3d3v_cef )
          if (sim%no_weights == 1) then
             call qpcef%init(sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%boundary_particles, &
                  betar=sim%plasma_betar(1:2), electrostatic=electrostatic, rhob = sim%rhob )
          else
             call qpcef%init(sim%maxwell_solver, &
                  sim%particle_mesh_coupling, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(:,1), sim%domain(:,3), sim%boundary_particles, &
                  betar=sim%plasma_betar(1:2), electrostatic=electrostatic, rhob = sim%rhob, control_variate=sim%control_variate, lindf = lindf )
          end if
       end select
    elseif (sim%splitting_case == sll_p_splitting_disgradE_trunc) then
       allocate( sll_t_time_propagator_pic_vm_3d3v_disgradE_trunc :: sim%propagator )
       select type( qpdisgradE=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_3d3v_disgradE_trunc )
          call qpdisgradE%init(sim%maxwell_solver, &
               sim%particle_mesh_coupling, sim%particle_group, &
               sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(:,1), sim%domain(:,3), sim%filter)
       end select
    else
       print*, 'Propagator not implemented.'
    end if

    call sll_s_set_time_mark( end_init )
    if (sim%rank == 0 ) then
       sim%time = sll_f_time_elapsed_between( start_init, end_init)
       write(*, "(A, F10.3)") "Init run time [s] = ", sim%time

       open(newunit=file_id, file=trim(sim%file_prefix)//'_parameters_used.dat', position = 'append', status='old', action='write', iostat=ierr)
       write(file_id, *) 'delta t:', sim%delta_t
       write(file_id, *) 'n_time_steps:', sim%n_time_steps
       write(file_id, *) 'charge:', sim%charge
       write(file_id, *) 'mass:', mass
       write(file_id, *) 'plasma betar:', sim%plasma_betar
       write(file_id, *) 'force sign simulation:', sim%force_sign
       write(file_id, *) 'electrostatic simulation:', electrostatic
       write(file_id, *) 'output filename:', sim%file_prefix
       write(file_id, *) 'output fields:', sim%output_fields
       write(file_id, *) 'output particles:', sim%output_particles
       write(file_id, *) 'n_cells:', sim%n_gcells
       write(file_id, *) 'domain:', sim%domain(:,1:2)
       write(file_id, *) 'delta x:', sim%delta_x
       write(file_id, *) 'jmean:', jmean
       write(file_id, *) 'n_particles:', sim%n_total_particles
       write(file_id, *) 'spline degree:', sim%degree_smoother
       write(file_id, *) 'no_weights:', sim%no_weights
       write(file_id, *) 'control_variate', with_control_variate
       write(file_id, *) 'linear delta f', lindf
       close(file_id)
    end if

  end subroutine init_pic_vm_3d3v


  !------------------------------------------------------------------------------!
  !> Run simulation
  subroutine run_pic_vm_3d3v (sim)
    class(sll_t_sim_pic_vm_3d3v_cart), intent(inout) :: sim !< Singlespecies simulation
    ! Local variables
    sll_int32  :: j, ierr
    sll_real64, allocatable :: rho(:), rho_local(:), scratch(:)
    sll_int32  :: th_diag_id, file_id
    character(len=4) :: crank
    character(len=6) :: step
    character(len=256) :: diag_file_name
    type(sll_t_time_mark) :: start_loop, end_loop

    ! Initialize file for diagnostics
    if (sim%rank == 0) then
       if(sim%restart) then
          open(newunit=th_diag_id, file=trim(sim%restart_file)//'_diag.dat', position = 'append', status='old', action='write', iostat=ierr)
       else
          diag_file_name = trim(sim%file_prefix)//"_diag.dat"
          call sll_s_ascii_file_create(trim(diag_file_name), th_diag_id, ierr)
       end if
    end if

    if ( sim%restart ) then
       call sll_s_int2string( sim%rank, crank )
       call sll_s_int2string( sim%restart_steps, step )
       call sim%particle_group%group(1)%read(trim(sim%restart_file)//step//'_particles_'//crank//'.dat')
    else
       if(sim%ct) then
          if (sim%no_weights == 1) then
             call sim%sampler%sample( sim%particle_group%group(1), sim%init_distrib_params, &
                  sim%map%get_x([0._f64,0._f64,0._f64]), sim%map%Lx, map=sim%map )
          else
             call sim%sampler%sample_cv( sim%particle_group%group(1), sim%init_distrib_params, &
                  sim%map%get_x([0._f64,0._f64,0._f64]), sim%map%Lx, sim%control_variate%cv(1), map=sim%map )
          end if
       else
          if (sim%no_weights == 1) then
             call sim%sampler%sample( sim%particle_group%group(1), sim%init_distrib_params, &
                  sim%domain(:,1), sim%domain(:,3) )
          else
             call sim%sampler%sample_cv( sim%particle_group%group(1), sim%init_distrib_params, &
                  sim%domain(:,1), sim%domain(:,3), sim%control_variate%cv(1) )
          end if
       end if
    end if

    ! Print particle array to file
    if ( sim%output_particles ) then
       call sll_s_int2string( sim%rank, crank )
       call sim%particle_group%group(1)%print(trim(sim%file_prefix)//'_particles_start_'//crank//'.dat')
    end if

    ! Set the initial fields
    SLL_ALLOCATE(rho_local(1:sim%n_totaldofs0), ierr)
    SLL_ALLOCATE(rho(1:sim%n_totaldofs0), ierr)
    SLL_ALLOCATE(scratch(1:sim%n_totaldofs0), ierr)
    scratch = 0.0_f64


    if ( sim%restart ) then
       call sll_s_int2string( sim%restart_steps, step )
       open(newunit=file_id, file=trim(sim%restart_file)//step//'_efield.dat', status='old', action='read', iostat=ierr)
       if (ierr /= 0 ) then
          SLL_ERROR("run", "Restart file for efield does not exist: "//trim(sim%restart_file)//step//'_efield.dat')
       end if
       read(file_id, *) sim%efield_dofs
       close(file_id)
       open(newunit=file_id, file=trim(sim%restart_file)//step//'_bfield.dat', status='old', action='read', iostat=ierr)
       if (ierr /= 0 ) then
          SLL_ERROR("run", "Restart file for bfield does not exist: "//trim(sim%restart_file)//step//'_bfield.dat')
       end if
       read(file_id, *) sim%bfield_dofs
       close(file_id)
    else
       ! Efield  by Poisson
       call solve_poisson( sim, rho_local, rho )

       if(sim%ct) then
          select case( sim%initial_bfield )
          case (  sll_p_bfield_cossum )
             call sim%maxwell_solver%L2projection( 2, 0, sim%bfield_dofs, bx_cossum_k, by_cossum_k, bz_cossum_k )
          case (  sll_p_bfield_sumcos )
             call sim%maxwell_solver%L2projection( 2, 0, sim%bfield_dofs, zero, zero, bfield_sumcos_k )
          case( sll_p_bfield_prodcos)
             call sim%maxwell_solver%L2projection( 2, 0, sim%bfield_dofs, bx_prodcos_k, by_prodcos_k, bz_prodcos_k )
          case (  sll_p_bfield_polynom )
             call sim%maxwell_solver%L2projection( 2, 0, sim%bfield_dofs, bx_polynom_k, by_polynom_k, bz_polynom_k )
          case (  sll_p_bfield_constant )
             call sim%maxwell_solver%L2projection( 2, 0, sim%bfield_dofs, zero, zero, bfield_constant_k )
          end select
       else
          select case( sim%initial_bfield )
          case (  sll_p_bfield_cossum )
             call sim%maxwell_solver%L2projection( 2, 1, sim%bfield_dofs(1:sim%n_totaldofs0), bx_cossum_k )
             call sim%maxwell_solver%L2projection( 2, 2, sim%bfield_dofs(sim%n_totaldofs0+1:sim%n_totaldofs0+sim%n_totaldofs1), by_cossum_k )
             call sim%maxwell_solver%L2projection( 2, 3, sim%bfield_dofs(sim%n_totaldofs0+sim%n_totaldofs1+1:sim%n_totaldofs0+sim%n_totaldofs1*2), bz_cossum_k )
          case (  sll_p_bfield_sumcos )
             call sim%maxwell_solver%L2projection( 2, 3, sim%bfield_dofs(sim%n_totaldofs0+sim%n_totaldofs1+1:sim%n_totaldofs0+sim%n_totaldofs1*2), bfield_sumcos_k )
          case( sll_p_bfield_prodcos)
             call sim%maxwell_solver%L2projection( 2, 1, sim%bfield_dofs(1:sim%n_totaldofs0), bx_prodcos_k )
             call sim%maxwell_solver%L2projection( 2, 2, sim%bfield_dofs(sim%n_totaldofs0+1:sim%n_totaldofs0+sim%n_totaldofs1), by_prodcos_k )
             call sim%maxwell_solver%L2projection( 2, 3, sim%bfield_dofs(sim%n_totaldofs0+sim%n_totaldofs1+1:sim%n_totaldofs0+sim%n_totaldofs1*2), bz_prodcos_k )
          case (  sll_p_bfield_constant )
             call sim%maxwell_solver%L2projection( 2, 3, sim%bfield_dofs(sim%n_totaldofs0+sim%n_totaldofs1+1:sim%n_totaldofs0+sim%n_totaldofs1*2), bfield_constant_k )
          end select
       end if
    end if

    ! Diagnostics
    if(sim%restart)then
    else
       call sll_s_time_history_diagnostics_pic_vm_3d3v( &
            sim,  0.0_f64, th_diag_id, rho_local, rho )
    end if
    if (sim%rank == 0 ) then
       call sll_s_set_time_mark( start_loop )
    end if


    select case (sim%splitting_type)
    case(sll_p_strang_splitting)
       do j=1+sim%restart_steps, sim%n_time_steps+sim%restart_steps
          call sim%propagator%strang_splitting(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_3d3v( &
               sim, sim%delta_t*real(j,f64), th_diag_id, rho, scratch)

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
       do j=1+sim%restart_steps, sim%n_time_steps+sim%restart_steps
          call sim%propagator%splitting_fourth(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_3d3v( &
               sim, sim%delta_t*real(j,f64), th_diag_id, rho, scratch)
       end do
    case(sll_p_lie_splitting)
       do j=1+sim%restart_steps, sim%n_time_steps+sim%restart_steps
          call sim%propagator%lie_splitting(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_3d3v( &
               sim, sim%delta_t*real(j,f64), th_diag_id, rho, scratch)
       end do
    case(sll_p_lie_splitting_back)
       do j=1+sim%restart_steps, sim%n_time_steps+sim%restart_steps
          call sim%propagator%lie_splitting_back(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_3d3v( &
               sim, sim%delta_t*real(j,f64), th_diag_id, rho, scratch)
       end do
    case(sll_p_splitting_fourth_10steps)
       do j=1+sim%restart_steps, sim%n_time_steps+sim%restart_steps
          call sim%propagator%splitting_fourth_10steps(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_3d3v( &
               sim, sim%delta_t*real(j,f64), th_diag_id, rho, scratch)
       end do
    case(sll_p_splitting_second_4steps)
       do j=1+sim%restart_steps, sim%n_time_steps+sim%restart_steps
          call sim%propagator%splitting_second_4steps(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_3d3v( &
               sim, sim%delta_t*real(j,f64), th_diag_id, rho, scratch)
       end do
    case default
       print*, 'this splitting type is not implemented'
    end select

    if (sim%rank == 0 ) then
       call sll_s_set_time_mark( end_loop )
       write(*, "(A, F10.3)") "Main loop run time [s] = ", sll_f_time_elapsed_between( start_loop, end_loop)
       close(th_diag_id)
       if ( sim%output_fields ) then
          call sll_s_int2string( sim%n_time_steps+sim%restart_steps, step )
          open(newunit=file_id, file=trim(sim%file_prefix)//step//'_efield.dat')
          write(file_id, *) sim%efield_dofs
          close(file_id)
          open(newunit=file_id, file=trim(sim%file_prefix)//step//'_bfield.dat')
          write(file_id, *) sim%bfield_dofs
          close(file_id)
       end if

       if ( sim%make_ctest .eqv. .true. ) then
          ! Check for ctest
          call sll_s_check_diagnostics(trim(sim%ctest_ref_file),"ctest_diag.dat" , 8E-13_f64, sim%ctest_passed)
       end if
    end if

    ! Print particle array to file
    if ( sim%output_particles ) then
       call sll_s_int2string( sim%rank, crank )
       call sll_s_int2string( sim%n_time_steps+sim%restart_steps, step )
       call sim%particle_group%group(1)%print(trim(sim%file_prefix)//step//'_particles_'//crank//'.dat')
    end if

  contains
    function bx_cossum_k(x)
      sll_real64             :: bx_cossum_k
      sll_real64, intent(in) :: x(3)

      bx_cossum_k = sim%bconstant(1) + sim%balpha(1) * cos( sum( sim%bkx * x ) - sim%bphase_shift(1) * sll_p_pi )

    end function bx_cossum_k

    function by_cossum_k(x)
      sll_real64             :: by_cossum_k
      sll_real64, intent(in) :: x(3)

      by_cossum_k = sim%bconstant(2) + sim%balpha(2) * cos( sum( sim%bkx * x ) - sim%bphase_shift(2) * sll_p_pi ) 

    end function by_cossum_k

    function bz_cossum_k(x)
      sll_real64             :: bz_cossum_k
      sll_real64, intent(in) :: x(3)

      bz_cossum_k = sim%bconstant(3) + sim%balpha(3) * cos( sum ( sim%bkx * x ) - sim%bphase_shift(3) * sll_p_pi )

    end function bz_cossum_k

    function bfield_sumcos_k(x)
      sll_real64             :: bfield_sumcos_k
      sll_real64, intent(in) :: x(3)

      bfield_sumcos_k = sim%bconstant(3) + sum( sim%balpha * cos(  sim%bkx * x - sim%bphase_shift * sll_p_pi ) ) 

    end function bfield_sumcos_k

    function bx_prodcos_k(x)
      sll_real64             :: bx_prodcos_k
      sll_real64, intent(in) :: x(3)

      bx_prodcos_k = sim%bconstant(1) + sim%balpha(1) * cos( sim%bkx(1) * x(1) - sim%bphase_shift(1) * sll_p_pi )* cos( sim%bkx(2) * x(2) - sim%bphase_shift(2) * sll_p_pi )* cos( sim%bkx(3) * x(3) - sim%bphase_shift(3) * sll_p_pi )

    end function bx_prodcos_k

    function by_prodcos_k(x)
      sll_real64             :: by_prodcos_k
      sll_real64, intent(in) :: x(3)

      by_prodcos_k = sim%bconstant(2) + sim%balpha(2) * cos( sim%bkx(1) * x(1) - sim%bphase_shift(1) * sll_p_pi )* cos( sim%bkx(2) * x(2) - sim%bphase_shift(2) * sll_p_pi )* cos( sim%bkx(3) * x(3) - sim%bphase_shift(3) * sll_p_pi )

    end function by_prodcos_k

    function bz_prodcos_k(x)
      sll_real64             :: bz_prodcos_k
      sll_real64, intent(in) :: x(3)

      bz_prodcos_k = sim%bconstant(3) + sim%balpha(3) * cos( sim%bkx(1) * x(1) - sim%bphase_shift(1) * sll_p_pi )* cos( sim%bkx(2) * x(2) - sim%bphase_shift(2) * sll_p_pi )* cos( sim%bkx(3) * x(3) - sim%bphase_shift(3) * sll_p_pi )

    end function bz_prodcos_k

    function bx_polynom_k(x)
      sll_real64             :: bx_polynom_k
      sll_real64, intent(in) :: x(3)

      bx_polynom_k = sim%bconstant(1) + sim%bkx(2) * x(2)

    end function bx_polynom_k

    function by_polynom_k(x)
      sll_real64             :: by_polynom_k
      sll_real64, intent(in) :: x(3)

      by_polynom_k = sim%bconstant(2) + sim%bkx(1) * x(1)

    end function by_polynom_k

    function bz_polynom_k(x)
      sll_real64             :: bz_polynom_k
      sll_real64, intent(in) :: x(3)

      bz_polynom_k = sim%bconstant(3) + sim%bkx(3) * x(3)

    end function bz_polynom_k

    function bfield_constant_k(x)
      sll_real64             :: bfield_constant_k
      sll_real64, intent(in) :: x(3)

      bfield_constant_k = sim%bconstant(3)

    end function bfield_constant_k

    function zero(x)
      sll_real64             :: zero
      sll_real64, intent(in) :: x(3)

      zero = 0._f64

    end function zero

  end subroutine run_pic_vm_3d3v


  !------------------------------------------------------------------------------!
  !> Finalize simulation
  subroutine delete_pic_vm_3d3v (sim)
    class(sll_t_sim_pic_vm_3d3v_cart), intent(inout) :: sim !< Singlespecies simulation

    SLL_ASSERT(storage_size(sim)>0)
    
    call sim%propagator%free()
    deallocate(sim%propagator)
    call sim%particle_group%group(1)%free()
    deallocate (sim%particle_group)
    call sim%maxwell_solver%free()
    !deallocate(sim%maxwell_solver)
    call sim%particle_mesh_coupling%free()
    deallocate(sim%particle_mesh_coupling)
    deallocate(sim%efield_dofs)
    deallocate(sim%bfield_dofs)
    if( sim%no_weights == 3) then
       call sim%control_variate%cv(1)%free()
    end if
    if ( sim%restart .eqv. .false. ) then
       call sim%init_distrib_params%free()
       deallocate(sim%init_distrib_params)
    end if
    call sim%sampler%free()

  end subroutine delete_pic_vm_3d3v


  !------------------------------------------------------------------------------!
  !> As a control variate, we use the equilibrium (v part of the initial distribution)
  function control_variate_equi( this, xi, vi, time) result(sll_f_control_variate)
    class(sll_t_control_variate) :: this !> control variate
    sll_real64, optional,  intent( in ) :: xi(:) !< particle position
    sll_real64, optional,  intent( in ) :: vi(:) !< particle velocity
    sll_real64, optional,  intent( in ) :: time  !< current time
    sll_real64               :: sll_f_control_variate

    sll_f_control_variate = &
         this%control_variate_distribution_params%eval_v_density( vi, xi, 1._f64 ) 

  end function control_variate_equi


  !> As a control variate, we use the equilibrium (v part of the initial distribution)
  function control_variate_xi( this, xi, vi, time) result(sll_f_control_variate)
    class(sll_t_control_variate) :: this !> control variate
    sll_real64, optional,  intent( in ) :: xi(:) !< particle position
    sll_real64, optional,  intent( in ) :: vi(:) !< particle velocity
    sll_real64, optional,  intent( in ) :: time  !< current time
    sll_real64               :: sll_f_control_variate

    sll_f_control_variate = &
         this%control_variate_distribution_params%eval_v_density( vi, xi, m=1._f64  )

  end function control_variate_xi


  !------------------------------------------------------------------------------!
  !Diagnostic functions and other helper functions
  !> Diagnostics for PIC Vlasov-Maxwell 3d3v 
  !> @todo (should be part of the library)
  subroutine sll_s_time_history_diagnostics_pic_vm_3d3v(&
       sim,&
       time, &
       file_id, &
       scratch1, scratch2)
    class(sll_t_sim_pic_vm_3d3v_cart), intent( inout ) :: sim !< Singlespecies simulation
    sll_real64,                        intent( in    ) :: time !< time
    sll_int32,                         intent( in    ) :: file_id !< file ide
    sll_real64,                        intent(   out ) :: scratch1(:) !< scratch data
    sll_real64,                        intent(   out ) :: scratch2(:) !< scratch data
    ! local variables
    sll_real64 :: diagnostics_local(6)
    sll_real64 :: diagnostics(6), phival
    sll_real64 :: potential_energy(6)
    sll_int32  :: i_part
    sll_real64 :: vi(3),  xi(3)
    sll_real64 :: wi(1)
    sll_real64 :: error_gauss

    diagnostics_local = 0.0_f64
    do i_part= 1, sim%particle_group%group(1)%n_particles
       vi = sim%particle_group%group(1)%get_v(i_part)
       xi = sim%particle_group%group(1)%get_x(i_part)
       wi = sim%particle_group%group(1)%get_mass(i_part)

       ! Kinetic energy
       diagnostics_local(1) = diagnostics_local(1) + &
            (vi(1)**2)*wi(1) !*0.5_f64
       diagnostics_local(2) = diagnostics_local(2) + &
            (vi(2)**2)*wi(1) !*0.5_f64
       diagnostics_local(3) = diagnostics_local(3) + &
            (vi(3)**2)*wi(1) !*0.5_f64
       ! Momentum 1
       diagnostics_local(4) = diagnostics_local(4) + &
            vi(1)*wi(1)
       ! Momentum 2
       diagnostics_local(5) = diagnostics_local(5) + &
            vi(2)*wi(1)
       ! Momentum 3
       diagnostics_local(6) = diagnostics_local(6) + &
            vi(3)*wi(1)
    end do
    diagnostics = 0.0_f64
    call sll_s_collective_reduce_real64(sll_v_world_collective, diagnostics_local, 6,&
         MPI_SUM, 0, diagnostics)
    diagnostics_local = 0._f64

    call compute_e_cross_b ( sim%maxwell_solver, scratch1, sim%efield_dofs, sim%bfield_dofs, diagnostics_local(4:6) )
    ! sum up the total momentum
    diagnostics(4:6) = diagnostics(4:6) + diagnostics_local(4:6)

    ! Check error in Gauss law
    call check_gauss_law ( sim, scratch1, scratch2, error_gauss )

    if (sim%rank == 0) then
       if(sim%ct) then
          potential_energy(1) = sim%maxwell_solver%inner_product &
               ( sim%efield_dofs, sim%efield_dofs, 1, 1 )/sim%plasma_betar(2)
          potential_energy(2) = sim%maxwell_solver%inner_product &
               ( sim%efield_dofs, sim%efield_dofs, 1, 2 )/sim%plasma_betar(2)
          potential_energy(3) = sim%maxwell_solver%inner_product &
               ( sim%efield_dofs, sim%efield_dofs, 1, 3 )/sim%plasma_betar(2)
          potential_energy(4) = sim%maxwell_solver%inner_product &
               ( sim%bfield_dofs, sim%bfield_dofs, 2, 1 )*sim%plasma_betar(3)
          potential_energy(5) = sim%maxwell_solver%inner_product &
               ( sim%bfield_dofs, sim%bfield_dofs, 2, 2 )*sim%plasma_betar(3)
          potential_energy(6) = sim%maxwell_solver%inner_product &
               ( sim%bfield_dofs, sim%bfield_dofs, 2, 3 )*sim%plasma_betar(3)
       else
          potential_energy(1) = sim%maxwell_solver%inner_product &
               ( sim%efield_dofs(1:sim%n_totaldofs1), sim%efield_dofs(1:sim%n_totaldofs1), 1, 1 )/sim%plasma_betar(2)
          potential_energy(2) = sim%maxwell_solver%inner_product &
               ( sim%efield_dofs(sim%n_totaldofs1+1:sim%n_totaldofs1+sim%n_totaldofs0), sim%efield_dofs(sim%n_totaldofs1+1:sim%n_totaldofs1+sim%n_totaldofs0), 1, 2 )/sim%plasma_betar(2)
          potential_energy(3) = sim%maxwell_solver%inner_product &
               ( sim%efield_dofs(sim%n_totaldofs1+sim%n_totaldofs0+1:sim%n_totaldofs1+sim%n_totaldofs0*2), sim%efield_dofs(sim%n_totaldofs1+sim%n_totaldofs0+1:sim%n_totaldofs1+sim%n_totaldofs0*2), 1, 3 )/sim%plasma_betar(2)
          potential_energy(4) = sim%maxwell_solver%inner_product &
               ( sim%bfield_dofs(1:sim%n_totaldofs0), sim%bfield_dofs(1:sim%n_totaldofs0), 2, 1 )*sim%plasma_betar(3)
          potential_energy(5) = sim%maxwell_solver%inner_product &
               ( sim%bfield_dofs(sim%n_totaldofs0+1:sim%n_totaldofs0+sim%n_totaldofs1), sim%bfield_dofs(sim%n_totaldofs0+1:sim%n_totaldofs0+sim%n_totaldofs1), 2, 2 )*sim%plasma_betar(3)
          potential_energy(6) =sim%maxwell_solver%inner_product &
               ( sim%bfield_dofs(sim%n_totaldofs0+sim%n_totaldofs1+1:sim%n_totaldofs0+sim%n_totaldofs1*2), sim%bfield_dofs(sim%n_totaldofs0+sim%n_totaldofs1+1:sim%n_totaldofs0+sim%n_totaldofs1*2), 2, 3 )*sim%plasma_betar(3)
       end if

       if( sim%adiabatic_electrons) then
          phival = sim%maxwell_solver%l2norm_squared( sim%phi_dofs, 0, 0 ) 
          write(file_id,'(f12.5,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16)' ) &
               time,  potential_energy, diagnostics(1:3), &
               sum(diagnostics(1:3)) + phival, diagnostics(4:6), &
               diagnostics_local(4:6), error_gauss, phival
       else
          write(file_id,'(f12.5,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16)' ) &
               time,  potential_energy, diagnostics(1:3), &
               sum(diagnostics(1:3)) +  sim%force_sign*sum(potential_energy), diagnostics(4:6), &
               diagnostics_local(4:6), error_gauss
       end if
    end if

  end subroutine sll_s_time_history_diagnostics_pic_vm_3d3v


  !> Accumulate rho and solve Poisson
  subroutine solve_poisson( sim, rho_local, rho )
    class(sll_t_sim_pic_vm_3d3v_cart), intent( inout ) :: sim !< Singlespecies simulation
    sll_real64,                        intent(   out ) :: rho_local(:) !< charge local to mpi processors
    sll_real64,                        intent(   out ) :: rho(:) !< charge
    !local variables
    sll_int32 :: i_part
    sll_real64 :: xi(3), wi(1), scratch(1:sim%n_totaldofs1+2*sim%n_totaldofs0)

    rho_local = 0.0_f64
    do i_part = 1, sim%particle_group%group(1)%n_particles
       xi = sim%particle_group%group(1)%get_x(i_part)
       ! Get charge for accumulation of rho
       wi(1) = sim%particle_group%group(1)%get_charge( i_part, sim%no_weights)
       call sim%particle_mesh_coupling%add_charge(xi, wi(1), &
            sim%degree_smoother, rho_local)
    end do
    ! MPI to sum up contributions from each processor
    rho = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, &
         rho_local, sim%n_totaldofs0, MPI_SUM, rho)

    call sim%filter%apply_inplace( rho )

    if (sim%no_weights == 1) then
       rho_local = 0._f64
       ! Add background distribution of neutralizing ions
       call add_background_charge(sim, rho_local)
       rho =  sim%force_sign*(rho - sim%charge * rho_local) 
    end if

    if( sim%adiabatic_electrons) then
       call sim%maxwell_solver%compute_phi_from_rho( rho, sim%phi_dofs, sim%efield_dofs)
    else
       rho = rho * sim%plasma_betar(2)
       call sim%maxwell_solver%compute_E_from_rho( rho, sim%efield_dofs )
    end if

    if( sim%external_field) then
       call sim%maxwell_solver%L2projection( 1, 0, scratch, E_x, E_y, zero )
       sim%efield_dofs = sim%efield_dofs + scratch
    end if

  contains
    function sin_1k(x)
      sll_real64             :: sin_1k
      sll_real64, intent(in) :: x(3)
      sll_int32  :: j
      select type ( q=>sim%init_distrib_params )
      type is ( sll_t_params_cos_gaussian )
         sin_1k = 0._f64
         do j = 1, q%n_cos
            sin_1k = sin_1k - q%alpha(j)/q%kx(1,j)*sin( sum(q%kx(:,j) * x) - q%phase_shift(j)*sll_p_pi )
         end do
      end select
    end function sin_1k

    function E_x(x)
      sll_real64             :: E_x
      sll_real64, intent(in) :: x(3)

      E_x = -x(1)
    end function E_x

    function E_y(x)
      sll_real64             :: E_y
      sll_real64, intent(in) :: x(3)

      E_y = -x(2)
    end function E_y

    function zero(x)
      sll_real64             :: zero
      sll_real64, intent(in) :: x(3)

      zero = 0._f64
    end function zero

  end subroutine solve_poisson


  !> Add background charge
  subroutine add_background_charge(sim, background_charge)
    class(sll_t_sim_pic_vm_3d3v_cart), intent(inout) :: sim !< Singlespecies simulation
    sll_real64,                        intent(  out) :: background_charge(sim%n_totaldofs0) !< Coefficient vector of the charge distribution
    !local variables
    sll_int32 :: i1, i2, i3, c1, c2, c3, k1, k2, k3, l1, l2, l3
    sll_int32 :: index1d, index1, index2, index3, q(3)
    sll_real64, allocatable :: quad_xw1(:,:), quad_xw2(:,:), quad_xw3(:,:)
    sll_real64, allocatable :: bspl1(:,:), bspl2(:,:), bspl3(:,:)
    sll_real64 :: jacobian, p(3)

    background_charge = 0._f64
    q = sim%degree_smoother+1
    allocate(quad_xw1(2, q(1)))
    allocate(bspl1(sim%degree_smoother(1)+1, q(1) ))
    bspl1 = 0._f64
    quad_xw1 = sll_f_gauss_legendre_points_and_weights( q(1) , 0.0_f64, 1.0_f64 )
    do i1 = 1, q(1)
       call sll_s_uniform_bsplines_eval_basis(sim%degree_smoother(1), quad_xw1(1,i1), bspl1(:,i1) )
    end do

    allocate(quad_xw2(2, q(2)))
    allocate(bspl2(sim%degree_smoother(2)+1, q(2) ))
    bspl2 = 0._f64
    quad_xw2 = sll_f_gauss_legendre_points_and_weights( q(2) , 0.0_f64, 1.0_f64 )
    do i2 = 1, q(2)
       call sll_s_uniform_bsplines_eval_basis(sim%degree_smoother(2), quad_xw2(1,i2), bspl2(:,i2) )
    end do

    allocate(quad_xw3(2, q(3)))
    allocate(bspl3(sim%degree_smoother(3)+1, q(3) ))
    bspl3 = 0._f64
    quad_xw3 = sll_f_gauss_legendre_points_and_weights( q(3) , 0.0_f64, 1.0_f64 )
    do i3 = 1, q(3)
       call sll_s_uniform_bsplines_eval_basis(sim%degree_smoother(3), quad_xw3(1,i3), bspl3(:,i3) )
    end do

    if( sim%ct) then
       if(sim%boundary) then
          if( sim%adiabatic_electrons ) then
             do i3 = 1, sim%n_gcells(3)
                do i2 = 1, sim%n_gcells(2)
                   do i1 = 1, sim%degree_smoother(1)-1
                      do k3 = 1, q(3)
                         p(3) = (quad_xw3(1, k3) + real(i3-1, f64))/real(sim%n_gcells(3), f64)
                         do k2 = 1, q(2)
                            p(2) = (quad_xw2(1, k2) + real(i2-1, f64))/real(sim%n_gcells(2), f64)
                            do k1 = 1, q(1)
                               p(1) = (quad_xw1(1, k1) + real(i1-1, f64))/real(sim%n_gcells(1), f64)
                               jacobian= sim%map%jacobian( p ) * sim%profile%rho_0( p(1) )
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1= i1+c1-1
                                        index1d = index1 + index2*(sim%n_gcells(1)+sim%degree_smoother(1)) + index3*(sim%n_gcells(1)+sim%degree_smoother(1))*sim%n_gcells(2)

                                        background_charge(index1d) = background_charge(index1d) + &
                                             sll_f_spline_pp_horner_1d( sim%degree_smoother(1), sim%particle_mesh_coupling%spline_pp_0%spline1%poly_coeffs_boundary_left(:,:,i1), quad_xw1(1,k1), c1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) * jacobian
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                   do i1 = sim%degree_smoother(1), sim%n_gcells(1)-sim%degree_smoother(1)+1
                      do k3 = 1, q(3)
                         p(3) = (quad_xw3(1, k3) + real(i3-1, f64))/real(sim%n_gcells(3), f64)
                         do k2 = 1, q(2)
                            p(2) = (quad_xw2(1, k2) + real(i2-1, f64))/real(sim%n_gcells(2), f64)
                            do k1 = 1, q(1)
                               p(1) = (quad_xw1(1, k1) + real(i1-1, f64))/real(sim%n_gcells(1), f64)
                               jacobian= sim%map%jacobian( p ) * sim%profile%rho_0( p(1) )
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1= i1+c1-1
                                        index1d = index1 + index2*(sim%n_gcells(1)+sim%degree_smoother(1)) + index3*(sim%n_gcells(1)+sim%degree_smoother(1))*sim%n_gcells(2)

                                        background_charge(index1d) = background_charge(index1d) + &
                                             bspl1(c1,k1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) * jacobian
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                   do i1 = sim%n_gcells(1)-sim%degree_smoother(1)+2, sim%n_gcells(1)
                      do k3 = 1, q(3)
                         p(3) = (quad_xw3(1, k3) + real(i3-1, f64))/real(sim%n_gcells(3), f64)
                         do k2 = 1, q(2)
                            p(2) = (quad_xw2(1, k2) + real(i2-1, f64))/real(sim%n_gcells(2), f64)
                            do k1 = 1, q(1)
                               p(1) = (quad_xw1(1, k1) + real(i1-1, f64))/real(sim%n_gcells(1), f64)
                               jacobian= sim%map%jacobian( p ) * sim%profile%rho_0( p(1) ) 
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1= i1+c1-1
                                        index1d = index1 + index2*(sim%n_gcells(1)+sim%degree_smoother(1)) + index3*(sim%n_gcells(1)+sim%degree_smoother(1))*sim%n_gcells(2)

                                        background_charge(index1d) = background_charge(index1d) + &
                                             sll_f_spline_pp_horner_1d( sim%degree_smoother(1), sim%particle_mesh_coupling%spline_pp_0%spline1%poly_coeffs_boundary_right(:,:,i1-sim%n_gcells(1)+sim%degree_smoother(1)-1), quad_xw1(1,k1), c1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) * jacobian
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
             background_charge = background_charge/real(sim%n_totalcells,f64)
          else
             do i3 = 1, sim%n_gcells(3)
                do i2 = 1, sim%n_gcells(2)
                   do i1 = 1, sim%degree_smoother(1)-1
                      do k3 = 1, q(3)
                         p(3) = (quad_xw3(1, k3) + real(i3-1, f64))/real(sim%n_gcells(3), f64)
                         do k2 = 1, q(2)
                            p(2) = (quad_xw2(1, k2) + real(i2-1, f64))/real(sim%n_gcells(2), f64)
                            do k1 = 1, q(1)
                               p(1) = (quad_xw1(1, k1) + real(i1-1, f64))/real(sim%n_gcells(1), f64)
                               jacobian= sim%map%jacobian( p ) 
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1= i1+c1-1
                                        index1d = index1 + index2*(sim%n_gcells(1)+sim%degree_smoother(1)) + index3*(sim%n_gcells(1)+sim%degree_smoother(1))*sim%n_gcells(2)

                                        background_charge(index1d) = background_charge(index1d) + &
                                             sll_f_spline_pp_horner_1d( sim%degree_smoother(1), sim%particle_mesh_coupling%spline_pp_0%spline1%poly_coeffs_boundary_left(:,:,i1), quad_xw1(1,k1), c1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) * jacobian
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                   do i1 = sim%degree_smoother(1), sim%n_gcells(1)-sim%degree_smoother(1)+1
                      do k3 = 1, q(3)
                         p(3) = (quad_xw3(1, k3) + real(i3-1, f64))/real(sim%n_gcells(3), f64)
                         do k2 = 1, q(2)
                            p(2) = (quad_xw2(1, k2) + real(i2-1, f64))/real(sim%n_gcells(2), f64)
                            do k1 = 1, q(1)
                               p(1) = (quad_xw1(1, k1) + real(i1-1, f64))/real(sim%n_gcells(1), f64)
                               jacobian= sim%map%jacobian( p ) 
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1= i1+c1-1
                                        index1d = index1 + index2*(sim%n_gcells(1)+sim%degree_smoother(1)) + index3*(sim%n_gcells(1)+sim%degree_smoother(1))*sim%n_gcells(2)

                                        background_charge(index1d) = background_charge(index1d) + &
                                             bspl1(c1,k1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) * jacobian
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                   do i1 = sim%n_gcells(1)-sim%degree_smoother(1)+2, sim%n_gcells(1)
                      do k3 = 1, q(3)
                         p(3) = (quad_xw3(1, k3) + real(i3-1, f64))/real(sim%n_gcells(3), f64)
                         do k2 = 1, q(2)
                            p(2) = (quad_xw2(1, k2) + real(i2-1, f64))/real(sim%n_gcells(2), f64)
                            do k1 = 1, q(1)
                               p(1) = (quad_xw1(1, k1) + real(i1-1, f64))/real(sim%n_gcells(1), f64)
                               jacobian= sim%map%jacobian( p )  
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1= i1+c1-1
                                        index1d = index1 + index2*(sim%n_gcells(1)+sim%degree_smoother(1)) + index3*(sim%n_gcells(1)+sim%degree_smoother(1))*sim%n_gcells(2)

                                        background_charge(index1d) = background_charge(index1d) + &
                                             sll_f_spline_pp_horner_1d( sim%degree_smoother(1), sim%particle_mesh_coupling%spline_pp_0%spline1%poly_coeffs_boundary_right(:,:,i1-sim%n_gcells(1)+sim%degree_smoother(1)-1), quad_xw1(1,k1), c1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) * jacobian
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
             background_charge = background_charge/real(sim%n_totalcells,f64)
          end if
       else
          if( sim%adiabatic_electrons ) then
             do i3 = 1, sim%n_gcells(3)
                do i2 = 1, sim%n_gcells(2)
                   do i1 = 1, sim%n_gcells(1)
                      do l3 = 1, q(3)
                         k3 = 1 - (-1)**l3 * l3/2 + modulo(l3+1,2)*q(3)
                         p(3) = (quad_xw3(1, k3) + real(i3-1, f64))/real(sim%n_gcells(3), f64)
                         do l2 = 1, q(2)
                            k2 = 1 - (-1)**l2 * l2/2 + modulo(l2+1,2)*q(2)
                            p(2) = (quad_xw2(1, k2) + real(i2-1, f64))/real(sim%n_gcells(2), f64)
                            do l1 = 1, q(1)
                               k1 = 1 - (-1)**l1 * l1/2 + modulo(l1+1,2)*q(1)
                               p(1) = (quad_xw1(1, k1) + real(i1-1, f64))/real(sim%n_gcells(1), f64)
                               jacobian= sim%map%jacobian( p ) * sim%profile%rho_0( p(1) ) 
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1 = modulo(i1-sim%degree_smoother(1)+c1-2, sim%n_gcells(1))
                                        index1d = index1+1 + index2*sim%n_gcells(1) + index3*sim%n_gcells(1)*sim%n_gcells(2)
                                        background_charge(index1d) = background_charge(index1d) + &
                                             bspl1(c1,k1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) * jacobian
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
             background_charge = background_charge/real(sim%n_totalcells,f64)
          else
             do i3 = 1, sim%n_gcells(3)
                do i2 = 1, sim%n_gcells(2)
                   do i1 = 1, sim%n_gcells(1)
                      do l3 = 1, q(3)
                         k3 = 1 - (-1)**l3 * l3/2 + modulo(l3+1,2)*q(3)
                         p(3) = (quad_xw3(1, k3) + real(i3-1, f64))/real(sim%n_gcells(3), f64)
                         do l2 = 1, q(2)
                            k2 = 1 - (-1)**l2 * l2/2 + modulo(l2+1,2)*q(2)
                            p(2) = (quad_xw2(1, k2) + real(i2-1, f64))/real(sim%n_gcells(2), f64)
                            do l1 = 1, q(1)
                               k1 = 1 - (-1)**l1 * l1/2 + modulo(l1+1,2)*q(1)
                               p(1) = (quad_xw1(1, k1) + real(i1-1, f64))/real(sim%n_gcells(1), f64)
                               jacobian= sim%map%jacobian( p )
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1 = modulo(i1-sim%degree_smoother(1)+c1-2, sim%n_gcells(1))
                                        index1d = index1+1 + index2*sim%n_gcells(1) + index3*sim%n_gcells(1)*sim%n_gcells(2)
                                        background_charge(index1d) = background_charge(index1d) + &
                                             bspl1(c1,k1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) * jacobian
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
             background_charge = background_charge/real(sim%n_totalcells,f64)
          end if
       end if
    else
       if(sim%boundary) then
          if( sim%adiabatic_electrons ) then
             do i3 = 1, sim%n_gcells(3)
                do i2 = 1, sim%n_gcells(2)
                   do i1 = 1, sim%degree_smoother(1)-1
                      do k3 = 1, q(3)
                         do k2 = 1, q(2)
                            do k1 = 1, q(1)
                               p(1) = (quad_xw1(1, k1) + real(i1-1, f64))/real(sim%n_gcells(1), f64)
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1= i1+c1-1
                                        index1d = index1 + index2*(sim%n_gcells(1)+sim%degree_smoother(1)) + index3*(sim%n_gcells(1)+sim%degree_smoother(1))*sim%n_gcells(2)
                                        background_charge(index1d) = background_charge(index1d) + &
                                             sll_f_spline_pp_horner_1d( sim%degree_smoother(1), sim%particle_mesh_coupling%spline_pp_0%spline1%poly_coeffs_boundary_left(:,:,i1), quad_xw1(1,k1), c1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) * sim%profile%rho_0( p(1) ) 
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                   do i1 = sim%degree_smoother(1), sim%n_gcells(1)-sim%degree_smoother(1)+1
                      do k3 = 1, q(3)
                         do k2 = 1, q(2)
                            do k1 = 1, q(1)
                               p(1) = (quad_xw1(1, k1) + real(i1-1, f64))/real(sim%n_gcells(1), f64)
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1= i1+c1-1
                                        index1d = index1 + index2*(sim%n_gcells(1)+sim%degree_smoother(1)) + index3*(sim%n_gcells(1)+sim%degree_smoother(1))*sim%n_gcells(2)
                                        background_charge(index1d) = background_charge(index1d) + &
                                             bspl1(c1,k1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) * sim%profile%rho_0( p(1) ) 
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                   do i1 = sim%n_gcells(1)-sim%degree_smoother(1)+2, sim%n_gcells(1)
                      do k3 = 1, q(3)
                         do k2 = 1, q(2)
                            do k1 = 1, q(1)
                               p(1) = (quad_xw1(1, k1) + real(i1-1, f64))/real(sim%n_gcells(1), f64)
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1= i1+c1-1
                                        index1d = index1 + index2*(sim%n_gcells(1)+sim%degree_smoother(1)) + index3*(sim%n_gcells(1)+sim%degree_smoother(1))*sim%n_gcells(2)
                                        background_charge(index1d) = background_charge(index1d) + &
                                             sll_f_spline_pp_horner_1d( sim%degree_smoother(1), sim%particle_mesh_coupling%spline_pp_0%spline1%poly_coeffs_boundary_right(:,:,i1-sim%n_gcells(1)+sim%degree_smoother(1)-1), quad_xw1(1,k1), c1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) * sim%profile%rho_0( p(1) ) 
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
             background_charge = background_charge * sim%volume
          else
             do i3 = 1, sim%n_gcells(3)
                do i2 = 1, sim%n_gcells(2)
                   do i1 = 1, sim%degree_smoother(1)-1
                      do k3 = 1, q(3)
                         do k2 = 1, q(2)
                            do k1 = 1, q(1)
                               p(1) = sim%delta_x(1)*(quad_xw1(1, k1) + real(i1-1,f64) ) 
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1= i1+c1-1
                                        index1d = index1 + index2*(sim%n_gcells(1)+sim%degree_smoother(1)) + index3*(sim%n_gcells(1)+sim%degree_smoother(1))*sim%n_gcells(2)
                                        background_charge(index1d) = background_charge(index1d) + &
                                             sll_f_spline_pp_horner_1d( sim%degree_smoother(1), sim%particle_mesh_coupling%spline_pp_0%spline1%poly_coeffs_boundary_left(:,:,i1), quad_xw1(1,k1), c1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3)
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                   do i1 = sim%degree_smoother(1), sim%n_gcells(1)-sim%degree_smoother(1)+1
                      do k3 = 1, q(3)
                         do k2 = 1, q(2)
                            do k1 = 1, q(1)
                               p(1) = sim%delta_x(1)*(quad_xw1(1, k1) + real(i1-1,f64) ) 
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1= i1+c1-1
                                        index1d = index1 + index2*(sim%n_gcells(1)+sim%degree_smoother(1)) + index3*(sim%n_gcells(1)+sim%degree_smoother(1))*sim%n_gcells(2)
                                        background_charge(index1d) = background_charge(index1d) + &
                                             bspl1(c1,k1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) 
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                   do i1 = sim%n_gcells(1)-sim%degree_smoother(1)+2, sim%n_gcells(1)
                      do k3 = 1, q(3)
                         do k2 = 1, q(2)
                            do k1 = 1, q(1)
                               p(1) = sim%delta_x(1)*(quad_xw1(1, k1) + real(i1-1,f64) ) 
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1= i1+c1-1
                                        index1d = index1 + index2*(sim%n_gcells(1)+sim%degree_smoother(1)) + index3*(sim%n_gcells(1)+sim%degree_smoother(1))*sim%n_gcells(2)
                                        background_charge(index1d) = background_charge(index1d) + &
                                             sll_f_spline_pp_horner_1d( sim%degree_smoother(1), sim%particle_mesh_coupling%spline_pp_0%spline1%poly_coeffs_boundary_right(:,:,i1-sim%n_gcells(1)+sim%degree_smoother(1)-1), quad_xw1(1,k1), c1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) 
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
             background_charge = background_charge * sim%volume
          end if
       else
          if( sim%adiabatic_electrons ) then
             do i3 = 1, sim%n_gcells(3)
                do i2 = 1, sim%n_gcells(2)
                   do i1 = 1, sim%n_gcells(1)
                      do l3 = 1, q(3)
                         k3 = 1 - (-1)**l3 * l3/2 + modulo(l3+1,2)*q(3)
                         do l2 = 1, q(2)
                            k2 = 1 - (-1)**l2 * l2/2 + modulo(l2+1,2)*q(2)
                            do l1 = 1, q(1)
                               k1 = 1 - (-1)**l1 * l1/2 + modulo(l1+1,2)*q(1)
                               p(1) = (quad_xw1(1, k1) + real(i1-1, f64))/real(sim%n_gcells(1), f64)
                               do c3 = 1,sim%degree_smoother(3)+1
                                  index3 = modulo(i3-sim%degree_smoother(3)+c3-2, sim%n_gcells(3))
                                  do c2 = 1,sim%degree_smoother(2)+1
                                     index2 = modulo(i2-sim%degree_smoother(2)+c2-2, sim%n_gcells(2))
                                     do c1 = 1,sim%degree_smoother(1)+1
                                        index1 = modulo(i1-sim%degree_smoother(1)+c1-2, sim%n_gcells(1))
                                        index1d = index1+1 + index2*sim%n_gcells(1) + index3*sim%n_gcells(1)*sim%n_gcells(2)
                                        background_charge(index1d) = background_charge(index1d) + &
                                             bspl1(c1,k1) * bspl2(c2,k2) * bspl3(c3,k3) * &
                                             quad_xw1(2, k1) * quad_xw2(2, k2) * quad_xw3(2, k3) *  sim%profile%rho_0( p(1) )
                                     end do
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
             background_charge = background_charge* sim%volume
          else
             background_charge = sim%volume
          end if
       end if
    end if

  end subroutine add_background_charge


  !> check Gauss' law
  subroutine check_gauss_law (sim, rho_gauss, rho, error )
    class(sll_t_sim_pic_vm_3d3v_cart), intent( inout ) :: sim !< Singlespecies simulation
    sll_real64,                        intent(   out ) :: rho_gauss(:) !< charge local to mpi processes
    sll_real64,                        intent(   out ) :: rho(:) !< charge
    sll_real64,                        intent(   out ) :: error !< error in Gauss' law
    !local variables
    sll_int32 :: i_part
    sll_real64 :: xi(3), wi(1)

    rho_gauss = 0.0_f64
    do i_part = 1, sim%particle_group%group(1)%n_particles
       xi = sim%particle_group%group(1)%get_x(i_part)
       if (xi(1) >= sim%domain(1,1) .and. xi(1) <= sim%domain(1,2)) then
          ! Get charge for accumulation of rho
          wi(1) = sim%particle_group%group(1)%get_charge( i_part, sim%no_weights)
          call sim%particle_mesh_coupling%add_charge(xi, wi(1), &
               sim%degree_smoother, rho_gauss)
       end if
    end do

    ! MPI to sum up contributions from each processor
    rho = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, &
         rho_gauss, sim%n_totaldofs0, MPI_SUM, rho)
    rho_gauss = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, &
         sim%rhob, sim%n_totaldofs0, MPI_SUM, rho_gauss)
    rho = rho + rho_gauss

    call sim%filter%apply_inplace( rho )

    if (sim%no_weights == 1) then
       ! Add neutralizing background distribution 
       call add_background_charge(sim, rho_gauss)
       rho = sim%force_sign*(rho - sim%charge * rho_gauss)
    end if

    ! Solve Gauss law
    if(sim%adiabatic_electrons) then
       call sim%maxwell_solver%multiply_mass( [0], sim%phi_dofs, rho_gauss )
    else
       if( .not. sim%boundary ) then
          error = 0._f64
          do i_part = 1, sim%n_totaldofs0
             error = error + rho(i_part)
          end do
          error = error/real(sim%n_totaldofs0,f64)
          rho = rho - error
          rho_gauss=0._f64
       end if

       rho = rho * sim%plasma_betar(2)
       call sim%maxwell_solver%compute_rho_from_e( sim%efield_dofs, rho_gauss )
    end if
    error = maxval(abs(rho - rho_gauss ))

  end subroutine check_gauss_law


  !> Compute ExB
  subroutine compute_e_cross_b ( maxwell, scratch, efield, bfield, e_cross_b )
    class(sll_c_maxwell_3d_base), intent(inout) :: maxwell !> Maxwell solver
    sll_real64, intent( out ) :: scratch(:) !> scratch data
    sll_real64, intent( in  ) :: efield(:) !> E
    sll_real64, intent( in  ) :: bfield(:) !> B
    sll_real64, intent( out ) :: e_cross_b(3)  !< E cross B

    call  maxwell%multiply_mass( [3, 2, 1], bfield(maxwell%n_total0+maxwell%n_total1+1:maxwell%n_total0+maxwell%n_total1*2), scratch )
    e_cross_b(1) = sum(efield(maxwell%n_total1+1:maxwell%n_total1+maxwell%n_total0)*scratch)
    call  maxwell%multiply_mass([3, 1, 2], bfield(maxwell%n_total0+1:maxwell%n_total0+maxwell%n_total1), scratch )
    e_cross_b(1) = e_cross_b(1) - sum(efield(maxwell%n_total1+maxwell%n_total0+1:maxwell%n_total1+maxwell%n_total0*2)*scratch)

    call maxwell%multiply_mass( [1, 3, 2], bfield(1:maxwell%n_total0), scratch )
    e_cross_b(2) = sum(efield(maxwell%n_total1+maxwell%n_total0+1:maxwell%n_total1+maxwell%n_total0*2)*scratch)
    call maxwell%multiply_mass( [2, 3, 1], bfield(maxwell%n_total0+maxwell%n_total1+1:maxwell%n_total0+maxwell%n_total1*2), scratch(1:maxwell%n_total1) )
    e_cross_b(2) = e_cross_b(2) - sum(efield(1:maxwell%n_total1)*scratch(1:maxwell%n_total1))


    call  maxwell%multiply_mass( [2, 1, 3], bfield(maxwell%n_total0+1:maxwell%n_total0+maxwell%n_total1), scratch(1:maxwell%n_total1) )
    e_cross_b(3) = sum(efield(1:maxwell%n_total1)*scratch(1:maxwell%n_total1))
    call  maxwell%multiply_mass( [1, 2, 3], bfield(1:maxwell%n_total0), scratch )
    e_cross_b(3) = e_cross_b(3) - sum(efield(maxwell%n_total1+1:maxwell%n_total1+maxwell%n_total0)*scratch)

  end subroutine compute_e_cross_b


  !> Compute ExB with coordinate transformation
  subroutine compute_e_cross_b_curvilinear( particle_mesh_coupling, deg, map, efield_dofs, bfield_dofs, n_total0, n_total1, ecb)
    class(sll_c_particle_mesh_coupling_3d), intent(inout) :: particle_mesh_coupling !> Particle mesh coupling
    sll_int32,  intent( in    )             :: deg(3)     !< maximal spline deg
    type(sll_t_mapping_3d), intent( inout    ) :: map        !< coordinate transformation
    sll_real64, intent( in    )             :: efield_dofs(:) !< DoFs describing the two components of the electric field
    sll_real64, intent( in    )             :: bfield_dofs(:)   !< DoFs describing the magnetic field
    sll_int32,  intent( in    )             :: n_total0, n_total1 !< total number of DoFs for 0- and 1-form
    sll_real64, intent(   out )             :: ecb(3) !< E cross B 
    !local variables
    sll_int32  :: i, k3, k2, k1, j3, j2, j1, q(3)
    sll_real64, allocatable :: xw_gauss_d1(:,:), xw_gauss_d2(:,:), xw_gauss_d3(:,:)
    sll_real64 :: xi(3), N(3,3), DF(3,3)
    sll_real64 :: e_phys(3), b_phys(3), efield(3), bfield(3)

    efield = 0.0_f64
    bfield = 0.0_f64
    e_phys = 0.0_f64
    b_phys = 0.0_f64
    ecb = 0.0_f64

    q = 2*maxval(deg)+1
    allocate( xw_gauss_d1(1:2, 1:q(1)) )
    allocate( xw_gauss_d2(1:2, 1:q(2)) )
    allocate( xw_gauss_d3(1:2, 1:q(3)) )

    xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights( q(1), 0._f64, 1._f64 )
    xw_gauss_d2 = sll_f_gauss_legendre_points_and_weights( q(2), 0._f64, 1._f64 )
    xw_gauss_d3 = sll_f_gauss_legendre_points_and_weights( q(3), 0._f64, 1._f64 )

    do j3 = 1, particle_mesh_coupling%n_cells(3) 
       do j2 = 1, particle_mesh_coupling%n_cells(2)
          do j1 = 1, particle_mesh_coupling%n_cells(1)
             ! loop over Gauss points
             do k3 = 1, q(3) 
                xi(3) = particle_mesh_coupling%delta_x(3)*(xw_gauss_d3(1,k3) + real(j3 - 1,f64))
                do k2 = 1, q(2)
                   xi(2) = particle_mesh_coupling%delta_x(2)*(xw_gauss_d2(1,k2) + real(j2 - 1,f64))
                   do k1 = 1, q(1)
                      xi(1) = particle_mesh_coupling%delta_x(1)*(xw_gauss_d1(1,k1) + real(j1 - 1,f64))
                      
                      call particle_mesh_coupling%evaluate( xi, [deg(1)-1, deg(2), deg(3)], &
                           efield_dofs(1:n_total1), efield(1))
                      call particle_mesh_coupling%evaluate( xi, [deg(1), deg(2)-1, deg(3)], &
                           efield_dofs(1+n_total1:n_total1+n_total0), efield(2))
                      call particle_mesh_coupling%evaluate( xi, [deg(1), deg(2), deg(3)-1], &
                           efield_dofs(1+n_total1+n_total0:n_total1+2*n_total0), efield(3))

                      call particle_mesh_coupling%evaluate( xi, [deg(1), deg(2)-1, deg(3)-1], &
                           bfield_dofs(1:n_total0), bfield(1))
                      call particle_mesh_coupling%evaluate( xi, [deg(1)-1, deg(2), deg(3)-1], &
                           bfield_dofs(1+n_total0:n_total0+n_total1), bfield(2))
                      call particle_mesh_coupling%evaluate( xi, [deg(1)-1, deg(2)-1, deg(3)], &
                           bfield_dofs(1+n_total0+n_total1:n_total0+2*n_total1), bfield(3))

                      N = map%jacobian_matrix_inverse( xi )
                      DF = map%jacobian_matrix( xi ) 

                      do i = 1, 3
                         e_phys(i) = (N(i,1)* efield(1)+N(i,2)* efield(2)+N(i,3)* efield(3)) *&
                              xw_gauss_d1(2,k1)* xw_gauss_d2(2,k2)* xw_gauss_d3(2,k3)
                         b_phys(i) = (DF(i,1)* bfield(1)+DF(i,2)* bfield(2)+DF(i,3)* bfield(3)) *&
                              xw_gauss_d1(2,k1)* xw_gauss_d2(2,k2)* xw_gauss_d3(2,k3)
                      end do
                      
                      ecb(1) = ecb(1) + e_phys(2) * b_phys(3) - e_phys(3) * b_phys(2)
                      ecb(2) = ecb(2) + e_phys(3) * b_phys(1) - e_phys(1) * b_phys(3)
                      ecb(3) = ecb(3) + e_phys(1) * b_phys(2) - e_phys(2) * b_phys(1)
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine compute_e_cross_b_curvilinear


  !> Check diagnostics
  subroutine sll_s_check_diagnostics(reffile, simfile, tol_error, passed)
    character(*), intent(in) :: reffile !< Name of reference file (stored in same folder as source file)
    character(*), intent(in) :: simfile !< Name of file with simulation result
    sll_real64, intent(in)   :: tol_error !< tolerance 
    logical, intent(out)     :: passed !< true if diagnostics checks out
    !local variables
    sll_real64 :: error
    sll_real64 :: data_sim(4,18)
    sll_real64 :: data_ref(4,18)
    sll_int32  :: file_id

    ! Read simulation result
    open(newunit=file_id, file=simfile, status='old', action='read')
    read(unit=file_id,fmt=*) data_sim
    close(file_id)

    ! Read reference
    open(newunit=file_id, file=reffile, status='old', action='read')
    read(unit=file_id,fmt=*) data_ref
    close(file_id)

    ! Compare
    data_sim = data_sim - data_ref
    error = maxval(abs(data_sim))
    print*, 'Max error in time history diagnostics: ', error
    if (error < tol_error) then
       passed = .true.
    else
       passed = .false.
    end if

  end subroutine sll_s_check_diagnostics


end module sll_m_sim_pic_vm_3d3v_cart
