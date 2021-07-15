! Simulation of 1d2v Vlasov-Maxwell with simple PIC method, periodic boundary conditions, Weibel instability. FEM with splines, degree 3 for B and 2 for E

! Species 1: electrons, species 2: ions

! author: Katharina Kormann, IPP

module sll_m_sim_pic_vm_1d2v_cart_multispecies

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_ascii_io, only: &
       sll_s_ascii_file_create, &
       sll_s_ascii_write_array_1d, &
       sll_s_ascii_write_array_1d_as_row

  use sll_m_filter_base_1d, only: &
       sll_c_filter_base_1d

  use sll_m_fft_filter_1d, only: &
       sll_t_fft_filter_1d

  use sll_m_binomial_filter, only: &
       sll_t_binomial_filter

  use sll_m_cartesian_meshes, only: &
       sll_f_new_cartesian_mesh_1d, &
       sll_t_cartesian_mesh_1d

  use sll_m_collective, only: &
       sll_o_collective_allreduce, &
       sll_s_collective_reduce_real64, &
       sll_f_get_collective_rank, &
       sll_f_get_collective_size, &
       sll_v_world_collective

  use sll_m_constants, only: &
       sll_p_pi, sll_p_twopi

  use sll_m_control_variate, only: &
       sll_t_control_variate, &
       sll_t_control_variates

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

  use sll_m_time_propagator_pic_vm_1d2v_trafo, only: &
       sll_t_time_propagator_pic_vm_1d2v_trafo

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

  use sll_m_initial_distribution, only : &
       sll_c_distribution_params, &
       sll_s_initial_distribution_file_new

  use sll_m_io_utilities, only : &
       sll_s_read_data_real_array, &
       sll_s_concatenate_filename_and_path

  use sll_m_particle_mesh_coupling_base_1d, only: &
       sll_p_galerkin, &
       sll_c_particle_mesh_coupling_1d

  use sll_m_particle_mesh_coupling_spline_1d, only: &
       sll_t_particle_mesh_coupling_spline_1d, &
       sll_s_new_particle_mesh_coupling_spline_1d, &
       sll_s_new_particle_mesh_coupling_spline_1d_ptr

  use sll_m_particle_mesh_coupling_spline_strong_1d, only: &
       sll_t_particle_mesh_coupling_spline_strong_1d, &
       sll_s_new_particle_mesh_coupling_spline_strong_1d

  use sll_m_particle_mesh_coupling_spline_cl_1d, only: &
       sll_t_particle_mesh_coupling_spline_cl_1d, &
       sll_s_new_particle_mesh_coupling_spline_cl_1d, &
       sll_s_new_particle_mesh_coupling_spline_cl_1d_ptr

  use sll_m_particle_mesh_coupling_spline_smooth_1d, only: &
       sll_t_particle_mesh_coupling_spline_smooth_1d, &
       sll_s_new_particle_mesh_coupling_spline_smooth_1d, &
       sll_s_new_particle_mesh_coupling_spline_smooth_1d_ptr

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

  use sll_m_particle_group_1d2v, only: &
       sll_t_particle_group_1d2v

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base, &
       sll_t_particle_array

  use sll_m_particle_sampling, only: &
       sll_t_particle_sampling, &
       sll_s_particle_sampling_randomized_weights

  use sll_m_sim_base, only: &
       sll_c_simulation_base_class

  use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_t_time_mark

  use sll_m_utilities, only : &
       sll_s_int2string

  use sll_mpi, only: &
       mpi_sum, &
       mpi_max

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_splines_pp

  implicit none

  public :: &
       sll_t_sim_pic_vm_1d2v_cart_multispecies

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: sll_p_splitting_hs = 0
  sll_int32, parameter :: sll_p_splitting_boris = 1
  sll_int32, parameter :: sll_p_splitting_disgradE = 2
  sll_int32, parameter :: sll_p_splitting_cef = 4
  sll_int32, parameter :: sll_p_splitting_ecsim = 6
  sll_int32, parameter :: sll_p_splitting_ecsim2o = 7
  sll_int32, parameter :: sll_p_splitting_disgradEC = 10
  sll_int32, parameter :: sll_p_splitting_disgradEC_sub = 13
  sll_int32, parameter :: sll_p_splitting_trafo = 14
  sll_int32, parameter :: sll_p_splitting_momentum = 20

  sll_int32, parameter :: sll_p_onegaussian = 0
  sll_int32, parameter :: sll_p_twogaussian = 1

  sll_int32, parameter :: sll_p_bfield_cos = 0
  sll_int32, parameter :: sll_p_bfield_sin = 1
  sll_int32, parameter :: sll_p_bfield_constant = 2

  sll_int32, parameter :: sll_p_strang_splitting=0
  sll_int32, parameter :: sll_p_splitting_fourth=1
  sll_int32, parameter :: sll_p_lie_splitting=2
  sll_int32, parameter :: sll_p_lie_splitting_back=3
  sll_int32, parameter :: sll_p_splitting_fourth_10steps=4
  sll_int32, parameter :: sll_p_splitting_second_4steps=5


  type, extends(sll_c_simulation_base_class) :: sll_t_sim_pic_vm_1d2v_cart_multispecies

     ! Abstract particle group
     class(sll_t_particle_array), pointer :: particle_group !< Particle group

     !
     sll_real64, pointer :: phi_dofs(:)
     sll_real64, pointer :: efield_dofs(:,:)
     sll_real64, pointer :: efield_dofs_n(:,:)
     sll_real64, allocatable :: bfield_dofs(:)

     sll_real64, allocatable :: x_array(:)
     sll_real64, allocatable :: field_grid(:)

     ! Cartesian mesh
     type(sll_t_cartesian_mesh_1d), pointer    :: mesh
     sll_real64 :: delta_x

     sll_real64 :: time

     ! Maxwell solver 
     ! Abstract 
     class(sll_c_maxwell_1d_base), allocatable :: maxwell_solver
     type(sll_t_maxwell_1d_fem_sm) :: maxwell_norm ! Sparse matrix based Maxwell solver to compute the norm

     ! Abstract kernel smoothers
     class(sll_c_particle_mesh_coupling_1d), allocatable :: kernel_smoother_0     
     class(sll_c_particle_mesh_coupling_1d), allocatable :: kernel_smoother_1


     ! Specific operator splitting
     class(sll_c_time_propagator_base), allocatable :: propagator
     sll_int32 :: splitting_case
     sll_int32 :: splitting_type

     ! Fields on the grid
     sll_real64, allocatable :: fields_grid(:,:)

     ! Control variate
     type(sll_t_control_variates) :: control_variate
     sll_int32  :: no_weights

     ! Physical parameters
     class(sll_c_distribution_params), allocatable :: init_distrib_params_sp1
     class(sll_c_distribution_params), allocatable :: init_distrib_params_sp2
     sll_real64 :: beta
     sll_real64 :: delta(2)
     sll_real64 :: domain(3) ! x_min, x_max, Lx
     sll_real64 :: domain_logical(3)
     type(sll_t_particle_sampling) :: sampler
     sll_real64 :: plasma_betar(3) = 1.0_f64
     sll_real64 :: force_sign 
     logical :: electrostatic 


     ! Simulation parameters
     sll_real64 :: delta_t
     sll_int32  :: n_time_steps
     sll_int32  :: n_particles
     sll_int32  :: n_total_particles
     sll_int32  :: degree_smoother
     sll_int32  :: degree_fem
     sll_int32  :: n_gcells
     sll_int32  :: n_total0
     sll_int32  :: n_total1
     logical    :: boundary = .false.
     sll_int32  :: boundary_fields = 100
     sll_int32  :: boundary_particles = 100


     ! Parameters for MPI
     sll_int32  :: rank
     sll_int32  :: world_size

     ! Case definitions
     sll_int32  :: initial_bfield

     ! Output
     character(len=256)   :: file_prefix
     logical              :: output_fields
     logical              :: output_particles

     ! For ctest
     logical    :: ctest_passed = .false.
     logical    :: make_ctest = .false.
     character(len=256)   :: ctest_ref_file!_rho, ctest_ref_file_thdiag


     logical :: strong_ampere

     !coordinate transformation
     type(sll_t_mapping_3d) :: map
     logical                :: ct = .false.

     !spline_pp
     type(sll_t_spline_pp_1d) :: spline0_pp

     ! Filter
     class(sll_c_filter_base_1d), allocatable :: filter
     type( sll_t_binomial_filter ) :: bfilter

     ! For restart
     logical    :: restart = .false.
     character(len=256) :: restart_file
     sll_int32 :: restart_steps = 0

     ! For random perturbation in initial function
     logical :: randomize_weights(2)
     sll_real64 :: randomize_weights_alpha(2)

   contains
     procedure :: init_from_file => init_pic_vm_1d2v
     procedure :: run => run_pic_vm_1d2v
     procedure :: delete => delete_pic_vm_1d2v

  end type sll_t_sim_pic_vm_1d2v_cart_multispecies


contains
  !------------------------------------------------------------------------------!
  ! Read in the simulation parameters from input file
  subroutine init_pic_vm_1d2v (sim, filename)
    class(sll_t_sim_pic_vm_1d2v_cart_multispecies), intent(inout) :: sim
    character(len=*),                  intent(in)    :: filename

    sll_int32   :: io_stat
    sll_int32   :: input_file, file_id
    sll_int32   :: ierr, j
    type(sll_t_time_mark) :: start_init, end_init
    sll_real64         :: delta_t
    sll_int32          :: n_time_steps
    sll_real64         :: beta
    sll_real64         :: delta(2)= [1._f64, 1._f64]
    character(len=256) :: initial_distrib_sp1
    character(len=256) :: initial_distrib_sp2
    character(len=256) :: initial_bfield
    sll_real64         :: charge(2) = [-1._f64, 1._f64] 
    sll_real64         :: mass(2) = [1._f64, 1._f64]
    sll_real64         :: plasma_beta(3)  = 1._f64
    character(len=256) :: particle_force
    logical            :: electrostatic = .false.
    logical            :: restart = .false.
    character(len=256) :: restart_file
    sll_int32          :: restart_steps = 0
    character(len=256) :: file_prefix
    logical            :: output_fields = .false.
    logical            :: output_particles = .false.
    sll_int32          :: ng_x
    sll_real64         :: x1_min, x1_max
    sll_int32          :: n_particles
    character(len=256) :: sampling_case
    logical            :: delta_perturb = .false.
    sll_real64         :: delta_eps(6)
    character(len=256) :: splitting_case
    sll_int32          :: spline_degree 
    character(len=256) :: splitting_type
    character(len=256) :: boundary_fields = "none"
    character(len=256) :: boundary_particles = "none"
    logical            :: smoothing = .false.
    character(len=256) :: ctest_case
    logical            :: jmean = .false.
    sll_int32          :: degree_fem = -2
    character(len=256) :: filtering
    sll_int32          :: filter_iter = 0
    sll_int32          :: mode = 2
    logical            :: with_control_variate = .false.
    logical            :: strong_ampere = .false.
    logical            :: randomize_weights(2) = .false.
    sll_real64         :: randomize_weights_alpha(2)
    logical            :: eval_grid_points


    namelist /sim_params/         delta_t, n_time_steps, beta, delta, initial_distrib_sp1,initial_distrib_sp2, initial_bfield, charge, mass, plasma_beta, particle_force, electrostatic, restart, restart_file, restart_steps

    namelist /output/             file_prefix, output_fields, output_particles

    namelist /grid_dims/          ng_x, x1_min, x1_max, jmean, degree_fem

    namelist /pic_params/         n_particles, sampling_case, delta_perturb, delta_eps, splitting_case, spline_degree, splitting_type, boundary_fields, boundary_particles, smoothing, filtering, filter_iter, mode, with_control_variate, strong_ampere, randomize_weights, randomize_weights_alpha

    namelist /ctest/              ctest_case

    call sll_s_set_time_mark( start_init )

    ! Read parameters from file
    open(newunit = input_file, file=trim(filename), status='old', IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'init_pic_1d2v() failed to open file ', filename
       STOP
    end if

    read(input_file, sim_params)
    if( restart ) then
       sim%restart = .true.
       sim%restart_file = restart_file
       sim%restart_steps = restart_steps
    end if
    call sll_s_initial_distribution_file_new( [1,2], initial_distrib_sp1, sim%init_distrib_params_sp1 )
    call sll_s_initial_distribution_file_new( [1,2], initial_distrib_sp2, sim%init_distrib_params_sp2 )
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
    sim%delta = delta
    sim%beta = beta
    sim%plasma_betar = 1._f64/plasma_beta
    sim%electrostatic = electrostatic

    select case( particle_force)
    case( "attraction" )
       sim%force_sign = -1._f64
    case( "repulsion" )
       sim%force_sign = 1._f64
    case default
       sim%force_sign = 1._f64
    end select

    select case ( initial_bfield )
    case( "cos")
       sim%initial_bfield = sll_p_bfield_cos
    case( "sin" )
       sim%initial_bfield = sll_p_bfield_sin
    case( "constant" )
       sim%initial_bfield = sll_p_bfield_constant
    case default
       print*, '#initial bfield must be either sin or cos or constant.'
    end select

    ! Output  
    sim%file_prefix = file_prefix
    sim%output_fields = output_fields
    sim%output_particles = output_particles

    sim%n_gcells = ng_x
    sim%mesh => sll_f_new_cartesian_mesh_1d( ng_x, &
         x1_min, x1_max)
    sim%domain = [x1_min, x1_max, x1_max - x1_min ]
    sim%domain_logical = sim%domain
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

    call sim%sampler%init( trim(sampling_case), [1,2], sim%n_particles, sim%rank, delta_perturb, delta_eps )
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
       sim%boundary_fields=sll_p_boundary_clamped
       sim%n_total0= ng_x+spline_degree
       sim%n_total1= ng_x+spline_degree-1
       call sll_s_spline_pp_init_1d( sim%spline0_pp, spline_degree, ng_x, sim%boundary_fields)
    case("clampeddiri")
       sim%boundary = .true.
       sim%boundary_fields=sll_p_boundary_clampeddiri
       sim%n_total0= ng_x+spline_degree
       sim%n_total1= ng_x+spline_degree-1
       call sll_s_spline_pp_init_1d( sim%spline0_pp, spline_degree, ng_x, sim%boundary_fields)
    case("periodic")
       sim%boundary = .true.
       sim%boundary_fields=sll_p_boundary_periodic
       sim%n_total0= ng_x+spline_degree
       sim%n_total1= ng_x+spline_degree-1
       call sll_s_spline_pp_init_1d( sim%spline0_pp, spline_degree, ng_x, sim%boundary_fields)
    case default
       sim%boundary = .false.
       sim%n_total0= ng_x
       sim%n_total1= ng_x
    end select

    select case(boundary_particles)
    case( "reflection" )
       sim%boundary_particles = sll_p_boundary_particles_reflection
    case( "absorption" )
       sim%boundary_particles = sll_p_boundary_particles_absorption
    case( "periodic" )
       sim%boundary_particles = sll_p_boundary_particles_periodic
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
    case("splitting_trafo")
       sim%splitting_case = sll_p_splitting_trafo
       sim%ct=.true.
       sim%domain_logical = [0._f64, 1._f64, x1_max - x1_min ]
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
       sim%ctest_ref_file = "reffile_pic_vm_1d2v_cart_multispecies.dat"
!!$    case("ecsim")
!!$       sim%make_ctest = .true.
!!$       sim%ctest_ref_file = "reffile_ecsim_pic_vm_1d2v_cart_rho.dat"
    end select


    ! Initialize the particles   (mass set to 1.0 and charge set to -1.0)
    allocate( sim%particle_group )
    sim%particle_group%n_species = 2
    allocate( sll_t_particle_group_1d2v :: sim%particle_group%group(sim%particle_group%n_species) )
    ! electrons
    select type ( qp => sim%particle_group%group(1) )
    type is (  sll_t_particle_group_1d2v )
       ! Note: This call produces a segmentation fault with the INTEL 17 compiler
       ! Therefore we manually initialize here
       ! TODO: Fix the problem with the init function
       call qp%init(sim%n_particles, &
            sim%n_total_particles, charge(1), mass(1), sim%no_weights)
       call qp%set_common_weight(sim%delta(1))

       !qp%n_particles = sim%n_particles
       !qp%n_total_particles = sim%n_total_particles
       !SLL_ALLOCATE( qp%particle_array(3+sim%no_weights,sim%n_particles), ierr)
       !allocate( qp%species, stat=ierr)
       !SLL_ASSERT(ierr==0)
       !call qp%species%init( charge(1), mass(1) )
       !qp%n_weights = sim%no_weights
    end select
    ! ions
    select type ( qp => sim%particle_group%group(2) )
    type is (  sll_t_particle_group_1d2v )
       ! Note: This call produces a segmentation fault with the INTEL 17 compiler
       ! Therefore we manually initialize here
       ! TODO: Fix the problem with the init function
       call qp%init(sim%n_particles, &
            sim%n_total_particles, charge(2), mass(2), sim%no_weights)
       call qp%set_common_weight(sim%delta(2))
       !qp%n_particles = sim%n_particles
       !qp%n_total_particles = sim%n_total_particles
       !SLL_ALLOCATE( qp%particle_array(3+sim%no_weights,sim%n_particles), ierr)
       !allocate( qp%species, stat=ierr)
       !SLL_ASSERT(ierr==0)
       !call qp%species%init( charge(2), mass(2) )
       !qp%n_weights = sim%no_weights
       !call qp%init(sim%n_particles, &
       !     sim%n_total_particles, 1.0_f64, mass_ion, sim%no_weights)
    end select

    if (sim%rank == 0 ) then
       open(newunit=file_id, file=trim(filename)//'_used.dat')
       close(file_id)
    end if

    ! Initialize control variate
    allocate(sim%control_variate%cv(2) )
    call sim%control_variate%cv(1)%init(control_variate_equi, &
         distribution_params=sim%init_distrib_params_sp1)
    call sim%control_variate%cv(2)%init(control_variate_equi, &
         distribution_params=sim%init_distrib_params_sp2)


    if (sim%rank == 0 ) then
       open(newunit=file_id, file=trim(filename)//'_used.dat')
       close(file_id)
    end if

    ! Initialize kernel smoother
    if( sim%boundary ) then
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
                  sim%degree_smoother, sim%boundary_fields , trim(filename) )
          end select
       end if
       call sim%maxwell_norm%init_from_file( sim%domain(1:2), sim%n_gcells, sim%degree_smoother, trim(filename) )

       call sll_s_new_particle_mesh_coupling_spline_cl_1d(sim%kernel_smoother_1, &
            sim%domain_logical(1:2), sim%n_gcells, &
            sim%n_particles, sim%degree_smoother-1, sll_p_galerkin, sim%boundary_fields) 
       call sll_s_new_particle_mesh_coupling_spline_cl_1d(sim%kernel_smoother_0, &
            sim%domain_logical(1:2), sim%n_gcells, &
            sim%n_particles, sim%degree_smoother, sll_p_galerkin, sim%boundary_fields)
    else
       ! Initialize the field solver
       if (sim%ct) then
          allocate( sll_t_maxwell_1d_trafo :: sim%maxwell_solver )
          select type ( q=>sim%maxwell_solver )
          type is ( sll_t_maxwell_1d_trafo )
             call q%init_from_file( sim%domain(1:2), sim%n_gcells, sim%degree_smoother, sim%map,  trim(filename) )
          end select
       else
          if ( sim%degree_fem > -1 ) then
             allocate( sll_t_maxwell_1d_fem :: sim%maxwell_solver )
             select type ( q=>sim%maxwell_solver )
             type is ( sll_t_maxwell_1d_fem )
                call q%init( sim%domain(1:2), sim%n_gcells, &
                     sim%degree_fem, delta_t*0.5_f64, strong_ampere = sim%strong_ampere )
             end select
!!$          allocate( sll_t_maxwell_1d_fem_sm :: sim%maxwell_solver )
!!$          select type ( q=>sim%maxwell_solver )
!!$          type is ( sll_t_maxwell_1d_fem_sm )
!!$             call q%init_from_file( sim%domain(1:2), sim%n_gcells, &
!!$                  sim%degree_fem, trim(filename) )
!!$          end select
          else
             allocate( sll_t_maxwell_1d_ps :: sim%maxwell_solver )
             select type ( q=>sim%maxwell_solver )
             type is ( sll_t_maxwell_1d_ps )
                call q%init( sim%domain(1:2), sim%n_gcells )
             end select
          end if
       end if
       if ( sim%degree_fem > -1) then
          call sim%maxwell_norm%init_from_file( sim%domain(1:2), sim%n_gcells, sim%degree_fem, trim(filename) )
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
       end if
    end if

    ! Initialize the arrays for the spline coefficients of the fields
    SLL_ALLOCATE(sim%phi_dofs(sim%n_gcells), ierr)
    SLL_ALLOCATE(sim%efield_dofs(sim%n_gcells,2), ierr)
    SLL_ALLOCATE(sim%bfield_dofs(sim%n_gcells), ierr)
    sim%efield_dofs = 0._f64
    sim%bfield_dofs = 0._f64

    ! Initialize the time-splitting propagator
    if (sim%splitting_case == sll_p_splitting_hs) then
       if (sim%no_weights == 1) then
          call sll_s_new_time_propagator_pic_vm_1d2v_hs(&
               sim%propagator, sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), sim%filter, force_sign=sim%force_sign, jmean=jmean, betar = sim%plasma_betar(1:2), electrostatic=sim%electrostatic)
       else
          call sll_s_new_time_propagator_pic_vm_1d2v_hs(&
               sim%propagator, sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), sim%filter, force_sign=sim%force_sign, jmean=jmean, control_variate=sim%control_variate, i_weight=sim%no_weights, betar = sim%plasma_betar(1:2), electrostatic=sim%electrostatic)
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
                  sim%domain(1), sim%domain(3), sim%filter, trim(filename), boundary_particles=sim%boundary_particles, force_sign=sim%force_sign, betar=sim%plasma_betar(1:2), electrostatic=sim%electrostatic, jmean=jmean  )
          else
             call qpdisgradE%init_from_file( sim%maxwell_solver, &
                  sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
                  sim%phi_dofs, sim%efield_dofs, sim%bfield_dofs, &
                  sim%domain(1), sim%domain(3), sim%filter, trim(filename), sim%boundary_particles, sim%force_sign, sim%control_variate, sim%no_weights, sim%plasma_betar(1:2), sim%electrostatic, jmean  )
          end if
          sim%efield_dofs_n => qpdisgradE%efield_dofs
       end select
    elseif  (sim%splitting_case == sll_p_splitting_cef) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_cef :: sim%propagator )
       select type( qpcef=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_cef )
          call qpcef%init( sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, sim%domain(1), sim%domain(3) )
          sim%efield_dofs_n => sim%efield_dofs
       end select
    elseif( sim%splitting_case == sll_p_splitting_ecsim) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_ecsim :: sim%propagator )
       select type( qpecsim=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_ecsim )
          call qpecsim%init_from_file(sim%kernel_smoother_0, sim%kernel_smoother_1, &
               sim%particle_group, sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), trim(filename) )
          sim%efield_dofs_n => qpecsim%efield_dofs
       end select
    elseif( sim%splitting_case == sll_p_splitting_ecsim2o) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_ecsim2o :: sim%propagator )
       select type( qpecsim2o=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_ecsim2o )
          call qpecsim2o%init(sim%kernel_smoother_0, sim%kernel_smoother_1, &
               sim%particle_group, sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3))
          sim%efield_dofs_n => qpecsim2o%efield_dofs
       end select
    elseif  (sim%splitting_case == sll_p_splitting_disgradEC) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_disgradEC :: sim%propagator )
       select type( qpdisgradEC=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_disgradEC )
          call qpdisgradEC%init( sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), sim%filter, trim(filename) )
          sim%efield_dofs_n => qpdisgradEC%helper%efield_dofs
       end select
    elseif  (sim%splitting_case == sll_p_splitting_disgradEC_sub) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub :: sim%propagator )
       select type( qpdgs=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_disgradEC_sub )
          call qpdgs%init( sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), sim%bfilter, trim(filename) )
          sim%efield_dofs_n => qpdgs%helper%efield_dofs
       end select
    elseif( sim%splitting_case == sll_p_splitting_trafo) then
       allocate( sll_t_time_propagator_pic_vm_1d2v_trafo :: sim%propagator )
       select type( qptrafo=>sim%propagator )
       type is ( sll_t_time_propagator_pic_vm_1d2v_trafo )
          call qptrafo%init_from_file( sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3), sim%map, trim(filename), sim%boundary_particles, force_sign=sim%force_sign, electrostatic=sim%electrostatic, jmean=jmean  )!,betar=sim%plasma_betar(1:2))
          sim%efield_dofs_n => qptrafo%efield_dofs
       end select
    end if

    ! Allocate the vector holding the values of the fields at the grid points
    SLL_ALLOCATE(sim%fields_grid(sim%n_gcells,3), ierr)

    allocate(sim%x_array(sim%n_gcells))
    allocate(sim%field_grid(sim%n_gcells))

    sim%x_array(1) = sim%domain(1)
    do j=2,sim%n_gcells
       sim%x_array(j) = sim%x_array(j-1) + (sim%domain(3)/real(sim%n_gcells, f64))
    end do

    sim%randomize_weights = randomize_weights
    sim%randomize_weights_alpha = randomize_weights_alpha


    call sll_s_set_time_mark( end_init )
    if (sim%rank == 0 ) then
       sim%time = sll_f_time_elapsed_between( start_init, end_init)
       write(*, "(A, F10.3)") "Init run time [s] = ", sim%time

       open(newunit=file_id, file=trim(filename)//'_used.dat', position = 'append', status='old', action='write', iostat=ierr)
       write(file_id, *) 'delta t:', sim%delta_t
       write(file_id, *) 'n_time_steps:', sim%n_time_steps
       write(file_id, *) 'beta:', sim%beta
       write(file_id, *) 'delta:', delta
       write(file_id, *) 'charge:', charge
       write(file_id, *) 'mass:', mass
       write(file_id, *) 'plasma betar:', sim%plasma_betar
       write(file_id, *) 'force sign simulation:', sim%force_sign
       write(file_id, *) 'electrostatic simulation:', sim%electrostatic
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

  subroutine run_pic_vm_1d2v (sim)
    class(sll_t_sim_pic_vm_1d2v_cart_multispecies), intent(inout) :: sim

    ! Local variables
    sll_int32 :: j, ierr, i_part, i, i_steps
    sll_real64, allocatable :: rho(:), rho_local(:), efield_poisson(:)
    sll_int32 :: th_diag_id, dfield_id, efield_id, bfield_id, rho_id
    character(len=4) :: crank
    character(len=4) :: step
    character(len=256) :: diag_file_name
    sll_real64 :: wi(1)
    sll_real64 :: xi(3)
    sll_int32 :: file_id
    type(sll_t_time_mark) :: start_loop, end_loop
    sll_int32 :: degreeb

    print*, 'Resolution dt, dx', sim%delta_t, sim%delta_x


    ! Initialize file for diagnostics
    if (sim%rank == 0) then
       diag_file_name = trim(sim%file_prefix)//"_thdiag.dat"
       call sll_s_ascii_file_create(trim(diag_file_name), th_diag_id, ierr)
       if ( sim%output_fields ) then
          call sll_s_ascii_file_create(trim(sim%file_prefix)//'_dfield.dat', dfield_id, ierr)
          call sll_s_ascii_file_create(trim(sim%file_prefix)//'_efield.dat', efield_id, ierr)
          call sll_s_ascii_file_create(trim(sim%file_prefix)//'_bfield.dat', bfield_id, ierr)
          call sll_s_ascii_file_create(trim(sim%file_prefix)//'_rho.dat', rho_id, ierr)
       end if
    end if


    call sll_s_int2string( sim%rank, crank )
    if ( sim%restart ) then
       i_steps = real(sim%restart_steps, f64)*sim%delta_t
       call sll_s_int2string( i_steps, step )
       call sim%particle_group%group(1)%read(trim(sim%restart_file)//step//'_particles_end_'//crank//'_sp1.dat')
       call sim%particle_group%group(2)%read(trim(sim%restart_file)//step//'_particles_end_'//crank//'_sp2.dat')
    else
       if(sim%ct) then
          call sim%sampler%sample( sim%particle_group%group(1), sim%init_distrib_params_sp1, sim%domain(1:1), sim%domain(3:3), sim%map )
          call sim%sampler%reset_seed_jump ( sim%n_total_particles )

          call sim%sampler%sample( sim%particle_group%group(2), sim%init_distrib_params_sp2, sim%domain(1:1), sim%domain(3:3), sim%map )
       else
          if (sim%no_weights == 1) then
             call sim%sampler%sample( sim%particle_group%group(1), sim%init_distrib_params_sp1, sim%domain(1:1), sim%domain(3:3) )
             call sim%sampler%reset_seed_jump ( sim%n_total_particles )
             call sim%sampler%sample( sim%particle_group%group(2), sim%init_distrib_params_sp2, sim%domain(1:1), sim%domain(3:3) )
          else
             call sim%sampler%sample_cv( sim%particle_group%group(1), sim%init_distrib_params_sp1, sim%domain(1:1), sim%domain(3:3), sim%control_variate%cv(1) )
             call sim%sampler%reset_seed_jump ( sim%n_total_particles )
             call sim%sampler%sample_cv( sim%particle_group%group(2), sim%init_distrib_params_sp2, sim%domain(1:1), sim%domain(3:3), sim%control_variate%cv(2) )
          end if
       end if
    end if

    if ( sim%restart .eqv. .false. ) then
       if ( sim%randomize_weights(1) .eqv. .true. ) then
          call sll_s_particle_sampling_randomized_weights(sim%particle_group%group(1), sim%randomize_weights_alpha(1))
       end if
       if ( sim%randomize_weights(2) .eqv. .true. ) then
          call sll_s_particle_sampling_randomized_weights(sim%particle_group%group(2), sim%randomize_weights_alpha(2))
       end if
    end if


    ! Print particle array to file
    if ( sim%output_particles ) then
       call sll_s_int2string( sim%rank, crank )
       call sim%particle_group%group(1)%print(trim(sim%file_prefix)//'_particles_start_1_'//crank//'.dat')
       call sim%particle_group%group(2)%print(trim(sim%file_prefix)//'_particles_start_2_'//crank//'.dat')
    end if

    ! Set the initial fields
    SLL_ALLOCATE(rho_local(sim%n_total0), ierr)
    SLL_ALLOCATE(rho(sim%n_total0), ierr)
    SLL_ALLOCATE(efield_poisson(sim%n_total0), ierr)


    if ( sim%restart ) then
       open(newunit=file_id, file=trim(sim%restart_file)//step//'_efield.dat', status='old', action='read', iostat=ierr)
       if (ierr /= 0 ) then
          SLL_ERROR("run", "Restart file for efield does not exist: "//trim(sim%restart_file)//step//'_efield.dat')
       end if
       read(file_id, *) sim%efield_dofs
       close(file_id)
       open(newunit=file_id, file=trim(sim%restart_file)//step//'_bfield.dat', status='old', action='read', iostat=ierr)
       if (ierr /= 0 ) then
          SLL_ERROR("run", "Restart file for efield does not exist: "//trim(sim%restart_file)//step//'_bfield.dat')
       end if
       read(file_id, *) sim%bfield_dofs
       close(file_id)
    else
       ! Efield 1 by Poisson
       call solve_poisson( sim, rho_local, rho )
       ! Efield 2 to zero
       sim%efield_dofs(:,2) = 0.0_f64

       if ( sim%strong_ampere .eqv. .true. ) then
          degreeb = sim%degree_fem
       else
          degreeb = sim%degree_fem-1
       end if

       ! Bfield = beta*cos(kx): Use b = M{-1}(N_i,beta*cos(kx))    
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
       end select
    end if

    call sim%propagator%reinit_fields()

    ! In case we use the Boris propagator, we need to initialize the staggering used in the scheme.
    select type( qp=>sim%propagator )
    type is ( sll_t_time_propagator_pic_vm_1d2v_boris)
       call qp%staggering( sim%delta_t )
    end select

    ! End field initialization

    !call solve_poisson( sim%particle_group, sim%kernel_smoother_0, sim%maxwell_solver, sim%delta_x, rho_local, rho, efield_poisson )
    ! Diagnostics
    if(sim%restart .eqv. .false.)then
       call sll_s_time_history_diagnostics_pic_vm_1d2v( &
            sim, 0.0_f64, th_diag_id, rho, efield_poisson)
    end if
    if (sim%rank == 0 ) then
       call sll_s_set_time_mark(start_loop )
       if ( sim%output_fields ) then  
          call sll_s_ascii_write_array_1d_as_row( dfield_id, &
               sim%efield_dofs(:,1),  &
               sim%n_gcells )
          call sll_s_ascii_write_array_1d_as_row( efield_id, &
               sim%efield_dofs(:,2),  &
               sim%n_gcells )
          call sll_s_ascii_write_array_1d_as_row( bfield_id, &
               sim%bfield_dofs,  &
               sim%n_gcells )
          call sll_s_ascii_write_array_1d_as_row( rho_id, &
               rho,  &
               sim%n_gcells )
       end if
    end if

    ! Time loop
    select case (sim%splitting_type)
    case(sll_p_strang_splitting)
       do j=1+sim%restart_steps, sim%n_time_steps+sim%restart_steps
          !print*, 'TIME STEP', j
          call sim%propagator%strang_splitting(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_1d2v( &
               sim,  sim%delta_t*real(j,f64), th_diag_id, rho, efield_poisson)

          if( sim%output_particles ) then
             if( modulo(j, 100) == 0 ) then
                call sll_s_int2string( j, crank )
                call sim%particle_group%group(1)%print(trim(sim%file_prefix)//'_particles_1_'//crank//'.dat')
                call sim%particle_group%group(2)%print(trim(sim%file_prefix)//'_particles_2_'//crank//'.dat')
             end if
          end if

!!$          if (sim%rank == 0 ) then     
!!$             if ( sim%output_fields ) then     
!!$                call sll_s_ascii_write_array_1d_as_row( dfield_id, &
!!$                     sim%efield_dofs(:,1),  &
!!$                     sim%n_gcells )
!!$                call sll_s_ascii_write_array_1d_as_row( efield_id, &
!!$                     sim%efield_dofs(:,2),  &
!!$                     sim%n_gcells )
!!$                call sll_s_ascii_write_array_1d_as_row( bfield_id, &
!!$                     sim%bfield_dofs,  &
!!$                     sim%n_gcells )
!!$             end if
!!$          end if
       end do
    case(sll_p_splitting_fourth)
       do j=1+sim%restart_steps, sim%n_time_steps+sim%restart_steps
          !print*, 'TIME STEP', j
          call sim%propagator%splitting_fourth(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_1d2v( &
               sim,  sim%delta_t*real(j,f64), th_diag_id, rho, efield_poisson)

!!$          if (sim%rank == 0 ) then    
!!$             if ( sim%output_fields ) then      
!!$                call sll_s_ascii_write_array_1d_as_row( dfield_id, &
!!$                     sim%efield_dofs(:,1),  &
!!$                     sim%n_gcells )
!!$                call sll_s_ascii_write_array_1d_as_row( efield_id, &
!!$                     sim%efield_dofs(:,2),  &
!!$                     sim%n_gcells )
!!$                call sll_s_ascii_write_array_1d_as_row( bfield_id, &
!!$                     sim%bfield_dofs,  &
!!$                     sim%n_gcells )
!!$             end if
!!$          end if
       end do
    case(sll_p_lie_splitting)
       do j=1+sim%restart_steps, sim%n_time_steps+sim%restart_steps
          call sim%propagator%lie_splitting(sim%delta_t,1)
          !print*, 'TIME STEP', j
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_1d2v( &
               sim,  sim%delta_t*real(j,f64), th_diag_id, rho, efield_poisson)

!!$          if (sim%rank == 0 ) then     
!!$             if ( sim%output_fields ) then     
!!$                call sll_s_ascii_write_array_1d_as_row( dfield_id, &
!!$                     sim%efield_dofs(:,1),  &
!!$                     sim%n_gcells )
!!$                call sll_s_ascii_write_array_1d_as_row( efield_id, &
!!$                     sim%efield_dofs(:,2),  &
!!$                     sim%n_gcells )
!!$                call sll_s_ascii_write_array_1d_as_row( bfield_id, &
!!$                     sim%bfield_dofs,  &
!!$                     sim%n_gcells )
!!$             end if
!!$          end if
       end do
    case(sll_p_lie_splitting_back)
       do j=1+sim%restart_steps, sim%n_time_steps+sim%restart_steps
          !print*, 'TIME STEP', j
          call sim%propagator%lie_splitting_back(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_1d2v( &
               sim,  sim%delta_t*real(j,f64), th_diag_id, rho, efield_poisson)

!!$          if (sim%rank == 0 ) then      
!!$             if ( sim%output_fields ) then    
!!$                call sll_s_ascii_write_array_1d_as_row( dfield_id, &
!!$                     sim%efield_dofs(:,1),  &
!!$                     sim%n_gcells )
!!$                call sll_s_ascii_write_array_1d_as_row( efield_id, &
!!$                     sim%efield_dofs(:,2),  &
!!$                     sim%n_gcells )
!!$                call sll_s_ascii_write_array_1d_as_row( bfield_id, &
!!$                     sim%bfield_dofs,  &
!!$                     sim%n_gcells )
!!$             end if
!!$          end if
       end do
    case(sll_p_splitting_fourth_10steps)
       do j=1+sim%restart_steps, sim%n_time_steps+sim%restart_steps
          !print*, 'TIME STEP', j
          call sim%propagator%splitting_fourth_10steps(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_1d2v( &
               sim,  sim%delta_t*real(j,f64), th_diag_id, rho, efield_poisson)

!!$          if (sim%rank == 0 ) then    
!!$             if ( sim%output_fields ) then      
!!$                call sll_s_ascii_write_array_1d_as_row( dfield_id, &
!!$                     sim%efield_dofs(:,1),  &
!!$                     sim%n_gcells )
!!$                call sll_s_ascii_write_array_1d_as_row( efield_id, &
!!$                     sim%efield_dofs(:,2),  &
!!$                     sim%n_gcells )
!!$                call sll_s_ascii_write_array_1d_as_row( bfield_id, &
!!$                     sim%bfield_dofs,  &
!!$                     sim%n_gcells )
!!$             end if
!!$          end if
       end do
    case(sll_p_splitting_second_4steps)
       do j=1+sim%restart_steps, sim%n_time_steps+sim%restart_steps
          !print*, 'TIME STEP', j
          call sim%propagator%splitting_second_4steps(sim%delta_t,1)
          ! Diagnostics
          call sll_s_time_history_diagnostics_pic_vm_1d2v( &
               sim,  sim%delta_t*real(j,f64), th_diag_id, rho, efield_poisson)

!!$          if (sim%rank == 0 ) then    
!!$             if ( sim%output_fields ) then      
!!$                call sll_s_ascii_write_array_1d_as_row( dfield_id, &
!!$                     sim%efield_dofs(:,1),  &
!!$                     sim%n_gcells )
!!$                call sll_s_ascii_write_array_1d_as_row( efield_id, &
!!$                     sim%efield_dofs(:,2),  &
!!$                     sim%n_gcells )
!!$                call sll_s_ascii_write_array_1d_as_row( bfield_id, &
!!$                     sim%bfield_dofs,  &
!!$                     sim%n_gcells )
!!$             end if
!!$          end if
       end do
    case default
       print*, 'this splitting type is not implemented'
    end select

    if (sim%rank == 0 ) then
       call sll_s_set_time_mark( end_loop )
       write(*, "(A, F10.3)") "Main loop run time [s] = ", sll_f_time_elapsed_between( start_loop, end_loop)
       close(th_diag_id)
       if(sim%output_fields)then
          i_steps = real(sim%n_time_steps+sim%restart_steps, f64)*sim%delta_t
          call sll_s_int2string( i_steps, step)
          open(newunit=file_id, file=trim(sim%file_prefix)//step//'efield.dat')
          write(file_id, *) sim%efield_dofs
          close(file_id)
          open(newunit=file_id, file=trim(sim%file_prefix)//step//'bfield.dat')
          write(file_id, *) sim%bfield_dofs
          close(file_id)
       end if
    end if

    ! Print particle array to file
    if ( sim%output_particles ) then
       call sll_s_int2string( sim%rank, crank )
       i_steps = real(sim%n_time_steps+sim%restart_steps, f64)*sim%delta_t
       call sll_s_int2string( i_steps, step)
       call sim%particle_group%group(1)%print(trim(sim%file_prefix)//step//'_particles_end_1_'//crank//'.dat')
       call sim%particle_group%group(2)%print(trim(sim%file_prefix)//step//'_particles_end_2_'//crank//'_.dat')
    end if

!!! Part for ctest
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
            sim%n_total0, MPI_SUM, rho)
       !write(95,*) rho
       if (sim%rank == 0) then
          call ctest( rho, rho_local, trim(sim%ctest_ref_file), sim%ctest_passed )
          !call sll_s_check_diagnostics(trim(sim%ctest_ref_file_thdiag),'ctest.dat', 1E-13_f64, sim%ctest_passed)
       end if
    end if
!!! Part for ctest end


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



  end subroutine run_pic_vm_1d2v

  !------------------------------------------------------------------------------!
  ! local subroutine to handle ctest
  subroutine ctest(rho_simulated, rho_ref, ctest_ref_file, passed)
    sll_real64, intent(in   ) :: rho_simulated(:)
    sll_real64, intent(inout) :: rho_ref(:)
    character(*), intent(in   ) :: ctest_ref_file
    logical,    intent(  out) :: passed     

    ! For testing
    character(len=256) :: reffile
    sll_real64 :: error

    call sll_s_concatenate_filename_and_path( trim(ctest_ref_file), __FILE__,&
         reffile)
    call sll_s_read_data_real_array( reffile, rho_ref)

    rho_ref = rho_ref -  rho_simulated
    error = maxval(rho_ref)
    print*, 'Maximum error in rho is', error, '.'
    if (abs(error)> 1E-14) then
       passed = .FALSE.
    else
       passed = .TRUE.
    end if

  end subroutine ctest


  !------------------------------------------------------------------------------!

  subroutine delete_pic_vm_1d2v (sim)
    class(sll_t_sim_pic_vm_1d2v_cart_multispecies), intent(inout) :: sim
    SLL_ASSERT(storage_size(sim)>0)

    call sim%propagator%free()
    deallocate(sim%propagator)
    call sim%particle_group%group(1)%free()
    call sim%particle_group%group(2)%free()
    deallocate (sim%particle_group)
    call sim%mesh%delete()
    deallocate(sim%mesh)
    call sim%maxwell_solver%free()
    deallocate(sim%maxwell_solver)
    call sim%kernel_smoother_0%free()
    deallocate(sim%kernel_smoother_0)
    call sim%kernel_smoother_1%free()
    deallocate(sim%kernel_smoother_1)
    deallocate(sim%fields_grid)
    call sim%control_variate%cv(1)%free()
    call sim%control_variate%cv(2)%free()

    if ( sim%restart .eqv. .false. ) then
       call sim%init_distrib_params_sp1%free()
       deallocate(sim%init_distrib_params_sp1)
       call sim%init_distrib_params_sp2%free()
       deallocate(sim%init_distrib_params_sp2)
    end if
    call sim%sampler%free()

  end subroutine delete_pic_vm_1d2v


  !> As a control variate, we use the equilibrium (v part of the initial distribution)
  function control_variate_equi( this, xi, vi, time) result(sll_f_control_variate)
    class(sll_t_control_variate) :: this
    sll_real64, optional,  intent( in ) :: xi(:) !< particle position
    sll_real64, optional,  intent( in ) :: vi(:) !< particle velocity
    sll_real64, optional,  intent( in ) :: time  !< current time
    sll_real64               :: sll_f_control_variate


    sll_f_control_variate = &
         this%control_variate_distribution_params%eval_v_density( vi(1:2) )


  end function control_variate_equi


  !------------------------------------------------------------------------------!
  !Diagnostic functions and other helper functions
  !> Diagnostics for PIC Vlasov-Maxwell 1d2v 
  !> @todo (should be part of the library)
  subroutine sll_s_time_history_diagnostics_pic_vm_1d2v(&
       sim,&
       time, &
       file_id, &
       rho, scratch)
    class(sll_t_sim_pic_vm_1d2v_cart_multispecies), intent(inout) :: sim
    sll_real64,                        intent(in   ) :: time
    sll_int32,                         intent(in   ) :: file_id
    sll_real64,                        intent(  out) :: rho(:)
    sll_real64,                        intent(  out) :: scratch(:)

    ! local variables
    sll_real64 :: diagnostics_local(6)
    sll_real64 :: diagnostics(6)
    sll_real64 :: potential_energy(3)
    sll_int32  :: i_part, i_sp
    sll_int32  :: degree
    sll_real64 :: vi(3),  xi(3)
    sll_real64 :: wi(1)
    sll_real64 :: transfer(1), vvb(1), poynting
    sll_real64 :: efield(2), bfield, error

    degree = sim%degree_fem

    diagnostics_local = 0.0_f64
    do i_sp = 1, sim%particle_group%n_species
       do i_part=1, sim%particle_group%group(i_sp)%n_particles
          vi = sim%particle_group%group(i_sp)%get_v(i_part)
          xi = sim%particle_group%group(i_sp)%get_x(i_part)
          wi = sim%particle_group%group(i_sp)%get_mass(i_part)

          ! Kinetic energy
          diagnostics_local(2*i_sp-1) = diagnostics_local(2*i_sp-1) + &
               0.5_f64*(vi(1)**2)*wi(1)
          diagnostics_local(2*i_sp) = diagnostics_local(2*i_sp) + &
               0.5_f64*(vi(2)**2)*wi(1)
          ! Momentum 1
          diagnostics_local(5) = diagnostics_local(5) + &
               vi(1)*wi(1)
          ! Momentum 2
          diagnostics_local(6) = diagnostics_local(6) + &
               vi(2)*wi(1)
       end do
    end do
    diagnostics = 0.0_f64
    call sll_s_collective_reduce_real64(sll_v_world_collective, diagnostics_local, 6,&
         MPI_SUM, 0, diagnostics)
    ! Add ExB part
    if( sim%boundary )then
       !TODO?
       diagnostics(5) = diagnostics(5) + sim%maxwell_solver%inner_product( sim%efield_dofs(:,2), sim%bfield_dofs, degree, degree-1 )

       diagnostics(5) = diagnostics(6) - sim%maxwell_solver%inner_product( sim%efield_dofs(1:sim%n_total1,1), sim%bfield_dofs, degree-1 )
    else
       if ( sim%strong_ampere .eqv. .false. ) then
          diagnostics(5) = diagnostics(5) + sim%maxwell_norm%inner_product( sim%efield_dofs(:,2), sim%bfield_dofs, degree, degree-1 )
          diagnostics(6) = diagnostics(6) - sim%maxwell_solver%inner_product( sim%efield_dofs(:,1), sim%bfield_dofs, degree-1 )
       else
          if ( degree == -1) then
             diagnostics(5) = diagnostics(5) + sim%maxwell_solver%inner_product( sim%efield_dofs(:,2), sim%bfield_dofs, degree-1, degree )
             diagnostics(6) = diagnostics(6) - sim%maxwell_solver%inner_product( sim%efield_dofs(:,1), sim%bfield_dofs, degree )
          else
             !scratch(2:sim%n_gcells) = 0.5_f64 * (sim%efield_dofs(1:sim%n_gcells-1,2)+sim%efield_dofs(2:sim%n_gcells,2) )
             !scratch(1) = 0.5_f64 * (sim%efield_dofs(1,2)+sim%efield_dofs(sim%n_gcells,2) )
             scratch(1:sim%n_gcells-1) = 0.5_f64 * (sim%efield_dofs(1:sim%n_gcells-1,2)+sim%efield_dofs(2:sim%n_gcells,2) )
             scratch(sim%n_gcells) = 0.5_f64 * (sim%efield_dofs(1,2)+sim%efield_dofs(sim%n_gcells,2) )

             diagnostics(5) = diagnostics(5) + sim%maxwell_solver%inner_product( scratch, sim%bfield_dofs, degree )
             !scratch(2:sim%n_gcells) = 0.5_f64 * (sim%efield_dofs(1:sim%n_gcells-1,1)+sim%efield_dofs(2:sim%n_gcells,1) )
             !scratch(1) = 0.5_f64 * (sim%efield_dofs(1,1)+sim%efield_dofs(sim%n_gcells,1) )

             diagnostics(6) = diagnostics(6) - sim%maxwell_solver%inner_product( sim%efield_dofs(:,1), sim%bfield_dofs, degree )
          end if
       end if
    end if

    ! Check error in Gauss law
    call check_gauss_law( sim, rho, scratch, error )

    if (sim%rank == 0) then
       if( sim%boundary )then
          potential_energy(1) = sim%maxwell_solver%L2norm_squared( sim%efield_dofs(1:sim%n_total1,1), degree-1 )/sim%plasma_betar(2)
          potential_energy(2) = sim%maxwell_solver%L2norm_squared( sim%efield_dofs(:,2), degree )/sim%plasma_betar(2)
          potential_energy(3) = sim%maxwell_solver%L2norm_squared&
               ( sim%bfield_dofs, degree-1 )*sim%plasma_betar(3)
       else
          if ( sim%strong_ampere .eqv. .false. ) then
             potential_energy(1) = sim%maxwell_solver%inner_product&
                  ( sim%efield_dofs(:,1), sim%efield_dofs_n(:,1), degree-1 )/sim%plasma_betar(2)
             potential_energy(2) = sim%maxwell_solver%inner_product&
                  ( sim%efield_dofs(:,2), sim%efield_dofs_n(:,2), degree )/sim%plasma_betar(2)
             potential_energy(3) = sim%maxwell_solver%L2norm_squared&
                  ( sim%bfield_dofs, degree-1 )*sim%plasma_betar(3)
          else
             potential_energy(1) = sim%maxwell_solver%inner_product&
                  ( sim%efield_dofs(:,1), sim%efield_dofs_n(:,1), degree )/sim%plasma_betar(2)
             potential_energy(2) = sim%maxwell_solver%inner_product&
                  ( sim%efield_dofs(:,2), sim%efield_dofs_n(:,2), degree-1 )/sim%plasma_betar(2)
             potential_energy(3) = sim%maxwell_solver%L2norm_squared&
                  ( sim%bfield_dofs, degree )*sim%plasma_betar(3)
          end if
       end if
       write(file_id,'(f12.5,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16)' ) &
            time,  potential_energy, diagnostics(1:4), &
            sum(diagnostics(1:4)) + sim%force_sign * 0.5_f64*sum(potential_energy), diagnostics(5:6), &
            error
    end if

  end subroutine sll_s_time_history_diagnostics_pic_vm_1d2v


  subroutine sll_s_diagnostics_fields( field_dofs,  field_grid, xi, n_dofs, kernel_smoother, file_id )
    sll_real64,                       intent(in) :: field_dofs(:)
    sll_real64,                       intent(inout) :: field_grid(:)
    sll_real64,                       intent(in) :: xi(:)
    sll_int32, intent(in) :: n_dofs
    class(sll_c_particle_mesh_coupling_1d),     intent(inout) :: kernel_smoother
    sll_int32, intent(in) :: file_id


    sll_int32 :: j

    do j=1,n_dofs
       call kernel_smoother%evaluate( [xi(j)], field_dofs, field_grid(j))
    end do

    write(file_id, *) field_grid


  end subroutine sll_s_diagnostics_fields
  !> compute v(index)-part of kinetic energy
  subroutine sll_s_pic_diagnostics_Hpi ( particle_group,  index, kinetic )
    class(sll_c_particle_group_base), intent(in)  :: particle_group !< particle group
    sll_int32,                        intent(in)  :: index !< velocity component
    sll_real64,                       intent(out) :: kinetic(1) !< value of \a index part of kinetic energy


    sll_real64 :: kinetic_local(1)
    sll_int32  :: i_part
    sll_real64 :: vi(3)
    sll_real64 :: wi(1)

    kinetic_local(1) = 0.0_f64
    do i_part = 1, particle_group%n_particles
       vi = particle_group%get_v(i_part)
       wi = particle_group%get_mass(i_part)
       ! Kinetic energy
       kinetic_local(1) = kinetic_local(1) + &
            (vi(index)**2)*wi(1)
    end do
    kinetic = 0.0_f64
    call sll_s_collective_reduce_real64( sll_v_world_collective, kinetic_local, 1, &
         MPI_SUM, 0, kinetic )


  end subroutine sll_s_pic_diagnostics_Hpi

  !> Compute the spline coefficient of the derivative of some given spline expansion
  subroutine sll_s_pic_diagnostics_eval_derivative_spline( position, xmin, delta_x, n_grid, field_dofs, degree, derivative )
    sll_real64, intent( in    ) :: position(:) !< particle position
    sll_real64, intent( in    ) :: xmin !< lower boundary of the domain
    sll_real64, intent( in    ) :: delta_x !< time step 
    sll_int32,  intent( in    ) :: n_grid !< number of grid points
    sll_real64, intent( in    ) :: field_dofs(:) !< coefficients of spline representation of the field
    sll_int32,  intent( in    ) :: degree !< degree of spline
    sll_real64, intent(   out ) :: derivative !< value of the derivative

    sll_int32 :: i1, der_degree, ind, index
    sll_real64 :: spline_val(degree)
    sll_real64 :: xi(3)

    der_degree = degree-1

    xi(1) = (position(1) - xmin)/delta_x
    index = ceiling(xi(1))
    xi(1) = xi(1) - real(index-1, f64)
    index = index - der_degree

    !call sll_s_uniform_b_splines_at_x( der_degree, xi(1), spline_val )
    call sll_s_uniform_bsplines_eval_basis( der_degree, xi(1), spline_val )

    derivative = 0.0_f64

    do i1 = 1, degree
       ind = modulo(index+i1-2, n_grid)+1
       derivative = derivative + spline_val(i1)*&
            (field_dofs(ind)-field_dofs(modulo(ind-2, n_grid)+1))
    end do

    derivative = derivative/delta_x


  end subroutine sll_s_pic_diagnostics_eval_derivative_spline


  !> Compute \sum(particles)w_p( v_1,p e_1(x_p) + v_2,p e_2(x_p))
  subroutine sll_s_pic_diagnostics_transfer ( particle_group, kernel_smoother_0, kernel_smoother_1, efield_dofs, transfer)
    class(sll_c_particle_group_base), intent( in   )  :: particle_group   
    class(sll_c_particle_mesh_coupling_1d) :: kernel_smoother_0  !< Kernel smoother (order p+1)
    class(sll_c_particle_mesh_coupling_1d) :: kernel_smoother_1  !< Kernel smoother (order p)   
    sll_real64, intent( in    ) :: efield_dofs(:,:) !< coefficients of efield
    sll_real64, intent(   out ) :: transfer(1) !< result

    sll_int32 :: i_part
    sll_real64 :: xi(3), vi(3), wi, efield(2), transfer_local(1)

    transfer_local = 0.0_f64
    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x( i_part )
       wi = particle_group%get_charge( i_part )
       vi = particle_group%get_v( i_part )

       call kernel_smoother_1%evaluate &
            (xi(1), efield_dofs(:,1), efield(1))
       call kernel_smoother_0%evaluate &
            (xi(1), efield_dofs(:,2), efield(2))

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

    poynting =  maxwell_solver%inner_product( efield_dofs, scratch, degree )

  end subroutine sll_s_pic_diagnostics_poynting



  !> Accumulate rho and solve Poisson
  subroutine solve_poisson( sim, rho_local, rho)
    class(sll_t_sim_pic_vm_1d2v_cart_multispecies), intent(inout) :: sim
    sll_real64,                        intent(  out) :: rho_local(:)
    sll_real64,                        intent(  out) :: rho(:)
    !loca variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: xi(3), wi(1)
    sll_real64 :: rho_cv(sim%n_gcells)

    rho_local = 0.0_f64
    do i_sp = 1, sim%particle_group%n_species
       do i_part = 1, sim%particle_group%group(i_sp)%n_particles
          xi = sim%particle_group%group(i_sp)%get_x(i_part)
          wi(1) = sim%particle_group%group(i_sp)%get_charge( i_part, sim%no_weights)
          call sim%kernel_smoother_0%add_charge(xi(1), wi(1), rho_local)
       end do
    end do
    ! MPI to sum up contributions from each processor
    rho = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, &
         rho_local, &
         sim%n_total0, MPI_SUM, rho)

    call sim%filter%apply_inplace( rho )

    if ( sim%no_weights == 3 ) then
       call sim%maxwell_solver%compute_rhs_from_function( one, sim%degree_fem, rho_cv)
       rho = rho-rho_cv
    end if

    ! Not needed here since the two control variates should cancel out each other in the x space (both are constant in x)
    !if ( sim%no_weights == 3 ) then
    !   call sim%maxwell_solver%compute_rhs_from_function( one, sim%degree_smoother, rho_cv)
    !      rho = rho-rho_cv
    !   end if
    if(sim%delta(1) == 1.0_f64) then
       rho = sim%force_sign *rho
    else
       rho = sim%force_sign *(rho - sum(sim%delta)*sim%delta_x) 
    end if
    ! Solve Poisson problem
!!$    if(sim%boundary)then
!!$       if(sim%ct) then
!!$          rho(1)=rho(1)-sin_l(0._f64)/sim%map%jacobian( [0._f64, 0._f64, 0._f64] )
!!$          rho(sim%n_total0)=rho(sim%n_total0) +sin_l(1._f64)/sim%map%jacobian( [1._f64, 0.0_f64, 0._f64])!+ rad_l(1._f64)!
!!$       else
!!$          rho(1) = rho(1) - sin_l(0._f64) !homogeneous Neumann boundary
!!$          rho(sim%n_total0) = rho(sim%n_total0) + sin_l(1._f64)!+ rad_l(1._f64)!
!!$       end if
!!$    end if

    call sim%maxwell_solver%compute_E_from_rho( rho, sim%efield_dofs(1:sim%n_total1,1) )

  contains

    function one( x)
      sll_real64 :: one
      sll_real64, intent(in) :: x

      one = 1.0_f64

    end function one

    function sin_l(eta)
      sll_real64             :: sin_l
      sll_real64, intent(in) :: eta
      !local variable
      sll_real64             :: x1
      sll_real64             :: k=0.5_f64
      sll_real64             :: alpha=0.5_f64

      if( sim%ct )then
         x1=sim%map%get_x1( [eta, 0._f64, 0._f64] )
      else
         x1 = sim%domain(3)*eta
      end if
      sin_l = -alpha/k*sin(k*x1)
    end function sin_l

  end subroutine solve_poisson


  !> Accumulate rho and solve Poisson
  subroutine check_gauss_law( sim, rho, rho_gauss, error )
    class(sll_t_sim_pic_vm_1d2v_cart_multispecies), intent(inout) :: sim
    sll_real64,                       intent(  out) :: rho(:)
    sll_real64,                       intent(  out) :: rho_gauss(:)
    sll_real64,                       intent(  out) :: error
    !local variables
    sll_int32 :: i_part, i_sp
    sll_real64 :: xi(3), wi(1)
    sll_real64 :: rho_cv(sim%n_gcells)

    rho_gauss = 0.0_f64
    do i_sp = 1, sim%particle_group%n_species
       do i_part = 1, sim%particle_group%group(i_sp)%n_particles
          xi = sim%particle_group%group(i_sp)%get_x(i_part)
          wi(1) = sim%particle_group%group(i_sp)%get_charge( i_part, sim%no_weights)
          if (wi(1) /= 0._f64) then
             call sim%kernel_smoother_0%add_charge(xi(1), wi(1), rho_gauss)
          end if
       end do
    end do
    ! MPI to sum up contributions from each processor
    rho = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, &
         rho_gauss, sim%n_total0, MPI_SUM, rho)

    call sim%filter%apply_inplace( rho )

    ! if ( sim%no_weights == 3 ) then
    !    call sim%maxwell_solver%compute_rhs_from_function( one, sim%degree_smoother, rho_cv)
    !    rho = rho-rho_cv
    ! end if


    rho = rho - sum(sim%delta)*sim%delta_x 

    call sim%maxwell_solver%compute_rho_from_e( sim%efield_dofs(1:sim%n_total1,1), rho_gauss )
    error = maxval(abs(rho-rho_gauss))

  contains

    function one( x)
      sll_real64 :: one
      sll_real64, intent(in) :: x

      one = 1.0_f64

    end function one
  end subroutine check_gauss_law



  subroutine compare_projected_bfield( kernel_smoother_0, kernel_smoother_1, maxwell_solver, maxwell_norm, particle_group, bfield_dofs )
    class(sll_c_particle_group_base), intent(in) :: particle_group
    class(sll_c_maxwell_1d_base),     intent(inout) :: maxwell_solver
    class(sll_t_maxwell_1d_fem_sm),     intent(in) :: maxwell_norm
    class(sll_c_particle_mesh_coupling_1d),     intent(inout) :: kernel_smoother_1
    class(sll_c_particle_mesh_coupling_1d),     intent(inout) :: kernel_smoother_0
    sll_real64,                       intent(in) :: bfield_dofs(:)

    sll_int32 :: i_part
    sll_real64 :: xi(3), wi(1), vi(3), vB(2), bfield
    sll_real64 :: bfield_proj(size(bfield_dofs)), scratch(size(bfield_dofs))

    vB = 0.0_f64

    call maxwell_norm%mixed_mass%dot( bfield_dofs, scratch )
    select type( maxwell_solver )
    type is ( sll_t_maxwell_1d_fem ) 
       call maxwell_solver%invert_mass( scratch,   bfield_proj, maxwell_solver%s_deg_0  )
    end select

    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x(i_part)
       vi = particle_group%get_v(i_part)
       wi(1) = particle_group%get_charge( i_part )

       call kernel_smoother_1%evaluate( xi, bfield_dofs, bfield )

       vB(1) = vB(1) + wi(1)*vi(2)*bfield 

       call kernel_smoother_0%evaluate( xi, bfield_proj, bfield )

       vB(2) = vB(2) + wi(1)*vi(2)*bfield 

    end do

    write(15,*) vB

  end subroutine compare_projected_bfield

end module sll_m_sim_pic_vm_1d2v_cart_multispecies
