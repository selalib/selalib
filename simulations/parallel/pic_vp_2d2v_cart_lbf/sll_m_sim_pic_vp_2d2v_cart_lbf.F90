!> @ingroup particle_methods
!> @author MCP ALH KK
!> @brief Simulation of 2d2v Vlasov-Poisson with several PIC methods (including the lbf-resampled pic method),
!>  periodic boundary conditions and Landau initial values along x1 only.

! MCP: I am writing this simulation from Katharina's sll_m_sim_pic_vp_2d2v_cart simulation, to add the lbf particle group

module sll_m_sim_pic_vp_2d2v_cart_lbf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_working_precision.h"

  use sll_m_ascii_io, only: &
    sll_s_ascii_file_create

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d

  use sll_m_collective, only: &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_control_variate, only: &
    sll_t_control_variate

  use sll_m_particle_mesh_coupling_base, only: &
    sll_p_collocation, &
    sll_c_particle_mesh_coupling

  use sll_m_particle_mesh_coupling_spline_2d, only: &
    sll_t_particle_mesh_coupling_spline_2d, &
    sll_s_new_particle_mesh_coupling_spline_2d_ptr

  use sll_m_operator_splitting_pic_vp_2d2v, only: &
    sll_t_operator_splitting_pic_vp_2d2v

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_particle_group_2d2v, only: &
    sll_s_new_particle_group_2d2v_ptr

  use sll_m_particle_group_2d2v_lbf, only: &
    sll_s_new_particle_group_2d2v_lbf_ptr

  use sll_m_initial_distribution, only : &
       sll_c_distribution_params, &
       sll_s_initial_distribution_new

  use sll_m_particle_visualization_interface, only : &
    sll_t_plotting_params_2d, &
    sll_s_visualize_particle_group

  use sll_m_particle_sampling_interface, only : &
    sll_s_sample_particle_group, &
    sll_s_resample_particle_group

  use sll_m_particle_sampling, only : &
       sll_t_particle_sampling

  use sll_m_pic_poisson_base, only : &
    sll_c_pic_poisson

  use  sll_m_pic_poisson_2d, only: &
    sll_s_new_pic_poisson_2d, &
    sll_t_pic_poisson_2d

  use sll_m_poisson_2d_periodic, only: &
    sll_f_new_poisson_2d_periodic

  use sll_m_poisson_2d_base, only: &
    sll_c_poisson_2d_base

  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  implicit none

  public :: &
    sll_t_sim_pic_vp_2d2v_cart_lbf, &
    sll_p_simple_particles,   &
    sll_p_lbf_particles

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: sll_p_simple_particles=0
  sll_int32, parameter :: sll_p_lbf_particles=1

  type, extends(sll_c_simulation_base_class) :: sll_t_sim_pic_vp_2d2v_cart_lbf

     ! Abstract particle group
     class(sll_c_particle_group_base), pointer :: particle_group

     ! Array for efield
     !sll_real64, pointer :: efield(:,:)

     ! Cartesian mesh
     type(sll_t_cartesian_mesh_2d), pointer    :: mesh  ! [[selalib:src/meshes/sll_m_cartesian_meshes.F90::sll_t_cartesian_mesh_2d]]

     ! Abstract kernel smoother
     class(sll_c_particle_mesh_coupling), pointer :: kernel_smoother

     ! Poisson solver
     class(sll_c_poisson_2d_base), pointer :: poisson_solver 

     ! PIC Poisson solver
     class(sll_c_pic_poisson), pointer :: pic_poisson

     ! Abstract operator splitting
     type(sll_t_operator_splitting_pic_vp_2d2v) :: propagator

     ! Control variate
     type(sll_t_control_variate), pointer :: control_variate
     sll_int32  :: no_weights
     
     ! Physical parameters
     class(sll_c_distribution_params), allocatable :: init_distrib_params

     ! Particle sampling (for random sampling of simple particles)
     type(sll_t_particle_sampling) :: sampler

     ! Simulation parameters
     sll_real64 :: delta_t
     sll_int32  :: n_time_steps
     sll_int32  :: resample_period
     sll_int32  :: n_particles
     sll_int32  :: n_total_particles
     sll_int32  :: degree_smoother
     sll_int32  :: particle_type
     !     sll_int32  :: init_case

     ! Parameters for MPI
     sll_int32  :: rank
     sll_int32  :: world_size

     ! 
     logical    :: ctest_passed = .false.
     
     
   contains
     procedure :: init_from_file => init_pic_2d2v
     procedure :: run => run_pic_2d2v
     procedure :: delete => delete_pic_2d2v

  end type sll_t_sim_pic_vp_2d2v_cart_lbf

  
contains
!------------------------------------------------------------------------------!
  ! Read in the simulation parameters from input file
  subroutine init_pic_2d2v (sim, filename)
    class(sll_t_sim_pic_vp_2d2v_cart_lbf), intent(inout) :: sim
    character(len=*),                  intent(in)    :: filename

    sll_int32   :: n_time_steps
    sll_real64  :: delta_t !, alpha, n_mode, thermal_v1, thermal_v2
    sll_real64  :: species_charge
    sll_real64  :: species_mass

    sll_int32   :: ng_x1, ng_x2
    sll_real64  :: x1_min, x1_max, x2_min, x2_max
    sll_int32   :: n_particles, degree_smoother
    character(len=256) :: particle_type_str  ! simple (std) pic, or deterministic particles with lbf resampling
    ! character(len=256) :: init_case_str

    character(len=256) :: initial_distrib
    character(len=256) :: sampling_case

    logical     :: with_control_variate

    logical   :: domain_is_x_periodic   ! needed for the LBF particles
    logical   :: domain_is_y_periodic   ! needed for the LBF particles

    sll_int32 :: input_file ! unit for nml file
    sll_int32 :: io_stat
    sll_real64 :: domain(2,2)

    ! parameters for the lbf particle group
    sll_int32   :: n_particles_x
    sll_int32   :: n_particles_y
    sll_int32   :: n_particles_vx
    sll_int32   :: n_particles_vy
    sll_int32   :: remap_period
    sll_int32   :: remap_degree
    sll_real64  :: remapping_grid_vx_min
    sll_real64  :: remapping_grid_vx_max
    sll_real64  :: remapping_grid_vy_min
    sll_real64  :: remapping_grid_vy_max
    sll_int32   :: remapping_sparse_grid_max_level_x
    sll_int32   :: remapping_sparse_grid_max_level_y
    sll_int32   :: remapping_sparse_grid_max_level_vx
    sll_int32   :: remapping_sparse_grid_max_level_vy
    sll_int32,  dimension(4)  :: remapping_sparse_grid_max_levels
    sll_real64, dimension(4)  :: remapping_grid_eta_min
    sll_real64, dimension(4)  :: remapping_grid_eta_max

    ! this namelist should follows std terminology
    namelist /sim_params/         delta_t, n_time_steps, initial_distrib

    ! this namelist should follows std terminology
    namelist /grid_dims/          ng_x1, ng_x2, x1_min, x2_min, x1_max, x2_max

    ! this namelist should follows std terminology
    namelist /pic_params/         n_particles, degree_smoother, sampling_case, with_control_variate   ! for the 'simple' pic method

    ! namelist specific to particle polymorphism
    namelist /particle_method/    particle_type_str

              !     namelist /sim_params/         delta_t, n_time_steps, alpha, n_mode, thermal_v1, thermal_v2
              ! namelist /pic_poisson_params/ ng_x1, ng_x2, x1_min, x2_min, x1_max, x2_max, degree_smoother   ! previous name was grid_dims
              ! namelist /simple_pic_params/  init_case_str, n_particles, with_control_variate

    ! namelist specific to lbf particles: 2 sets of parameters  -------------------------------------------------
    ! 1. remapping (sparse grid) parameters
    namelist /pic_lbf_remap_params/    &
                                  remap_period,                       &
                                  remap_degree,                       &
                                  remapping_grid_vx_min,              &
                                  remapping_grid_vx_max,              &
                                  remapping_grid_vy_min,              &
                                  remapping_grid_vy_max,              &
                                  remapping_sparse_grid_max_level_x,  &
                                  remapping_sparse_grid_max_level_y,  &
                                  remapping_sparse_grid_max_level_vx, &
                                  remapping_sparse_grid_max_level_vy
    ! 2. particle parameters (initially located on a 4d cartesian grid) used to approximate the charge and the flow
    namelist /pic_lbf_particles_params/                               &
                                  n_particles_x,                      &
                                  n_particles_y,                      &
                                  n_particles_vx,                     &
                                  n_particles_vy

    ! Read parameters from file
    open(newunit = input_file, file=trim(filename), IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'init_pic_2d2v() failed to open file ', filename
       STOP
    end if

    read(input_file, sim_params)
    read(input_file, grid_dims)
    call sll_s_initial_distribution_new( trim(initial_distrib), [2,2], input_file, sim%init_distrib_params )
    read(input_file, pic_params)

    read(input_file, particle_method)
    select case(particle_type_str)

    case("SLL_SIMPLE_PARTICLES")
       sim%particle_type = sll_p_simple_particles

    case("SLL_LBF_PARTICLES")
       sim%particle_type = sll_p_lbf_particles
       read(input_file, pic_lbf_remap_params)
       read(input_file, pic_lbf_particles_params)

    case default
       SLL_ERROR('sll_m_sim_pic_vp_2d2v_cart_lbf%init_pic_2d2v', 'particle_type' // particle_type_str // ' not implemented.')
    end select

    close (input_file)

    species_charge = 1.0_f64
    species_mass = 1.0_f64

    sim%world_size = sll_f_get_collective_size(sll_v_world_collective)
    sim%rank = sll_f_get_collective_rank(sll_v_world_collective)

    sim%delta_t = delta_t
    sim%n_time_steps = n_time_steps
        !    sim%landau_param = [alpha, n_mode*2.0_f64 * sll_p_pi/(x1_max - x1_min)]
        !     sim%thermal_velocity = [thermal_v1, thermal_v2]

    sim%mesh => sll_f_new_cartesian_mesh_2d( ng_x1, ng_x2, &
         x1_min, x1_max, x2_min, x2_max)

    sim%degree_smoother = degree_smoother

    ! Initialize the particle group
    select case(sim%particle_type)

    case( sll_p_simple_particles )

      ! here we follow Katharina's initialization in sll_m_sim_pic_vp_2d2v_cart

      call sim%sampler%init( trim(sampling_case),  [2,2], sim%n_particles, sim%rank  )

      sim%n_particles = n_particles/sim%world_size
      sim%n_total_particles = sim%n_particles * sim%world_size

      if (with_control_variate .EQV. .TRUE.) then
         sim%no_weights = 3
      else
         sim%no_weights = 1
      end if

      ! Initialize the particles
      call sll_s_new_particle_group_2d2v_ptr&
           (sim%particle_group, sim%n_particles, &
           sim%n_total_particles, species_charge, species_mass, sim%no_weights)

    case( sll_p_lbf_particles )

      sim%no_weights = 1   ! no control variate for the moment with lbf particles
      domain_is_x_periodic = .true.
      domain_is_y_periodic = .true.
      remapping_grid_eta_min(1) = x1_min
      remapping_grid_eta_min(2) = x2_min
      remapping_grid_eta_min(3) = remapping_grid_vx_min
      remapping_grid_eta_min(4) = remapping_grid_vy_min
      remapping_grid_eta_max(1) = x1_max
      remapping_grid_eta_max(2) = x2_max
      remapping_grid_eta_max(3) = remapping_grid_vx_max
      remapping_grid_eta_max(4) = remapping_grid_vy_max
      remapping_sparse_grid_max_levels(1) = remapping_sparse_grid_max_level_x
      remapping_sparse_grid_max_levels(2) = remapping_sparse_grid_max_level_y
      remapping_sparse_grid_max_levels(3) = remapping_sparse_grid_max_level_vx
      remapping_sparse_grid_max_levels(4) = remapping_sparse_grid_max_level_vy
      sim%resample_period = remap_period
      call sll_s_new_particle_group_2d2v_lbf_ptr( &
          sim%particle_group, &
          species_charge,    &
          species_mass,      &
          domain_is_x_periodic,    &
          domain_is_y_periodic,    &
          remap_degree,    &
          remapping_grid_eta_min, &
          remapping_grid_eta_max, &
          remapping_sparse_grid_max_levels, &
          n_particles_x,  &
          n_particles_y,  &
          n_particles_vx, &
          n_particles_vy &
      )
      SLL_ASSERT( sim%world_size == 1 )
      sim%n_particles = sim%particle_group%n_particles
      sim%n_total_particles = sim%n_particles * sim%world_size

    case default
      SLL_ERROR('sll_m_sim_pic_vp_2d2v_cart_lbf%init_pic_2d2v', 'this particle_type is not implemented.')
    end select

    ! Initialize control variate
    allocate(sim%control_variate)
    call sim%control_variate%init(control_variate_equi, &
         distribution_params=sim%init_distrib_params)

    ! Initialize the field solver
    sim%poisson_solver => sll_f_new_poisson_2d_periodic( &
         sim%mesh%eta1_min, sim%mesh%eta1_max, sim%mesh%num_cells1, &
         sim%mesh%eta2_min, sim%mesh%eta2_max, sim%mesh%num_cells2)

    ! Initialize the kernel smoother
    domain(:,1) = [sim%mesh%eta1_min, sim%mesh%eta2_min]
    domain(:,2) = [sim%mesh%eta1_max, sim%mesh%eta2_max]
    call sll_s_new_particle_mesh_coupling_spline_2d_ptr(sim%kernel_smoother, &
         domain, [sim%mesh%num_cells1, sim%mesh%num_cells2], sim%n_particles, &
         sim%degree_smoother, sll_p_collocation)

    ! Initialize the PIC field solver
    call sll_s_new_pic_poisson_2d(sim%pic_poisson, &
         [sim%mesh%num_cells1, sim%mesh%num_cells2], &
         sim%poisson_solver, sim%kernel_smoother)

    ! Initialize the time-splitting propagator
    if (sim%no_weights == 1) then
       call sim%propagator%init(sim%pic_poisson, sim%particle_group)
    elseif (sim%no_weights == 3) then
       SLL_ASSERT( sim%particle_type == sll_p_simple_particles )
       call sim%propagator%init( &
            sim%pic_poisson, sim%particle_group, sim%control_variate, 3)
    end if


  end subroutine init_pic_2d2v

!------------------------------------------------------------------------------!

  subroutine run_pic_2d2v (sim)
    class(sll_t_sim_pic_vp_2d2v_cart_lbf), intent(inout) :: sim
    type(sll_t_plotting_params_2d)          :: plotting_params_2d
    sll_int32 :: plot_np_x    !< nb of points in the x  plotting grid for a (x,vx) plot
    sll_int32 :: plot_np_vx   !< nb of points in the vx plotting grid for a (x,vx) plot
    sll_real64 :: slice_y
    sll_real64 :: slice_vy

    ! Loop variables
    sll_int32  :: nt, ierr
    sll_real64 :: eenergy
    sll_int32  :: th_diag_id

    if (sim%rank == 0) then
       call sll_s_ascii_file_create('thdiag.dat', th_diag_id, ierr)
    end if

    call sll_s_sample_particle_group( &
        sim%particle_group,   &
        distribution_params = sim%init_distrib_params, &
        random_sampler = sim%sampler, &
        nb_weights = sim%no_weights, &
        control_variate = sim%control_variate, &
        xmin = [sim%mesh%eta1_min, sim%mesh%eta2_min], &
        Lx = [sim%mesh%eta1_max - sim%mesh%eta1_min, sim%mesh%eta2_max -sim%mesh%eta2_min] )

    plot_np_x  = 100
    plot_np_vx = 100
    slice_y  = 0.0_f64
    slice_vy  = 0.0_f64

    call plotting_params_2d%reset_params( &
        'f_slice', &
        plot_np_x, &
        plot_np_vx, &
        slice_y, &
        slice_vy &
        )

    ! Time loop
    print*, 'Time loop'
    do nt=1, sim%n_time_steps

       print*, 'time step ', nt, '/', sim%n_time_steps, ' ... '

       ! particle resampling (conditional)
       if( (sim%particle_type == sll_p_lbf_particles)  &
           .and. (modulo(nt-1, sim%resample_period) == 0 ) )then

         print *, "-- plotting f slice in gnuplot format..."
         call sll_s_visualize_particle_group( sim%particle_group, plotting_params_2d, nt)
         print*, "-- plotting done."

         if( nt > 1 )then
           print *, "-- particle resampling with deterministic LBF method..."
           call sll_s_resample_particle_group( sim%particle_group )
           print*, "-- resampling done."
         end if

       end if

       ! time integration
       call sim%propagator%strang_splitting(sim%delta_t)

       ! Diagnostics
       if (sim%rank == 0) then
          eenergy = sim%pic_poisson%compute_field_energy(1)
          print*, "-- diagnostics: electric energy = ", eenergy
          write(th_diag_id,'(f12.5,2g20.12)' ) real(nt,f64)*sim%delta_t,  eenergy
       end if

    end do


    !!! Part for ctest
    if (sim%rank == 0) then
       if (abs(eenergy - 3.0503207170668825_f64) < 1d-13) then   !!! MCP: why this value ??
          sim%ctest_passed = .true.
       end if
    end if


  end subroutine run_pic_2d2v

!------------------------------------------------------------------------------!

  subroutine delete_pic_2d2v (sim)
    class(sll_t_sim_pic_vp_2d2v_cart_lbf), intent(inout) :: sim
    
    call sim%pic_poisson%free()
    deallocate(sim%pic_poisson)
    call sim%particle_group%free()
    deallocate (sim%particle_group)
    call sim%mesh%delete()
    deallocate(sim%mesh)
    call sim%poisson_solver%free()
    deallocate(sim%poisson_solver)
    call sim%kernel_smoother%free()
    deallocate(sim%kernel_smoother)
    call sim%control_variate%free()
    deallocate(sim%control_variate)
    call sim%sampler%free()
    call sim%init_distrib_params%free()
    deallocate(sim%init_distrib_params)

  end subroutine delete_pic_2d2v

!------------------------------------------------------------------------------!

  ! copied from Katharina's simulation:
  ! As a control variate, we use the equilibrium (v part of the initial distribution)
  function control_variate_equi( this, xi, vi, time) result(sll_f_control_variate)
    class(sll_t_control_variate) :: this
    sll_real64, optional,  intent( in ) :: xi(:) !< particle position
    sll_real64, optional,  intent( in ) :: vi(:) !< particle velocity
    sll_real64, optional,  intent( in ) :: time  !< current time
    sll_real64               :: sll_f_control_variate

    sll_f_control_variate = &
         this%control_variate_distribution_params%eval_v_density( vi(1:2) )

  end function control_variate_equi


end module sll_m_sim_pic_vp_2d2v_cart_lbf
