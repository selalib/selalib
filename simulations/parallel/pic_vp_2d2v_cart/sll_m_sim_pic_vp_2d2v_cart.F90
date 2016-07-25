! Simulation of 2d2v Vlasov-Poisson with simple PIC method, periodic boundary conditions and Landau initial values along x1 only.

module sll_m_sim_pic_vp_2d2v_cart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
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

  use sll_m_initial_distribution, only : &
       sll_c_distribution_params, &
       sll_s_initial_distribution_new

  use sll_m_particle_mesh_coupling_base, only: &
    sll_p_collocation, &
    sll_c_particle_mesh_coupling

  use sll_m_particle_mesh_coupling_spline_2d, only: &
    sll_t_particle_mesh_coupling_spline_2d, &
    sll_s_new_particle_mesh_coupling_spline_2d_ptr

  use sll_m_operator_splitting_pic_vp_2d2v, only: &
    sll_t_operator_splitting_pic_vp_2d2v

  use sll_m_particle_group_2d2v, only: &
    sll_s_new_particle_group_2d2v_ptr, &
    sll_t_particle_group_2d2v

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base
  
  use sll_m_particle_sampling, only : &
       sll_t_particle_sampling, &
       sll_p_particle_sampling_sobol_symmetric, &
       sll_p_particle_sampling_sobol, &      
       sll_p_particle_sampling_random_symmetric, &
       sll_p_particle_sampling_random

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
    sll_t_sim_pic_vp_2d2v_cart

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: SLL_INIT_RANDOM=0
  sll_int32, parameter :: SLL_INIT_SOBOL=1

  type, extends(sll_c_simulation_base_class) :: sll_t_sim_pic_vp_2d2v_cart

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
     class(sll_c_distribution_params), allocatable :: params
     type(sll_t_particle_sampling) :: sampler
     
     ! Simulation parameters
     sll_real64 :: delta_t
     sll_int32  :: n_time_steps
     sll_int32  :: n_particles
     sll_int32  :: n_total_particles
     sll_int32  :: degree_smoother
     
     ! Parameters for MPI
     sll_int32  :: rank
     sll_int32  :: world_size

     ! 
     logical    :: ctest_passed = .false.
     
     
   contains
     procedure :: init_from_file => init_pic_2d2v
     procedure :: run => run_pic_2d2v
     procedure :: delete => delete_pic_2d2v

  end type sll_t_sim_pic_vp_2d2v_cart

  
contains
!------------------------------------------------------------------------------!
  ! Read in the simulation parameters from input file
  subroutine init_pic_2d2v (sim, filename)
    class(sll_t_sim_pic_vp_2d2v_cart), intent(inout) :: sim
    character(len=*),                  intent(in)    :: filename

    sll_int32   :: n_time_steps
    sll_real64  :: delta_t
    sll_int32   :: ng_x1, ng_x2
    sll_real64  :: x1_min, x1_max, x2_min, x2_max
    sll_int32   :: n_particles, degree_smoother
    character(len=256) :: initial_distrib
    character(len=256) :: sampling_case
    logical     :: with_control_variate

    sll_int32 :: input_file ! unit for nml file
    sll_int32 :: io_stat
    sll_real64 :: domain(2,2)
    

    namelist /sim_params/         delta_t, n_time_steps, initial_distrib
    
    namelist /grid_dims/          ng_x1, ng_x2, x1_min, x2_min, x1_max, x2_max

    namelist /pic_params/         n_particles, degree_smoother, sampling_case, with_control_variate

    ! Read parameters from file
    open(newunit = input_file, file=trim(filename), IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'init_pic_2d2v() failed to open file ', filename
       STOP
    end if

    read(input_file, sim_params)
    read(input_file, grid_dims)
    call sll_s_initial_distribution_new( trim(initial_distrib), [2,2], input_file, sim%params )
    read(input_file, pic_params)
    
    close (input_file)

    sim%world_size = sll_f_get_collective_size(sll_v_world_collective)
    sim%rank = sll_f_get_collective_rank(sll_v_world_collective)

    sim%delta_t = delta_t
    sim%n_time_steps = n_time_steps

    sim%mesh => sll_f_new_cartesian_mesh_2d( ng_x1, ng_x2, &
         x1_min, x1_max, x2_min, x2_max)

    sim%n_particles = n_particles/sim%world_size
    sim%n_total_particles = sim%n_particles * sim%world_size
    sim%degree_smoother = degree_smoother

    select case(sampling_case)
    case("particle_sampling_random")
       call sim%sampler%init( sll_p_particle_sampling_random, [2,2], sim%n_particles, sim%rank  )
    case("particle_sampling_sobol")     
       call sim%sampler%init( sll_p_particle_sampling_sobol, [2,2], sim%n_particles, sim%rank  )
    case("particle_sampling_random_symmetric")
       call sim%sampler%init( sll_p_particle_sampling_random_symmetric, [2,2], sim%n_particles, sim%rank  )
    case("particle_sampling_sobol_symmetric")
       call sim%sampler%init( sll_p_particle_sampling_sobol_symmetric, [2,2], sim%n_particles, sim%rank  )
    case default
       print*, '#sampling_case ', sampling_case, ' not implemented.'
    end select

    if (with_control_variate .EQV. .TRUE.) then
       sim%no_weights = 3
    else
       sim%no_weights = 1
    end if

  ! Initialize the particles   
    call sll_s_new_particle_group_2d2v_ptr&
         (sim%particle_group, sim%n_particles, &
         sim%n_total_particles ,1.0_f64, 1.0_f64, sim%no_weights)
    
    
    ! Initialize control variate
    allocate(sim%control_variate)
    call sim%control_variate%init(control_variate_equi, &
         distribution_params=sim%params)



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
       call sim%propagator%init( &
            sim%pic_poisson, sim%particle_group, sim%control_variate, 3)
    end if


  end subroutine init_pic_2d2v

!------------------------------------------------------------------------------!

  subroutine run_pic_2d2v (sim)
    class(sll_t_sim_pic_vp_2d2v_cart), intent(inout) :: sim

    ! Local variables
    sll_int32 :: j, ierr
    sll_real64 :: eenergy
    sll_int32 :: th_diag_id


    if (sim%rank == 0) then
       call sll_s_ascii_file_create('thdiag.dat', th_diag_id, ierr)
    end if

    if (sim%no_weights == 1 ) then
       call sim%sampler%sample ( sim%particle_group, sim%params, [sim%mesh%eta1_min, sim%mesh%eta2_min] , &
            [sim%mesh%eta1_max - sim%mesh%eta1_min, sim%mesh%eta2_max -sim%mesh%eta2_min] )
    elseif ( sim%no_weights == 3 ) then
       call sim%sampler%sample_cv ( sim%particle_group, sim%params, [sim%mesh%eta1_min, sim%mesh%eta2_min] , &
            [sim%mesh%eta1_max - sim%mesh%eta1_min, sim%mesh%eta2_max -sim%mesh%eta2_min], sim%control_variate )
    end if

    print*, 'Time loop'
    ! Time loop
    do j=1, sim%n_time_steps

       call sim%propagator%strang_splitting(sim%delta_t)

       ! Diagnostics
       if (sim%rank == 0) then
          eenergy = sim%pic_poisson%compute_field_energy(1)
          write(th_diag_id,'(f12.5,2g20.12)' ) real(j,f64)*sim%delta_t,  eenergy
       end if
    end do




!!! Part for ctest
    if (sim%rank == 0) then
       if (abs(eenergy - 3.0503207170668825_f64) < 1d-13) then
          sim%ctest_passed = .true.
       end if
    end if

    

  end subroutine run_pic_2d2v

!------------------------------------------------------------------------------!

  subroutine delete_pic_2d2v (sim)
    class(sll_t_sim_pic_vp_2d2v_cart), intent(inout) :: sim
    
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

  end subroutine delete_pic_2d2v

!------------------------------------------------------------------------------!

  !> As a control variate, we use the equilibrium (v part of the initial distribution)
  function control_variate_equi( this, xi, vi, time) result(sll_f_control_variate)
    class(sll_t_control_variate) :: this
    sll_real64, optional,  intent( in ) :: xi(:) !< particle position
    sll_real64, optional,  intent( in ) :: vi(:) !< particle velocity
    sll_real64, optional,  intent( in ) :: time  !< current time
    sll_real64               :: sll_f_control_variate


    sll_f_control_variate = this%control_variate_distribution_params%evalv( vi(1:2) )


  end function control_variate_equi


end module sll_m_sim_pic_vp_2d2v_cart
