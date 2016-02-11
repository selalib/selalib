! Simulation of 2d2v Vlasov-Poisson with simple PIC method, periodic boundary conditions and Landau initial values along x1 only.

! TODO: Can be made more general

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

  use sll_m_kernel_smoother_base, only: &
    sll_p_collocation, &
    sll_c_kernel_smoother

  use sll_m_kernel_smoother_spline_2d, only: &
    sll_t_kernel_smoother_spline_2d, &
    sll_s_new_kernel_smoother_spline_2d_ptr

  use sll_m_operator_splitting_pic_vp_2d2v, only: &
    sll_t_operator_splitting_pic_vp_2d2v

  use sll_m_particle_group_2d2v, only: &
    sll_s_new_particle_group_2d2v_ptr, &
    sll_t_particle_group_2d2v

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_particle_initializer, only: &
    sll_s_particle_initialize_random_landau_2d2v, &
    sll_s_particle_initialize_sobol_landau_2d2v

  use sll_m_pic_poisson_base, only : &
    sll_c_pic_poisson

  use  sll_m_pic_poisson_2d, only: &
    sll_s_new_pic_poisson_2d, &
    sll_t_pic_poisson_2d

  use sll_m_poisson_2d_fft, only: &
    sll_f_new_poisson_2d_fft_solver, &
    sll_t_poisson_2d_fft_solver

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
     class(sll_c_kernel_smoother), pointer :: kernel_smoother

     ! Poisson solver
     class(sll_c_poisson_2d_base), pointer :: poisson_solver 

     ! PIC Poisson solver
     class(sll_c_pic_poisson), pointer :: solver

     ! Abstract operator splitting
     type(sll_t_operator_splitting_pic_vp_2d2v) :: propagator
     !class(sll_t_operator_splitting), pointer :: propagator

     ! Control variate
     type(sll_t_control_variate), pointer :: control_variate
     sll_int32  :: no_weights
     
     ! Physical parameters
     sll_real64 :: landau_param(2) ! (1+landau_param(1)*cos(landau_param(2)*x1) 
     sll_real64 :: thermal_velocity(2) ! 
     
     ! Simulation parameters
     sll_real64 :: delta_t
     sll_int32  :: n_time_steps
     sll_int32  :: n_particles
     sll_int32  :: n_total_particles
     sll_int32  :: degree_smoother
     sll_int32  :: init_case
     
     ! Parameters for MPI
     sll_int32  :: rank
     sll_int32  :: world_size

     ! 
     logical    :: ctest_passed
     
     
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
    character(len=*), intent(in)                                :: filename

    sll_int32   :: n_time_steps
    sll_real64  :: delta_t, alpha, n_mode, thermal_v1, thermal_v2
    sll_int32   :: ng_x1, ng_x2
    sll_real64  :: x1_min, x1_max, x2_min, x2_max
    sll_int32   :: n_particles, degree_smoother
    character(len=256) :: init_case
    logical     :: with_control_variate

    sll_int32, parameter :: input_file = 99
    sll_int32 :: io_stat
    

    namelist /sim_params/         delta_t, n_time_steps, alpha, n_mode, thermal_v1, thermal_v2
    
    namelist /grid_dims/          ng_x1, ng_x2, x1_min, x2_min, x1_max, x2_max

    namelist /pic_params/         init_case, n_particles, degree_smoother, with_control_variate

    ! Read parameters from file
    open(unit = input_file, file=trim(filename), IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'init_pic_2d2v() failed to open file ', filename
       STOP
    end if

    read(input_file, sim_params)
    read(input_file, grid_dims)
    read(input_file, pic_params)
    close (input_file)

    sim%world_size = sll_f_get_collective_size(sll_v_world_collective)
    sim%rank = sll_f_get_collective_rank(sll_v_world_collective)

    sim%delta_t = delta_t
    sim%n_time_steps = n_time_steps
    sim%landau_param = [alpha, n_mode*2.0_f64 * sll_p_pi/(x1_max - x1_min)]
    sim%thermal_velocity = [thermal_v1, thermal_v2]

    sim%mesh => sll_f_new_cartesian_mesh_2d( ng_x1, ng_x2, &
         x1_min, x1_max, x2_min, x2_max)

    sim%n_particles = n_particles/sim%world_size
    sim%n_total_particles = sim%n_particles * sim%world_size
    sim%degree_smoother = degree_smoother
    
    select case(init_case)
    case("SLL_INIT_RANDOM")
       sim%init_case = SLL_INIT_RANDOM
    case("SLL_INIT_SOBOL")
       sim%init_case = SLL_INIT_SOBOL
    case default
       print*, '#init case ', init_case, ' not implemented.'
    end select

    if (with_control_variate .EQV. .TRUE.) then
       sim%no_weights = 3
    else
       sim%no_weights = 1
    end if

  end subroutine init_pic_2d2v

!------------------------------------------------------------------------------!

  subroutine run_pic_2d2v (sim)
    class(sll_t_sim_pic_vp_2d2v_cart), intent(inout) :: sim

    ! Loop variables
    sll_int32, allocatable :: rnd_seed(:)
    sll_int32 :: j, ierr
    sll_real64 :: domain(2,2)
    sll_int32 :: rnd_seed_size
    sll_int64 :: sobol_seed
    sll_real64 :: eenergy
    sll_int32 :: th_diag_id
    sll_real64, pointer :: control_variate_parameter(:)


    if (sim%rank == 0) then
       call sll_s_ascii_file_create('thdiag.dat', th_diag_id, ierr)
    end if

    ! Initialize the particles   
    call sll_s_new_particle_group_2d2v_ptr&
         (sim%particle_group, sim%n_particles, &
         sim%n_total_particles ,1.0_f64, 1.0_f64, sim%no_weights)
    
    
    ! Initialize control variate
    SLL_ALLOCATE(control_variate_parameter(2), ierr)
    control_variate_parameter = sim%thermal_velocity
    allocate(sim%control_variate)
    call sim%control_variate%init(control_variate_equi, &
         control_variate_parameter)


    if (sim%init_case == SLL_INIT_RANDOM) then
       ! Set the seed for the random initialization
       call random_seed(size=rnd_seed_size)
       SLL_ALLOCATE(rnd_seed(rnd_seed_size), j)
       do j=1, rnd_seed_size
          rnd_seed(j) = (-1)**j*(100 + 15*j)*(2*sim%rank + 1)
       end do

       ! Initialize position and velocity of the particles.
       ! Random initialization
       call sll_s_particle_initialize_random_landau_2d2v &
            (sim%particle_group, sim%landau_param, &
            [sim%mesh%eta1_min, sim%mesh%eta2_min] , &
            [sim%mesh%eta1_max - sim%mesh%eta1_min, sim%mesh%eta2_max -sim%mesh%eta2_min], &
            sim%thermal_velocity, rnd_seed)
    elseif (sim%init_case == SLL_INIT_SOBOL) then
       sobol_seed = int(10 + sim%rank*sim%particle_group%n_particles, 8)
       ! Pseudorandom initialization with sobol numbers
       call sll_s_particle_initialize_sobol_landau_2d2v(sim%particle_group, &
            sim%landau_param,  [sim%mesh%eta1_min, sim%mesh%eta2_min] , &
            [sim%mesh%eta1_max - sim%mesh%eta1_min, &
            sim%mesh%eta2_max -sim%mesh%eta2_min], &
            sim%thermal_velocity, sobol_seed)
    end if


    ! Initialize the field solver
    sim%poisson_solver => sll_f_new_poisson_2d_fft_solver( &
         sim%mesh%eta1_min, sim%mesh%eta1_max, sim%mesh%num_cells1, &
         sim%mesh%eta2_min, sim%mesh%eta2_max, sim%mesh%num_cells2)

    ! Initialize the kernel smoother
    domain(:,1) = [sim%mesh%eta1_min, sim%mesh%eta2_min]
    domain(:,2) = [sim%mesh%eta1_max, sim%mesh%eta2_max]
    call sll_s_new_kernel_smoother_spline_2d_ptr(sim%kernel_smoother, &
         domain, [sim%mesh%num_cells1, sim%mesh%num_cells2], sim%n_particles, &
         sim%degree_smoother, sll_p_collocation)

    ! Initialize the PIC field solver
    call sll_s_new_pic_poisson_2d(sim%solver, &
         [sim%mesh%num_cells1, sim%mesh%num_cells2], &
         sim%poisson_solver, sim%kernel_smoother)


    ! Initialize the time-splitting propagator
    if (sim%no_weights == 1) then
       call sim%propagator%init(sim%solver, sim%particle_group)
    elseif (sim%no_weights == 3) then
       call sim%propagator%init( &
            sim%solver, sim%particle_group, sim%control_variate, 3)
    end if


    print*, 'Time loop'
    ! Time loop
    do j=1, sim%n_time_steps

       call sim%propagator%strang_splitting(sim%delta_t)

       ! Diagnostics
       if (sim%rank == 0) then
          eenergy = sim%solver%compute_field_energy(1)
          write(th_diag_id,'(f12.5,2g20.12)' ) real(j,f64)*sim%delta_t,  eenergy
       end if
    end do

!!$    if (sim%rank == 0) then
!!$       select type( q => sim%solver )
!!$       type is ( sll_t_pic_poisson_2d )
!!$          print*,  q%rho_dofs
!!$       end select
!!$    end if
    

  end subroutine run_pic_2d2v

!------------------------------------------------------------------------------!

  subroutine delete_pic_2d2v (sim)
    class(sll_t_sim_pic_vp_2d2v_cart), intent(inout) :: sim
    
    call sim%solver%free()
    deallocate(sim%solver)
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

  function control_variate_equi( this, xi, vi, time) result(sll_f_control_variate)
    class(sll_t_control_variate) :: this
    sll_real64, optional,  intent( in ) :: xi(:) !< particle position
    sll_real64, optional, intent( in ) :: vi(:) !< particle velocity
    sll_real64, optional, intent( in ) :: time  !< current time
    sll_real64               :: sll_f_control_variate


    sll_f_control_variate = exp(-0.5_f64*&
         ((vi(1)/this%control_variate_parameters(1))**2+&
         (vi(2)/this%control_variate_parameters(2))**2))/&
         (2.0_f64*sll_p_pi*product(this%control_variate_parameters))

  end function control_variate_equi



end module sll_m_sim_pic_vp_2d2v_cart
