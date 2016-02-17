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

  use sll_m_kernel_smoother_base, only: &
    sll_p_collocation, &
    sll_c_kernel_smoother_base

  use sll_m_kernel_smoother_spline_2d, only: &
    sll_t_kernel_smoother_spline_2d, &
    sll_f_new_smoother_spline_2d

  use sll_m_operator_splitting, only: &
    sll_t_operator_splitting

  use sll_m_operator_splitting_pic_vp_2d2v, only: &
    sll_f_new_hamiltonian_splitting_pic_vp_2d2v, &
    sll_t_operator_splitting_pic_vp_2d2v

  use sll_m_particle_group_2d2v, only: &
    sll_f_new_particle_group_2d2v, &
    sll_t_particle_group_2d2v

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_particle_initializer, only: &
    sll_s_particle_initialize_random_landau_2d2v, &
    sll_s_particle_initialize_sobol_landau_2d2v

  use sll_m_poisson_2d_periodic_fft, only: &
    sll_f_new_poisson_2d_periodic_fft, &
    sll_t_poisson_2d_periodic_fft

  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  implicit none

  public :: &
    sll_o_delete, &
    sll_t_sim_pic_vp_2d2v_cart

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: SLL_INIT_RANDOM=0
  sll_int32, parameter :: SLL_INIT_SOBOL=1

  type, extends(sll_c_simulation_base_class) :: sll_t_sim_pic_vp_2d2v_cart

     ! Abstract particle group
     class(sll_c_particle_group_base), pointer :: particle_group
     ! Specific particle group
     class(sll_t_particle_group_2d2v), pointer :: specific_particle_group 

     ! Array for efield
     sll_real64, pointer :: efield(:,:)

     ! Cartesian mesh
     type(sll_t_cartesian_mesh_2d), pointer    :: mesh  ! [[selalib:src/meshes/sll_m_cartesian_meshes.F90::sll_t_cartesian_mesh_2d]]

     ! Abstract kernel smoother
     class(sll_c_kernel_smoother_base), pointer :: kernel_smoother
     ! Specific kernel smoother
     class(sll_t_kernel_smoother_spline_2d), pointer :: specific_kernel_smoother

     ! Poisson solver
     class(sll_t_poisson_2d_periodic_fft), pointer :: poisson_solver 

     ! Abstract operator splitting
     class(sll_t_operator_splitting), pointer :: propagator
     ! Specific operator splitting
     class(sll_t_operator_splitting_pic_vp_2d2v), pointer :: specific_propagator
     
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
     
     
   contains
     procedure :: init_from_file => init_pic_2d2v
     procedure :: run => run_pic_2d2v

  end type sll_t_sim_pic_vp_2d2v_cart

  interface sll_o_delete
     module procedure delete_pic_2d2v
  end interface sll_o_delete
  
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

    sll_int32, parameter :: input_file = 99
    sll_int32 :: io_stat
    

    namelist /sim_params/         delta_t, n_time_steps, alpha, n_mode, thermal_v1, thermal_v2
    
    namelist /grid_dims/          ng_x1, ng_x2, x1_min, x2_min, x1_max, x2_max

    namelist /pic_params/         init_case, n_particles, degree_smoother

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


    if (sim%rank == 0) then
       call sll_s_ascii_file_create('thdiag.dat', th_diag_id, ierr)
    end if

    ! Initialize the particles   
     sim%specific_particle_group => sll_f_new_particle_group_2d2v(sim%n_particles, &
         sim%n_total_particles ,1.0_f64, 1.0_f64, 1)
    
    !print*, 'size', size(sim%specific_particle_group%particle_array,1)
    sim%particle_group => sim%specific_particle_group
    

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
       sobol_seed = 10_8 + sim%rank*sim%particle_group%n_particles
       ! Pseudorandom initialization with sobol numbers
       call sll_s_particle_initialize_sobol_landau_2d2v(sim%particle_group, &
            sim%landau_param,  [sim%mesh%eta1_min, sim%mesh%eta2_min] , &
            [sim%mesh%eta1_max - sim%mesh%eta1_min, sim%mesh%eta2_max -sim%mesh%eta2_min], &
            sim%thermal_velocity, sobol_seed)
    end if


    !print*, 'rd', rnd_seed_size

    ! Initialize the field solver
    sim%poisson_solver => sll_f_new_poisson_2d_periodic_fft( &
         sim%mesh%eta1_min, sim%mesh%eta1_max, sim%mesh%num_cells1, &
         sim%mesh%eta2_min, sim%mesh%eta2_max, sim%mesh%num_cells2)

    ! Initialize the kernel smoother
    domain(:,1) = [sim%mesh%eta1_min, sim%mesh%eta2_min]
    domain(:,2) = [sim%mesh%eta1_max, sim%mesh%eta2_max]
    sim%specific_kernel_smoother => sll_f_new_smoother_spline_2d(&
         domain, [sim%mesh%num_cells1, sim%mesh%num_cells2], sim%n_particles, &
         sim%degree_smoother, sll_p_collocation)
    sim%kernel_smoother => sim%specific_kernel_smoother


    ! Initialize the time-splitting propagator
    SLL_ALLOCATE(sim%efield(sim%kernel_smoother%n_dofs,2),ierr)
    sim%specific_propagator => sll_f_new_hamiltonian_splitting_pic_vp_2d2v(sim%poisson_solver, sim%kernel_smoother, sim%particle_group, sim%efield)
    sim%propagator => sim%specific_propagator


    print*, 'Time loop'
    ! Time loop
    do j=1, sim%n_time_steps

       call sim%specific_propagator%strang_splitting(sim%delta_t)

       ! Diagnostics
       if (sim%rank == 0) then
          eenergy = sum(sim%efield(:,1)**2)*&
               sim%mesh%delta_eta1*sim%mesh%delta_eta2
          write(th_diag_id,'(f12.5,2g20.12)' ) real(j,f64)*sim%delta_t,  eenergy
       end if
    end do
    

  end subroutine run_pic_2d2v

!------------------------------------------------------------------------------!

  subroutine delete_pic_2d2v (sim)
    class(sll_t_sim_pic_vp_2d2v_cart), intent(inout) :: sim
  end subroutine delete_pic_2d2v

!------------------------------------------------------------------------------!




end module sll_m_sim_pic_vp_2d2v_cart
