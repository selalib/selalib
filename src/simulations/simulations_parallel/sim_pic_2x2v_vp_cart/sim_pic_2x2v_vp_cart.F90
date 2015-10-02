! Simulation of 2x2v Vlasov-Poisson with simple PIC method, periodic boundary conditions and Landau initial values along x1 only.

! TODO: Can be made more general

module sll_m_sim_pic_2x2v_vp_cart

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_utilities.h"

  use sll_constants
  use sll_collective
  use sll_simulation_base
  use sll_cartesian_meshes
  use sll_m_pic_base
  use sll_m_particle_initializer
  use sll_m_particle_group_2x2v
  use sll_m_kernel_smoother_base
  use sll_m_kernel_smoother_spline_2d
  use sll_m_poisson_2d_fft ! TODO: Use generic interface
  use sll_m_poisson_2d_base
  use sll_operator_splitting
  !use sll_m_operator_splitting_pic_2x2v_vp
  use sll_ascii_io

  type, extends(sll_simulation_base_class) :: sll_sim_pic_2x2v_vp_cart

     ! Abstract particle group
     class(sll_particle_group_base), pointer :: particle_group
     ! Specific particle group
     class(sll_particle_group_2x2v), pointer :: specific_particle_group 

     ! Cartesian mesh
     type(sll_cartesian_mesh_2d), pointer    :: mesh  ! [[selalib:src/meshes/sll_cartesian_meshes.F90::sll_cartesian_mesh_2d]]

     ! Abstract kernel smoother
     class(sll_kernel_smoother_base), pointer :: kernel_smoother
     ! Specific kernel smoother
     class(sll_kernel_smoother_spline_2d), pointer :: specific_kernel_smoother

     ! Poisson solver
     class(poisson_2d_fft_solver), pointer :: poisson_solver 

     ! Abstract operator splitting
     class(operator_splitting), pointer :: propagator
     ! Specific operator splitting
     class(sll_operator_splitting_pic_2x2v_vp), pointer :: specific_propagator
     
     ! Physical parameters
     sll_real64 :: landau_param(2) ! (1+landau_param(1)*cos(landau_param(2)*x1) 
     sll_real64 :: thermal_velocity(2) ! 
     
     ! Simulation parameters
     sll_real64 :: delta_t
     sll_int32  :: n_time_steps
     sll_int32  :: n_particles
     sll_int32  :: n_total_particles
     sll_int32  :: degree_smoother
     
     ! Parameters for MPI
     sll_int32  :: rank
     sll_int32  :: world_size
     
     
   contains
     procedure :: init_from_file => init_pic_2x2v
     procedure :: run => run_pic_2x2v

  end type sll_sim_pic_2x2v_vp_cart

  interface sll_delete
     module procedure delete_pic_2x2v
  end interface sll_delete
  
contains
!------------------------------------------------------------------------------!
  ! Read in the simulation parameters from input file
  subroutine init_pic_2x2v (sim, filename)
    class(sll_sim_pic_2x2v_vp_cart), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename

    sll_int32   :: n_time_steps
    sll_real64  :: delta_t, alpha, n_mode, thermal_v1, thermal_v2
    sll_int32   :: ng_x1, ng_x2
    sll_real64  :: x1_min, x1_max, x2_min, x2_max
    sll_int32   :: n_particles, degree_smoother

    sll_int32, parameter :: input_file = 99
    

    namelist /sim_params/         delta_t, n_time_steps, alpha, n_mode, thermal_v1, thermal_v2
    
    namelist /grid_dims/          ng_x1, ng_x2, x1_min, x2_min, x1_max, x2_max

    namelist /pic_params/         n_particles, degree_smoother

    ! Read parameters from file
    open(unit = input_file, file=trim(filename), IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'init_pic_2x2v() failed to open file ', filename
       STOP
    end if

    read(input_file, sim_params)
    read(input_file, grid_dims)
    read(input_file, pic_params)
    close (input_file)

    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%rank = sll_get_collective_rank(sll_world_collective)

    sim%delta_t = delta_t
    sim%n_time_steps = n_time_steps
    sim%landau_param = [alpha, n_mode*2.0_f64 * sll_pi/(x1_max - x1_min)]
    sim%thermal_velocity = [thermal_v1, thermal_v2]

    sim%mesh => new_cartesian_mesh_2d( ng_x1, ng_x2, &
         x1_min, x1_max, x2_min, x2_max)

    sim%n_particles = n_particles/sim%world_size
    sim%n_total_particles = sim%n_particles * sim%world_size
    sim%degree_smoother = degree_smoother
    

  end subroutine init_pic_2x2v

!------------------------------------------------------------------------------!

  subroutine run_pic_2x2v (sim)
    class(sll_sim_pic_2x2v_vp_cart), intent(inout) :: sim

    ! Loop variables
    sll_int32, allocatable :: rnd_seed(:)
    sll_int32 :: j, ierr
    sll_real64 :: domain(2,2)
    sll_int32 :: rnd_seed_size
    sll_real64 :: eenergy
    sll_int32 :: th_diag_id


    if (sim%rank == 0) then
       call sll_ascii_file_create('thdiag.dat', th_diag_id, ierr)
    end if

    ! Initialize the particles   
     sim%specific_particle_group => sll_new_particle_group_2x2v(sim%n_particles, &
         sim%n_total_particles ,1.0_f64, 1.0_f64)
    
    !print*, 'size', size(sim%specific_particle_group%particle_array,1)
    sim%particle_group => sim%specific_particle_group
    
    ! Set the seed for the random initialization
    call random_seed(size=rnd_seed_size)
    SLL_ALLOCATE(rnd_seed(rnd_seed_size), j)
    do j=1, rnd_seed_size
       rnd_seed(j) = (-1)**j*(100 + 15*j)*(2*sim%rank + 1)
    end do

    !print*, 'rd', rnd_seed_size

    ! Initialize the field solver
    sim%poisson_solver => new_poisson_2d_fft_solver( &
         sim%mesh%eta1_min, sim%mesh%eta1_max, sim%mesh%num_cells1, &
         sim%mesh%eta2_min, sim%mesh%eta2_max, sim%mesh%num_cells2)

    ! Initialize the kernel smoother
    domain(:,1) = [sim%mesh%eta1_min, sim%mesh%eta2_min]
    domain(:,2) = [sim%mesh%eta1_max, sim%mesh%eta2_max]
    sim%specific_kernel_smoother => sll_new_smoother_spline_2d(&
         domain, [sim%mesh%num_cells1, sim%mesh%num_cells2], sim%n_particles, &
         sim%degree_smoother)
    sim%kernel_smoother => sim%specific_kernel_smoother

    ! Initialize position and velocity of the particles.
    call sll_particle_initialize_landaux1_2x2v (sim%particle_group, sim%landau_param, &
         [sim%mesh%eta1_min, sim%mesh%eta2_min] , &
         [sim%mesh%eta1_max - sim%mesh%eta1_min, sim%mesh%eta1_max -sim%mesh%eta2_min], &
         sim%thermal_velocity, rnd_seed)

    ! Initialize the time-splitting propagator
    sim%specific_propagator => sll_new_splitting_pic_2x2v_vp(sim%poisson_solver, sim%kernel_smoother, sim%particle_group)
    sim%propagator => sim%specific_propagator


    print*, 'Time loop'
    ! Time loop
    do j=1, sim%n_time_steps

       call sim%specific_propagator%strang_splitting(sim%delta_t)

       ! Diagnostics
       if (sim%rank == 0) then
          eenergy = sum(sim%specific_propagator%efield1**2)*&
               sim%mesh%delta_eta1*sim%mesh%delta_eta2
          write(th_diag_id,'(f12.5,2g20.12)' ) real(j,f64)*sim%delta_t,  eenergy
       end if
    end do
    

  end subroutine run_pic_2x2v

!------------------------------------------------------------------------------!

  subroutine delete_pic_2x2v (sim)
    class(sll_sim_pic_2x2v_vp_cart), intent(inout) :: sim
  end subroutine delete_pic_2x2v

!------------------------------------------------------------------------------!




end module sll_m_sim_pic_2x2v_vp_cart
