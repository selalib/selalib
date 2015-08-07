! Simulation of 1d2v Vlasov-Maxwell with simple PIC method, periodic boundary conditions, Weibel instability.

! TODO: Can be made more general

module sll_m_sim_pic_1d2v_vm_cart

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_utilities.h"

  use sll_constants
  use sll_collective
  use sll_simulation_base
  use sll_cartesian_meshes
  use sll_m_pic_maxwell_base
  use sll_m_pic_maxwell_1d_fourier_modes
  use sll_module_pic_base
  use sll_m_particle_initializer
  use sll_m_particle_group_1d2v
  use sll_m_operator_splitting_pic_1d2v_vm
  use sll_ascii_io

  type, extends(sll_simulation_base_class) :: sll_sim_pic_1d2v_vm_cart

     ! Abstract particle group
     class(sll_particle_group_base), pointer :: particle_group
     ! Specific particle group
     class(sll_particle_group_1d2v), pointer :: specific_particle_group 

     ! Cartesian mesh
     type(sll_cartesian_mesh_1d), pointer    :: mesh 

     ! PIC Maxwell solver with kernel smoothing 
     ! Abstract 
     class(sll_pic_maxwell_base), pointer :: maxwell_pic
     ! Specific
     class(sll_pic_maxwell_1d_fourier_modes), pointer :: specific_maxwell_pic


     ! Specific operator splitting
     class(sll_operator_splitting_pic_1d2v_vm), pointer :: propagator
     
     ! Physical parameters
     sll_real64 :: landau_param(2) ! (1+landau_param(1)*cos(landau_param(2)*x1) 
     sll_real64 :: thermal_velocity(2) ! 
     sll_real64 :: domain(3) ! x_min, x_max, Lx
     
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
     procedure :: init_from_file => init_pic_1d2v_vm
     procedure :: run => run_pic_1d2v_vm

  end type sll_sim_pic_1d2v_vm_cart

  interface sll_delete
     module procedure delete_pic_1d2v_vm
  end interface sll_delete
  
contains
!------------------------------------------------------------------------------!
  ! Read in the simulation parameters from input file
  subroutine init_pic_1d2v_vm (sim, filename)
    class(sll_sim_pic_1d2v_vm_cart), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename

    sll_int32   :: n_time_steps
    sll_real64  :: delta_t, alpha, n_mode, thermal_v, T_r
    sll_int32   :: ng_x
    sll_real64  :: x1_min, x1_max
    sll_int32   :: n_particles, degree_smoother

    sll_int32, parameter :: input_file = 99
    

    namelist /sim_params/         delta_t, n_time_steps, alpha, n_mode, thermal_v, T_r
    
    namelist /grid_dims/          ng_x, x1_min, x1_max

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
    sim%thermal_velocity = [thermal_v, thermal_v* sqrt(T_r)]

    sim%mesh => new_cartesian_mesh_1d( ng_x1, &
         x1_min, x1_max)
    sim%domain = [x1_min, x1_max, x1_max - x1_min ]

    sim%n_particles = n_particles/sim%world_size
    sim%n_total_particles = sim%n_particles * sim%world_size
    sim%degree_smoother = degree_smoother
    

  end subroutine init_pic_1d2v_vm

!------------------------------------------------------------------------------!

  subroutine run_pic_1d2v_vm (sim)
    class(sll_sim_pic_1d2v_vm_cart), intent(inout) :: sim

    ! Loop variables
    sll_int32, allocatable :: rnd_seed(:)
    sll_int32 :: j, ierr
    sll_int32 :: rnd_seed_size
    sll_real64 :: eenergy
    sll_int32 :: th_diag_id


    if (sim%rank == 0) then
       call sll_ascii_file_create('thdiag.dat', th_diag_id, ierr)
    end if

    ! Initialize the particles   
     sim%specific_particle_group => sll_new_particle_group_1d2v(sim%n_particles, &
         sim%n_total_particles ,1.0_f64, 1.0_f64)
    
    !print*, 'size', size(sim%specific_particle_group%particle_array,1)
    sim%particle_group => sim%specific_particle_group
    
    ! Set the seed for the random initialization
    call random_seed(size=rnd_seed_size)
    SLL_ALLOCATE(rnd_seed(rnd_seed_size), j)
    do j=1, rnd_seed_size
       rnd_seed(j) = (-1)**j*(100 + 15*j)*(2*sim%rank + 1)
    end do

    ! Initialize the field solver
    sim%specific_maxwell_pic => sll_new_pic_maxwell_1d_fourier_modes(&
         sim%degree_smoother, sim%domain(3), sim%n_particles)
    sim%maxwell_pic => sim%specific_maxwell_pic

    ! Initialize position and velocity of the particles.
    call sll_particle_initialize_landau_1d2v (sim%particle_group, sim%landau_param, &
         sim%domain(1) , &
         sim%domain(3), &
         sim%thermal_velocity, rnd_seed)

    ! Initialize the time-splitting propagator
    sim%propagator => sll_new_splitting_pic_1d2v_vm(sim%maxwell_pic, &
         sim%particle_group)


    write(th_diag_id,'(f12.5,2g20.12,2g20.12,2g20.12)' ) &
         0.0_f64,  sim%propagator%efield_dofs(2,:), &
         sim%propagator%bfield(2)
    print*, 'Time loop'
    ! Time loop
    do j=1, sim%n_time_steps
       ! Lee splitting
       call sim%propagator%operatorHf(sim%delta_t)
       call sim%propagator%operatorHE(sim%delta_t)
       call sim%propagator%operatorHB(sim%delta_t)

       ! Diagnostics
       !TODO
       if (sim%rank == 0) then
          eenergy = sum(sim%propagator%efield_dofs(:,1)**2)*&
               sim%mesh%delta_eta
          write(th_diag_id,'(f12.5,2g20.12,2g20.12,2g20.12)' ) &
               real(j,f64)*sim%delta_t,  sim%propagator%efield_dofs(2,:), &
               sim%propagator%bfield(2)
       end if
    end do
    

  end subroutine run_pic_1d2v_vm

!------------------------------------------------------------------------------!

  subroutine delete_pic_1d2v_vm (sim)
    class(sll_sim_pic_1d2v_vm_cart), intent(inout) :: sim
  end subroutine delete_pic_1d2v_vm

!------------------------------------------------------------------------------!




end module sll_m_sim_pic_1d2v_vm_cart
