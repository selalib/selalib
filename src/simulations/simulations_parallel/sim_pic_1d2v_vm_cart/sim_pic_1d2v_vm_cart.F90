! Simulation of 1d2v Vlasov-Maxwell with simple PIC method, periodic boundary conditions, Weibel instability. FEM with splines, degree 3 for B and 2 for E


module sll_m_sim_pic_1d2v_vm_cart

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_utilities.h"

  use sll_constants
  use sll_collective
  use sll_simulation_base
  use sll_cartesian_meshes
  use sll_module_pic_base
  use sll_m_particle_initializer
  use sll_m_particle_group_1d2v
  use sll_m_operator_splitting_pic_1d2v_vm
  use sll_ascii_io
  
  use sll_m_kernel_smoother_base
  use sll_m_kernel_smoother_spline_1d
  use sll_m_maxwell_1d_base
  use sll_m_maxwell_1d_fem
  use sll_arbitrary_degree_splines

    type, extends(sll_simulation_base_class) :: sll_sim_pic_1d2v_vm_cart

     ! Abstract particle group
     class(sll_particle_group_base), pointer :: particle_group
     ! Specific particle group
     class(sll_particle_group_1d2v), pointer :: specific_particle_group 

     ! Cartesian mesh
     type(sll_cartesian_mesh_1d), pointer    :: mesh 

     ! Maxwell solver 
     ! Abstract 
     class(sll_maxwell_1d_base), pointer :: maxwell_solver
     ! Specific
     class(sll_maxwell_1d_fem), pointer :: specific_maxwell_solver

     ! Abstract kernel smoothers
     class(sll_kernel_smoother_base), pointer :: kernel_smoother_0     
     class(sll_kernel_smoother_base), pointer :: kernel_smoother_1
     ! Specific kernel smoother
     class(sll_kernel_smoother_spline_1d), pointer :: specific_kernel_smoother_0
     class(sll_kernel_smoother_spline_1d), pointer :: specific_kernel_smoother_1


     ! Specific operator splitting
     class(sll_operator_splitting_pic_1d2v_vm), pointer :: propagator

     ! Fields on the grid
     sll_real64, allocatable :: fields_grid(:,:)
     
     ! Physical parameters
     sll_real64 :: landau_param(2) ! (1+landau_param(1)*cos(landau_param(2)*x1) 
     sll_real64 :: beta
     sll_real64 :: thermal_velocity(2) ! 
     sll_real64 :: domain(3) ! x_min, x_max, Lx
     
     ! Simulation parameters
     sll_real64 :: delta_t
     sll_int32  :: n_time_steps
     sll_int32  :: n_particles
     sll_int32  :: n_total_particles
     sll_int32  :: degree_smoother
     sll_int32  :: n_gcells

     
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
    sll_real64  :: delta_t, alpha, n_mode, thermal_v, T_r, beta
    sll_int32   :: ng_x
    sll_real64  :: x1_min, x1_max
    sll_int32   :: n_particles!, degree_smoother

    sll_int32, parameter :: input_file = 99
    

    namelist /sim_params/         delta_t, n_time_steps, alpha, n_mode, thermal_v, T_r, beta
    
    namelist /grid_dims/          ng_x, x1_min, x1_max

    namelist /pic_params/         n_particles!, degree_smoother

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

    ! Set MPI parameters
    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%rank = sll_get_collective_rank(sll_world_collective)

    ! Copy the read parameters into the simulation parameters
    sim%delta_t = delta_t
    sim%n_time_steps = n_time_steps
    sim%landau_param = [alpha, n_mode*2.0_f64 * sll_pi/(x1_max - x1_min)]
    sim%thermal_velocity = [thermal_v, thermal_v* sqrt(T_r)]
    sim%beta = beta

    sim%n_gcells = ng_x
    sim%mesh => new_cartesian_mesh_1d( ng_x, &
         x1_min, x1_max)
    sim%domain = [x1_min, x1_max, x1_max - x1_min ]

    sim%n_particles = n_particles/sim%world_size
    sim%n_total_particles = sim%n_particles * sim%world_size
    sim%degree_smoother = 3!degree_smoother
    

  end subroutine init_pic_1d2v_vm

!------------------------------------------------------------------------------!

  subroutine run_pic_1d2v_vm (sim)
    class(sll_sim_pic_1d2v_vm_cart), intent(inout) :: sim

    ! Local variables
    sll_int32, allocatable :: rnd_seed(:)
    sll_int32 :: j, ierr
    sll_int32 :: rnd_seed_size
    sll_real64 :: eenergy
    sll_int32 :: th_diag_id

    ! For diagnostics
    sll_real64 :: kinetic_energy(1)
    sll_real64 :: total_energy(1)
    sll_real64 :: potential_energy(3)
    sll_real64 :: vi(3)

    ! Initialize file for diagnostics
    if (sim%rank == 0) then
       call sll_ascii_file_create('thdiag.dat', th_diag_id, ierr)
    end if

    ! Initialize the particles   (mass and charge set to 1.0)
     sim%specific_particle_group => sll_new_particle_group_1d2v(sim%n_particles, &
         sim%n_total_particles ,1.0_f64, 1.0_f64)
    sim%particle_group => sim%specific_particle_group

    ! Initialize the field solver
    sim%specific_maxwell_solver => sll_new_maxwell_1d_fem(sim%domain(1:2), sim%n_gcells, &
         sim%degree_smoother)
    sim%maxwell_solver => sim%specific_maxwell_solver

    ! Initialize kernel smoother    
    sim%specific_kernel_smoother_1 => sll_new_smoother_spline_1d(&
         sim%domain, [sim%n_gcells], &
         sim%n_particles, sim%degree_smoother) 
    sim%kernel_smoother_1 => sim%specific_kernel_smoother_1
    sim%specific_kernel_smoother_0 => &
         sll_new_smoother_spline_1d(sim%domain(1:2), [sim%n_gcells], &
         sim%n_particles, sim%degree_smoother-1) 
    sim%kernel_smoother_0 => sim%specific_kernel_smoother_0
    
    ! Set the seed for the random initialization
    call random_seed(size=rnd_seed_size)
    SLL_ALLOCATE(rnd_seed(rnd_seed_size), j)
    do j=1, rnd_seed_size
       rnd_seed(j) = (-1)**j*(100 + 15*j)*(2*sim%rank + 1)
    end do

    ! Initialize position and velocity of the particles.
    call sll_particle_initialize_landau_1d2v (sim%particle_group, sim%landau_param, &
         sim%domain(1) , &
         sim%domain(3), &
         sim%thermal_velocity, rnd_seed)

    ! Initialize the time-splitting propagator
    sim%propagator => sll_new_splitting_pic_1d2v_vm(sim%maxwell_solver, &
         sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
         sim%domain(1), sim%domain(3))

    ! Set the initial fields
    ! Efield 1 by Poisson
    call sim%kernel_smoother_0%compute_shape_factors(sim%particle_group)
    sim%propagator%j_dofs_local = 0.0_f64
    call sim%kernel_smoother_0%accumulate_rho_from_klimontovich(sim%particle_group, &
         sim%propagator%j_dofs_local(:,1))
    ! MPI to sum up contributions from each processor
    sim%propagator%j_dofs = 0.0_f64
    call sll_collective_allreduce( sll_world_collective, &
         sim%propagator%j_dofs_local(:,1), &
         sim%n_gcells, MPI_SUM, sim%propagator%j_dofs(:,1))
    ! Solve Poisson problem
    call sim%maxwell_solver%compute_E_from_rho(sim%propagator%efield_dofs(:,1),&
         sim%propagator%j_dofs(:,1))

    ! Efield 2 to zero
    sim%propagator%efield_dofs(:,2) = 0.0_f64
    ! Bfield = beta*cos(kx): Use b = M{-1}(N_i,beta*cos(kx))
    call L2projection(sim%specific_maxwell_solver, beta_cos_k, sim%degree_smoother-1, &
         sim%propagator%bfield_dofs) 
    !print*, sim%propagator%bfield_dofs

    ! End field initialization

    ! Allocate the vector holding the values of the fields at the grid points
    SLL_ALLOCATE(sim%fields_grid(sim%n_gcells,3), ierr)

    ! Diagnostics
    ! Kinetic energy
    kinetic_energy = 0.0_f64
    do i_part=1,sim%particle_group%n_particles
       vi = sim%particle_group%get_v(i_part)
       kinetic_energy = kinetic_energy + &
            (vi(1)+vi(2))*sim%particle_group%get_charge(i_part)
    end do
    total_energy = 0.0_f64
    call sll_collective_reduce_real64(sll_world_collective, kinetic_energy, 1,&
         MPI_SUM, 0, total_energy)
    if (sim%rank == 0) then
       ! Compute fields at grid points from dofs
       sim%fields_grid(:,1) = eval_uniform_periodic_spline_curve(sim%degree_smoother-1, &
            sim%propagator%efield_dofs(:,1))
       sim%fields_grid(:,2) = eval_uniform_periodic_spline_curve(sim%degree_smoother, &
            sim%propagator%efield_dofs(:,2))
       sim%fields_grid(:,3) = eval_uniform_periodic_spline_curve(sim%degree_smoother-1, &
            sim%propagator%bfield_dofs) 
       potential_energy(1) = sum(sim%fields_grid(:,1)**2)*sim%mesh%delta_eta
       potential_energy(2) = sum(sim%fields_grid(:,2)**2)*sim%mesh%delta_eta
       potential_energy(3) = sum(sim%fields_grid(:,3)**2)*sim%mesh%delta_eta
       total_energy = total_energy + sum(potential_energy)
       write(th_diag_id,'(f12.5,2g20.12,2g20.12,2g20.12,2g20.12)' ) &
            0.0_f64,  &
            potential_energy, total_energy
       print*, 'Time loop'
    end if
    ! Time loop
    do j=1, sim%n_time_steps
       ! Lie splitting
!!$       call sim%propagator%operatorHf(sim%delta_t)
!!$       call sim%propagator%operatorHE(sim%delta_t)
!!$       call sim%propagator%operatorHB(sim%delta_t)
       ! Strang splitting
       call sim%propagator%operatorHf(0.5_f64*sim%delta_t)
       call sim%propagator%operatorHE(0.5_f64*sim%delta_t)
       call sim%propagator%operatorHB(sim%delta_t)
       call sim%propagator%operatorHE(0.5_f64*sim%delta_t)
       call sim%propagator%operatorHf(0.5_f64*sim%delta_t)

       ! Diagnostics

       ! Kinetic energy
       kinetic_energy = 0.0_f64
       do i_part=1,sim%particle_group%n_particles
          vi = sim%particle_group%get_v(i_part)
          kinetic_energy = kinetic_energy + &
               (vi(1)+vi(2))*sim%particle_group%get_charge(i_part)
       end do
       total_energy = 0.0_f64
       call sll_collective_reduce_real64(sll_world_collective, kinetic_energy, 1,&
            MPI_SUM, 0, total_energy)

       if (sim%rank == 0) then
          print*, 'Iteration=', j 
          ! Compute fields at grid points from dofs
          sim%fields_grid(:,1) = eval_uniform_periodic_spline_curve(&
               sim%degree_smoother-1, sim%propagator%efield_dofs(:,1))
          sim%fields_grid(:,2) = eval_uniform_periodic_spline_curve(&
               sim%degree_smoother, sim%propagator%efield_dofs(:,2))
          sim%fields_grid(:,3) = eval_uniform_periodic_spline_curve(&
               sim%degree_smoother-1, sim%propagator%bfield_dofs)
          potential_energy(1) = sum(sim%fields_grid(:,1)**2)*sim%mesh%delta_eta
          potential_energy(2) = sum(sim%fields_grid(:,2)**2)*sim%mesh%delta_eta
          potential_energy(3) = sum(sim%fields_grid(:,3)**2)*sim%mesh%delta_eta
          total_energy = total_energy + sum(potential_energy)
          write(th_diag_id,'(f12.5,2g20.12,2g20.12,2g20.12,2g20.12)' ) &
               real(j,f64)*sim%delta_t,  &
               potential_energy, total_energy
       end if
    end do
    
  contains
    function beta_cos_k(x)
      sll_real64             :: beta_cos_k
      sll_real64, intent(in) :: x

      beta_cos_k = sim%beta * cos(2*sll_pi*x/sim%domain(3)) 
    end function beta_cos_k
  end subroutine run_pic_1d2v_vm

!------------------------------------------------------------------------------!

  subroutine delete_pic_1d2v_vm (sim)
    class(sll_sim_pic_1d2v_vm_cart), intent(inout) :: sim
  end subroutine delete_pic_1d2v_vm

!------------------------------------------------------------------------------!




end module sll_m_sim_pic_1d2v_vm_cart
