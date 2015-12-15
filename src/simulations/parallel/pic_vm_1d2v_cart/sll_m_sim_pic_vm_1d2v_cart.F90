! Simulation of 1d2v Vlasov-Maxwell with simple PIC method, periodic boundary conditions, Weibel instability. FEM with splines, degree 3 for B and 2 for E


module sll_m_sim_pic_vm_1d2v_cart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_ascii_io, only: &
    sll_ascii_file_create

  use sll_m_cartesian_meshes, only: &
    new_cartesian_mesh_1d, &
    sll_cartesian_mesh_1d

  use sll_m_collective, only: &
    sll_collective_allreduce, &
    sll_collective_reduce_real64, &
    sll_get_collective_rank, &
    sll_get_collective_size, &
    sll_world_collective

  use sll_m_constants, only: &
    sll_pi

  use sll_m_hamiltonian_splitting_base, only: &
    sll_t_hamiltonian_splitting_base

  use sll_m_hamiltonian_splitting_cef_pic_vm_1d2v, only: &
    sll_new_hamiltonian_splitting_cef_pic_vm_1d2v, &
    sll_t_hamiltonian_splitting_cef_pic_vm_1d2v

  use sll_m_hamiltonian_splitting_pic_vm_1d2v, only: &
    sll_new_hamiltonian_splitting_pic_vm_1d2v, &
    sll_t_hamiltonian_splitting_pic_vm_1d2v

  use sll_m_kernel_smoother_base, only: &
    sll_galerkin, &
    sll_c_kernel_smoother

  use sll_m_kernel_smoother_spline_1d, only: &
    sll_t_kernel_smoother_spline_1d, &
    sll_new_smoother_spline_1d

  use sll_m_maxwell_1d_base, only: &
    sll_maxwell_1d_base

  use sll_m_maxwell_1d_fem, only: &
    sll_maxwell_1d_fem, &
    sll_new_maxwell_1d_fem

  use sll_m_particle_group_1d2v, only: &
    sll_new_particle_group_1d2v, &
    sll_particle_group_1d2v

  use sll_m_particle_group_base, only: &
    sll_particle_group_base

  use sll_m_particle_initializer, only: &
    sll_particle_initialize_random_landau_1d2v, &
    sll_particle_initialize_sobol_landau_1d2v

  use sll_m_sim_base, only: &
    sll_simulation_base_class

  use sll_mpi, only: &
    mpi_sum

  implicit none

  public :: &
    sll_delete, &
    sll_t_sim_pic_vm_1d2v_cart

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: SLL_INIT_RANDOM=0
  sll_int32, parameter :: SLL_INIT_SOBOL=1

  sll_int32, parameter :: SLL_SPLITTING_SYMPLECTIC=0
  sll_int32, parameter :: SLL_SPLITTING_CEF=1

    type, extends(sll_simulation_base_class) :: sll_t_sim_pic_vm_1d2v_cart

     ! Abstract particle group
     class(sll_particle_group_base), pointer :: particle_group

     ! 
     sll_real64, pointer :: efield_dofs(:,:)
     sll_real64, pointer :: bfield_dofs(:)

     ! Cartesian mesh
     type(sll_cartesian_mesh_1d), pointer    :: mesh 

     ! Maxwell solver 
     ! Abstract 
     class(sll_maxwell_1d_base), pointer :: maxwell_solver

     ! Abstract kernel smoothers
     class(sll_c_kernel_smoother), pointer :: kernel_smoother_0     
     class(sll_c_kernel_smoother), pointer :: kernel_smoother_1


     ! Specific operator splitting
     class(sll_t_hamiltonian_splitting_base), pointer :: propagator
     sll_int32 :: splitting_case

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
     
     ! Case definitions
     sll_int32  :: init_case
     
   contains
     procedure :: init_from_file => init_pic_vm_1d2v
     procedure :: run => run_pic_vm_1d2v

  end type sll_t_sim_pic_vm_1d2v_cart

  interface sll_delete
     module procedure delete_pic_vm_1d2v
  end interface sll_delete
  
contains
!------------------------------------------------------------------------------!
  ! Read in the simulation parameters from input file
  subroutine init_pic_vm_1d2v (sim, filename)
    class(sll_t_sim_pic_vm_1d2v_cart), intent(inout) :: sim
    character(len=*), intent(in)                                :: filename

    sll_int32   :: io_stat
    sll_int32   :: n_time_steps
    sll_real64  :: delta_t, alpha, n_mode, thermal_v, T_r, beta
    sll_int32   :: ng_x
    sll_real64  :: x1_min, x1_max
    sll_int32   :: n_particles!, degree_smoother
    character(len=256)   :: init_case
    character(len=256)   :: splitting_case
    sll_int32   :: spline_degree 

    sll_int32, parameter :: input_file = 99
    

    namelist /sim_params/         delta_t, n_time_steps, alpha, n_mode, thermal_v, T_r, beta
    
    namelist /grid_dims/          ng_x, x1_min, x1_max

    namelist /pic_params/         n_particles, init_case, splitting_case, spline_degree!, degree_smoother

    ! Read parameters from file
    open(unit = input_file, file=trim(filename), IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'init_pic_1d2v() failed to open file ', filename
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
    sim%degree_smoother = spline_degree
    
    select case(init_case)
    case("SLL_INIT_RANDOM")
       sim%init_case = SLL_INIT_RANDOM
    case("SLL_INIT_SOBOL")
       sim%init_case = SLL_INIT_SOBOL
    case default
       print*, '#init case ', init_case, ' not implemented.'
    end select

    select case(splitting_case)
    case("SLL_SPLITTING_CEF")
       sim%splitting_case = SLL_SPLITTING_CEF
    case("SLL_SPLITTING_SYMPLECTIC")
       sim%splitting_case = SLL_SPLITTING_SYMPLECTIC
    case default
       print*, '#splitting case ', splitting_case, ' not implemented.'
    end select

  end subroutine init_pic_vm_1d2v

!------------------------------------------------------------------------------!

  subroutine run_pic_vm_1d2v (sim)
    class(sll_t_sim_pic_vm_1d2v_cart), intent(inout) :: sim

    ! Local variables
    sll_int32, allocatable :: rnd_seed(:)
    sll_int32 :: rnd_seed_size
    sll_int64 :: sobol_seed
    sll_int32 :: j, ierr, i_part
    sll_real64, allocatable :: rho(:), rho_local(:)
    sll_int32 :: th_diag_id

    ! For diagnostics
    sll_real64 :: kinetic_energy(1)
    sll_real64 :: total_energy(1)
    sll_real64 :: potential_energy(3)
    sll_real64 :: vi(3)
    sll_real64 :: wi(1)
    sll_real64 :: xi(3)
   
    ! For testing
    sll_int32  :: reffile_id 
    sll_real64 :: error

    ! Initialize file for diagnostics
    if (sim%rank == 0) then
       call sll_ascii_file_create('thdiag5.dat', th_diag_id, ierr)
    end if

    ! Initialize the particles   (mass and charge set to 1.0)
     sim%particle_group => sll_new_particle_group_1d2v(sim%n_particles, &
         sim%n_total_particles ,1.0_f64, 1.0_f64, 1)

    ! Initialize the field solver
    sim%maxwell_solver => sll_new_maxwell_1d_fem(sim%domain(1:2), sim%n_gcells, &
         sim%degree_smoother)

    ! Initialize kernel smoother    
    sim%kernel_smoother_1 => sll_new_smoother_spline_1d(&
         sim%domain(1:2), [sim%n_gcells], &
         sim%n_particles, sim%degree_smoother-1, SLL_GALERKIN) 
    sim%kernel_smoother_0 => &
         sll_new_smoother_spline_1d(sim%domain(1:2), [sim%n_gcells], &
         sim%n_particles, sim%degree_smoother, SLL_GALERKIN) 
    

    if (sim%init_case == SLL_INIT_RANDOM) then
       ! Set the seed for the random initialization
       call random_seed(size=rnd_seed_size)
       SLL_ALLOCATE(rnd_seed(rnd_seed_size), j)
       do j=1, rnd_seed_size
          rnd_seed(j) = (-1)**j*(100 + 15*j)*(2*sim%rank + 1)
       end do

       ! Initialize position and velocity of the particles.
       ! Random initialization
       call sll_particle_initialize_random_landau_1d2v &
            (sim%particle_group, sim%landau_param, &
            sim%domain(1) , &
            sim%domain(3), &
            sim%thermal_velocity, rnd_seed)
    elseif (sim%init_case == SLL_INIT_SOBOL) then
       sobol_seed = int(10 + sim%rank*sim%particle_group%n_particles, 8)
       ! Pseudorandom initialization with sobol numbers
       !sim%thermal_velocity = 0.1_f64
       call sll_particle_initialize_sobol_landau_1d2v(sim%particle_group, &
            sim%landau_param, sim%domain(1),sim%domain(3), &
            sim%thermal_velocity, sobol_seed)
    end if

    ! Initialize the arrays for the spline coefficients of the fields
    SLL_ALLOCATE(sim%efield_dofs(sim%n_gcells,2), ierr)
    SLL_ALLOCATE(sim%bfield_dofs(sim%n_gcells), ierr)

    ! Initialize the time-splitting propagator
    if (sim%splitting_case == SLL_SPLITTING_SYMPLECTIC) then
       sim%propagator => sll_new_hamiltonian_splitting_pic_vm_1d2v(sim%maxwell_solver, &
            sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
            sim%efield_dofs, sim%bfield_dofs, &
            sim%domain(1), sim%domain(3))
    elseif (sim%splitting_case == SLL_SPLITTING_CEF) then
       sim%propagator =>  sll_new_hamiltonian_splitting_cef_pic_vm_1d2v(sim%maxwell_solver, &
            sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
            sim%efield_dofs, sim%bfield_dofs, &
            sim%domain(1), sim%domain(3))
    end if

    ! Set the initial fields
    SLL_ALLOCATE(rho_local(sim%n_gcells), ierr)
    SLL_ALLOCATE(rho(sim%n_gcells), ierr)
    ! Efield 1 by Poisson
    rho_local = 0.0_f64
    do i_part = 1, sim%particle_group%n_particles
       xi = sim%particle_group%get_x(i_part)
       wi(1) = sim%particle_group%get_charge( i_part)
       call sim%kernel_smoother_0%add_charge(xi(1), wi(1), rho_local)
    end do
    ! MPI to sum up contributions from each processor
    rho = 0.0_f64
    call sll_collective_allreduce( sll_world_collective, &
         rho_local, &
         sim%n_gcells, MPI_SUM, rho)
    ! Solve Poisson problem
    call sim%maxwell_solver%compute_E_from_rho(sim%efield_dofs(:,1),&
         rho)

    ! Efield 2 to zero
    sim%efield_dofs(:,2) = 0.0_f64
    ! Bfield = beta*cos(kx): Use b = M{-1}(N_i,beta*cos(kx))
    call sim%maxwell_solver%L2projection( beta_cos_k, sim%degree_smoother-1, &
         sim%bfield_dofs) 
    !print*, sim%propagator%bfield_dofs

    ! End field initialization

    ! Allocate the vector holding the values of the fields at the grid points
    SLL_ALLOCATE(sim%fields_grid(sim%n_gcells,3), ierr)

    ! Diagnostics
    ! Kinetic energy
    kinetic_energy = 0.0_f64
    do i_part=1,sim%particle_group%n_particles
       vi = sim%particle_group%get_v(i_part)
       wi = sim%particle_group%get_mass(i_part)
       kinetic_energy = kinetic_energy + &
            (vi(1)**2+vi(2)**2)*wi(1)
    end do
    total_energy = 0.0_f64
    call sll_collective_reduce_real64(sll_world_collective, kinetic_energy, 1,&
         MPI_SUM, 0, total_energy)
    if (sim%rank == 0) then
       potential_energy(1) = sim%maxwell_solver%L2norm_squared&
            (sim%efield_dofs(:,1), sim%degree_smoother-1)
       potential_energy(2) = sim%maxwell_solver%L2norm_squared&
            (sim%efield_dofs(:,2), sim%degree_smoother)
       potential_energy(3) = sim%maxwell_solver%L2norm_squared&
            ( sim%bfield_dofs, sim%degree_smoother-1)
       write(th_diag_id,'(f12.5,2g20.12,2g20.12,2g20.12,2g20.12,2g20.12)' ) &
            0.0_f64,  &
            potential_energy, total_energy, total_energy + sum(potential_energy)
       print*, 'Time loop'
    end if
    ! Time loop
    do j=1, sim%n_time_steps
       ! Strang splitting
       call sim%propagator%strang_splitting(sim%delta_t,1)

       ! Diagnostics
       ! Kinetic energy
       kinetic_energy = 0.0_f64
       do i_part=1,sim%particle_group%n_particles
          vi = sim%particle_group%get_v(i_part)
          wi = sim%particle_group%get_mass(i_part)
          kinetic_energy = kinetic_energy + &
               (vi(1)**2+vi(2)**2)*wi(1)
       end do
       total_energy = 0.0_f64
       call sll_collective_reduce_real64(sll_world_collective, kinetic_energy, 1,&
            MPI_SUM, 0, total_energy)

       if (sim%rank == 0) then
          print*, 'Iteration=', j 
          potential_energy(1) = sim%maxwell_solver%L2norm_squared&
               (sim%efield_dofs(:,1), sim%degree_smoother-1)
          potential_energy(2) = sim%maxwell_solver%L2norm_squared&
               (sim%efield_dofs(:,2), sim%degree_smoother)
          potential_energy(3) = sim%maxwell_solver%L2norm_squared&
               ( sim%bfield_dofs, sim%degree_smoother-1)
          write(th_diag_id,'(f12.5,2g20.12,2g20.12,2g20.12,2g20.12,2g20.12)' ) &
               real(j,f64)*sim%delta_t,  &
               potential_energy, total_energy, total_energy + sum(potential_energy)
       end if
    end do
    
    ! Compute final rho
    rho_local = 0.0_f64
    do i_part = 1, sim%particle_group%n_particles
       xi = sim%particle_group%get_x(i_part)
       wi(1) = sim%particle_group%get_charge( i_part)
       call sim%kernel_smoother_0%add_charge(xi(1), wi(1), rho_local)
    end do
    ! MPI to sum up contributions from each processor
    rho = 0.0_f64
    call sll_collective_allreduce( sll_world_collective, &
         rho_local, &
         sim%n_gcells, MPI_SUM, rho)
    write(32,'(D16.25)') rho
    write(33,*) sim%efield_dofs(:,1)
    write(34,*) sim%efield_dofs(:,2)
    write(35,*) sim%bfield_dofs
    close(32)
    close(33)
    close(34)
    close(35)

    if (sim%rank == 0) then
       reffile_id = 21
       rho_local = 0.0_f64
       open(unit=reffile_id, &
            file='reffile_pic_vm_1d2v_cart.dat', &
            IOStat=ierr)
       read(reffile_id,*,IOSTAT=ierr) rho_local
       !call sll_binary_read_array_1d(reffile_id, rho_local, ierr)
       rho_local = rho_local -  rho
       error = maxval(rho_local)
       print*, error
       if (abs(error)> 1E-14) then
          print*, 'FAILED'
       else
          print*, 'PASSED'
       end if
    end if
    
  contains
    function beta_cos_k(x)
      sll_real64             :: beta_cos_k
      sll_real64, intent(in) :: x

      beta_cos_k = sim%beta * cos(2*sll_pi*x/sim%domain(3)) 
    end function beta_cos_k
  end subroutine run_pic_vm_1d2v

!------------------------------------------------------------------------------!

  subroutine delete_pic_vm_1d2v (sim)
    class(sll_t_sim_pic_vm_1d2v_cart), intent(inout) :: sim
    SLL_ASSERT(storage_size(sim)>0)
  end subroutine delete_pic_vm_1d2v

!------------------------------------------------------------------------------!




end module sll_m_sim_pic_vm_1d2v_cart
