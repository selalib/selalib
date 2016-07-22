! Simulation of 1d2v Vlasov-Maxwell with simple PIC method, periodic boundary conditions, Weibel instability. FEM with splines, degree 3 for B and 2 for E

! author: Katharina Kormann, IPP

module sll_m_sim_pic_vm_1d2v_cart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_ascii_io, only: &
    sll_s_ascii_file_create

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
    sll_p_pi

  use sll_m_hamiltonian_splitting_base, only: &
    sll_c_hamiltonian_splitting_base

  use sll_m_hamiltonian_splitting_pic_vm_1d2v, only: &
    sll_s_new_hamiltonian_splitting_pic_vm_1d2v, &
    sll_t_hamiltonian_splitting_pic_vm_1d2v
  
  use sll_m_hamiltonian_splitting_pic_vm_1d2v_boris, only: &
       sll_t_hamiltonian_splitting_pic_vm_1d2v_boris

  use sll_m_io_utilities, only : &
    sll_s_read_data_real_array, &
    sll_s_concatenate_filename_and_path

  use sll_m_particle_mesh_coupling_base, only: &
    sll_p_galerkin, &
    sll_c_particle_mesh_coupling

  use sll_m_particle_mesh_coupling_spline_1d, only: &
    sll_t_particle_mesh_coupling_spline_1d, &
    sll_s_new_particle_mesh_coupling_spline_1d_ptr

  use sll_m_maxwell_1d_base, only: &
    sll_c_maxwell_1d_base

  use sll_m_maxwell_1d_fem, only: &
    sll_t_maxwell_1d_fem

  use sll_m_particle_group_1d2v, only: &
    sll_s_new_particle_group_1d2v_ptr, &
    sll_t_particle_group_1d2v

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  use sll_m_particle_initializer, only: &
    sll_s_particle_initialize_random_landau_1d2v, &
    sll_s_particle_initialize_sobol_landau_1d2v, &
    sll_s_particle_initialize_random_landau_symmetric_1d2v, &
    sll_s_particle_initialize_sobol_landau_symmetric_1d2v, &
    sll_s_particle_initialize_sobol_streaming_weibel_symmetric_1d2v, &
    sll_s_particle_initialize_twogaussian_sobol_symmetric_1d2v

  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_t_time_mark

  use sll_mpi, only: &
       mpi_sum

  use  sll_m_pic_diagnostics, only : &
       sll_s_pic_diagnostics_transfer, &
       sll_s_pic_diagnostics_vvb, &
       sll_s_pic_diagnostics_poynting

  implicit none

  public :: &
    sll_t_sim_pic_vm_1d2v_cart

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: sll_p_init_random=0
  sll_int32, parameter :: sll_p_init_sobol=1
  sll_int32, parameter :: sll_p_init_random_sym=2
  sll_int32, parameter :: sll_p_init_sobol_sym=3

  sll_int32, parameter :: sll_p_splitting_symplectic = 0
  sll_int32, parameter :: sll_p_splitting_boris = 1

  sll_int32, parameter :: sll_p_onegaussian = 0
  sll_int32, parameter :: sll_p_twogaussian = 1

  sll_int32, parameter :: sll_p_bfield_cos = 0
  sll_int32, parameter :: sll_p_bfield_sin = 1
  
    type, extends(sll_c_simulation_base_class) :: sll_t_sim_pic_vm_1d2v_cart

     ! Abstract particle group
     class(sll_c_particle_group_base), pointer :: particle_group

     ! 
     sll_real64, pointer :: efield_dofs(:,:)
     sll_real64, pointer :: efield_dofs_n(:,:)
     sll_real64, pointer :: bfield_dofs(:)

     sll_real64, allocatable :: x_array(:)
     sll_real64, allocatable :: field_grid(:)

     ! Cartesian mesh
     type(sll_t_cartesian_mesh_1d), pointer    :: mesh 

     ! Maxwell solver 
     ! Abstract 
     class(sll_c_maxwell_1d_base), pointer :: maxwell_solver

     ! Abstract kernel smoothers
     class(sll_c_particle_mesh_coupling), pointer :: kernel_smoother_0     
     class(sll_c_particle_mesh_coupling), pointer :: kernel_smoother_1


     ! Specific operator splitting
     class(sll_c_hamiltonian_splitting_base), allocatable :: propagator
     sll_int32 :: splitting_case

     ! Fields on the grid
     sll_real64, allocatable :: fields_grid(:,:)
     
     ! Physical parameters
     sll_real64 :: landau_param(2) ! (1+landau_param(1)*cos(landau_param(2)*x1)
     sll_real64 :: velocity_params(2,5)
     
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
     sll_int32  :: test_case
     sll_int32  :: initial_bfield
     
     ! For ctest
     logical    :: ctest_passed
     
   contains
     procedure :: init_from_file => init_pic_vm_1d2v
     procedure :: run => run_pic_vm_1d2v
     procedure :: delete => delete_pic_vm_1d2v

  end type sll_t_sim_pic_vm_1d2v_cart

  
contains
!------------------------------------------------------------------------------!
  ! Read in the simulation parameters from input file
  subroutine init_pic_vm_1d2v (sim, filename)
    class(sll_t_sim_pic_vm_1d2v_cart), intent(inout) :: sim
    character(len=*),                  intent(in)    :: filename

    sll_int32   :: io_stat
    sll_int32   :: n_time_steps
    sll_real64  :: delta_t, alpha, n_mode, beta
    sll_real64  :: sig1, sig2
    sll_real64  :: portion1, mu11, sig11, mu12, sig12
    sll_real64  :: portion2, mu21, sig21, mu22, sig22
    sll_int32   :: ng_x
    sll_real64  :: x1_min, x1_max
    sll_int32   :: n_particles!, degree_smoother
    character(len=256)   :: init_case
    character(len=256)   :: splitting_case
    character(len=256)   :: test_case
    character(len=256)   :: initial_bfield
    sll_int32   :: spline_degree 

    sll_int32   :: input_file
    sll_int32   :: ierr, j
    

    namelist /sim_params/         delta_t, n_time_steps, alpha, n_mode, beta, test_case, initial_bfield
    namelist /twogaussian_params/ portion1, mu11, sig11, mu12, sig12, portion2, mu21, sig21, mu22, sig22
    namelist /onegaussian_params/ sig1, sig2
    
    namelist /grid_dims/          ng_x, x1_min, x1_max

    namelist /pic_params/         n_particles, init_case, splitting_case, spline_degree!, degree_smoother

    ! Read parameters from file
    open(newunit = input_file, file=trim(filename), IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'init_pic_1d2v() failed to open file ', filename
       STOP
    end if

    read(input_file, sim_params)

    select case (test_case )
    case( "onegaussian" )
       sim%test_case = sll_p_onegaussian
       read(input_file, onegaussian_params)
       sim%thermal_velocity = [ sig1, sig2 ]
    case( "twogaussian")
       sim%test_case = sll_p_twogaussian
       read(input_file, twogaussian_params)
       SLL_ASSERT( portion1 <= 1.0_f64 )
       SLL_ASSERT( portion1 >= 0.0_f64 )
       SLL_ASSERT( portion2 <= 1.0_f64 )
       SLL_ASSERT( portion2 >= 0.0_f64 )
       
       sim%velocity_params(1,:) = [portion1, mu11, sig11, mu12, sig12]
       sim%velocity_params(2,:) = [portion2, mu21, sig21, mu22, sig22]
    case default
       print*, '#test case ', test_case, ' not implemented.'
    end select
       
    read(input_file, grid_dims)
    read(input_file, pic_params)
    close (input_file)

    ! Set MPI parameters
    sim%world_size = sll_f_get_collective_size(sll_v_world_collective)
    sim%rank = sll_f_get_collective_rank(sll_v_world_collective)

    ! Copy the read parameters into the simulation parameters
    sim%delta_t = delta_t
    sim%n_time_steps = n_time_steps
    sim%landau_param = [alpha, n_mode*2.0_f64 * sll_p_pi/(x1_max - x1_min)]

    sim%beta = beta

    select case ( initial_bfield )
    case( "cos")
       sim%initial_bfield = sll_p_bfield_cos
    case( "sin" )
       sim%initial_bfield = sll_p_bfield_sin
    case default
       print*, '#initial bfield must be either sin or cos.'
    end select

    
    sim%n_gcells = ng_x
    sim%mesh => sll_f_new_cartesian_mesh_1d( ng_x, &
         x1_min, x1_max)
    sim%domain = [x1_min, x1_max, x1_max - x1_min ]

    sim%n_particles = n_particles/sim%world_size
    sim%n_total_particles = sim%n_particles * sim%world_size
    sim%degree_smoother = spline_degree
    
    select case(init_case)
    case("SLL_INIT_RANDOM")
       sim%init_case = sll_p_init_random
    case("SLL_INIT_SOBOL")
       sim%init_case = sll_p_init_sobol
    case("SLL_INIT_RANDOM_SYM")
       sim%init_case = sll_p_init_random_sym
    case("SLL_INIT_SOBOL_SYM")
       sim%init_case = sll_p_init_sobol_sym
    case default
       print*, '#init case ', init_case, ' not implemented.'
    end select

    select case(splitting_case)
    case("SLL_SPLITTING_SYMPLECTIC")
       sim%splitting_case = sll_p_splitting_symplectic
    case("SLL_SPLITTING_BORIS")
       sim%splitting_case = sll_p_splitting_boris
    case default
       print*, '#splitting case ', splitting_case, ' not implemented.'
    end select

    ! Initialize the particles   (mass and charge set to 1.0)
     call sll_s_new_particle_group_1d2v_ptr(sim%particle_group, sim%n_particles, &
         sim%n_total_particles, 1.0_f64, 1.0_f64, 1)

    ! Initialize the field solver
     allocate( sll_t_maxwell_1d_fem :: sim%maxwell_solver )
     select type ( q=>sim%maxwell_solver )
     type is ( sll_t_maxwell_1d_fem )
        call q%init( sim%domain(1:2), sim%n_gcells, &
             sim%degree_smoother)
     end select

    ! Initialize kernel smoother    
    call sll_s_new_particle_mesh_coupling_spline_1d_ptr(sim%kernel_smoother_1, &
         sim%domain(1:2), [sim%n_gcells], &
         sim%n_particles, sim%degree_smoother-1, sll_p_galerkin) 
    call sll_s_new_particle_mesh_coupling_spline_1d_ptr(sim%kernel_smoother_0, &
         sim%domain(1:2), [sim%n_gcells], &
         sim%n_particles, sim%degree_smoother, sll_p_galerkin) 
   

    ! Initialize the arrays for the spline coefficients of the fields
    SLL_ALLOCATE(sim%efield_dofs(sim%n_gcells,2), ierr)
    SLL_ALLOCATE(sim%bfield_dofs(sim%n_gcells), ierr)

    ! Initialize the time-splitting propagator
    if (sim%splitting_case == sll_p_splitting_symplectic) then
       call sll_s_new_hamiltonian_splitting_pic_vm_1d2v(&
            sim%propagator, sim%maxwell_solver, &
            sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
            sim%efield_dofs, sim%bfield_dofs, &
            sim%domain(1), sim%domain(3))
       sim%efield_dofs_n => sim%efield_dofs
    elseif( sim%splitting_case == sll_p_splitting_boris) then
       allocate( sll_t_hamiltonian_splitting_pic_vm_1d2v_boris :: sim%propagator )
       select type( qp=>sim%propagator )
       type is ( sll_t_hamiltonian_splitting_pic_vm_1d2v_boris)
          call qp%init( sim%maxwell_solver, &
               sim%kernel_smoother_0, sim%kernel_smoother_1, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(1), sim%domain(3))
          sim%efield_dofs_n => qp%efield_dofs_mid
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
    sll_int32 :: th_diag_id, dfield_id, efield_id, bfield_id

    sll_real64 :: wi(1)
    sll_real64 :: xi(3)

    type(sll_t_time_mark) :: start_loop, end_loop
 
    ! Initialize file for diagnostics
    if (sim%rank == 0) then
       call sll_s_ascii_file_create('thdiag5.dat', th_diag_id, ierr)
       call sll_s_ascii_file_create('dfield.dat', dfield_id, ierr)
       call sll_s_ascii_file_create('efield.dat', efield_id, ierr)
       call sll_s_ascii_file_create('bfield.dat', bfield_id, ierr)
    end if

    if ( sim%test_case == sll_p_onegaussian ) then
       if (sim%init_case == sll_p_init_random) then
          ! Set the seed for the random initialization
          call random_seed(size=rnd_seed_size)
          SLL_ALLOCATE(rnd_seed(rnd_seed_size), j)
          do j=1, rnd_seed_size
             rnd_seed(j) = (-1)**j*(100 + 15*j)*(2*sim%rank + 1)
          end do
          
          ! Initialize position and velocity of the particles.
          ! Random initialization
          call sll_s_particle_initialize_random_landau_1d2v &
               (sim%particle_group, sim%landau_param, &
               sim%domain(1) , &
               sim%domain(3), &
               sim%thermal_velocity, rnd_seed)
       elseif (sim%init_case == sll_p_init_random_sym) then
          ! Set the seed for the random initialization
          call random_seed(size=rnd_seed_size)
          SLL_ALLOCATE(rnd_seed(rnd_seed_size), j)
          do j=1, rnd_seed_size
             rnd_seed(j) = (-1)**j*(100 + 15*j)*(2*sim%rank + 1)
          end do
          
          ! Initialize position and velocity of the particles.
          ! Random initialization
          call sll_s_particle_initialize_random_landau_symmetric_1d2v &
               (sim%particle_group, sim%landau_param, &
               sim%domain(1) , &
               sim%domain(3), &
               sim%thermal_velocity, rnd_seed)
       elseif (sim%init_case == sll_p_init_sobol) then
          sobol_seed = int(10 + sim%rank*sim%particle_group%n_particles/8, 8)
          ! Pseudorandom initialization with sobol numbers
          call sll_s_particle_initialize_sobol_landau_1d2v(sim%particle_group, &
               sim%landau_param, sim%domain(1),sim%domain(3), &
               sim%thermal_velocity, sobol_seed)
       elseif (sim%init_case == sll_p_init_sobol_sym) then
          sobol_seed = int(10 + sim%rank*sim%particle_group%n_particles/8, 8)
          ! Pseudorandom initialization with sobol numbers
          call sll_s_particle_initialize_sobol_landau_symmetric_1d2v(sim%particle_group, &
               sim%landau_param, sim%domain(1),sim%domain(3), &
               sim%thermal_velocity, sobol_seed)
       end if
    else
       if ( sim%init_case == sll_p_init_sobol_sym ) then
          sobol_seed = int(10 + sim%rank*sim%particle_group%n_particles/8, 8)
          call sll_s_particle_initialize_twogaussian_sobol_symmetric_1d2v(&
               sim%particle_group, &
               sim%landau_param, sim%domain(1),sim%domain(3), &
               sim%velocity_params, sobol_seed)
       else
          SLL_ERROR("run_pic_vm_1d2v_cart","init_case not implemented for twogaussian case.")
       end if
    end if

    
    ! Set the initial fields
    SLL_ALLOCATE(rho_local(sim%n_gcells), ierr)
    SLL_ALLOCATE(rho(sim%n_gcells), ierr)

    ! Efield 1 by Poisson
    call solve_poisson( sim%particle_group, sim%kernel_smoother_0, sim%maxwell_solver, rho_local, rho, sim%efield_dofs(:,1) )

    ! Efield 2 to zero
    sim%efield_dofs(:,2) = 0.0_f64
    ! Bfield = beta*cos(kx): Use b = M{-1}(N_i,beta*cos(kx))
    if ( sim%initial_bfield == sll_p_bfield_cos ) then
       call sim%maxwell_solver%L2projection( beta_cos_k, sim%degree_smoother-1, &
            sim%bfield_dofs)
    else
       call sim%maxwell_solver%L2projection( beta_sin_k, sim%degree_smoother-1, &
            sim%bfield_dofs)
    end if
       

    ! In case we use the Boris propagator, we need to initialize the staggering used in the scheme.
    select type( qp=>sim%propagator )
    type is ( sll_t_hamiltonian_splitting_pic_vm_1d2v_boris)
       call qp%staggering( sim%delta_t )
    end select

    ! End field initialization

    call solve_poisson( sim%particle_group, sim%kernel_smoother_0, sim%maxwell_solver, rho_local, rho, rho_local )
    ! Diagnostics
    call sll_s_time_history_diagnostics_pic_vm_1d2v( &
         sim%particle_group, sim%maxwell_solver, &
         sim%kernel_smoother_0, sim%kernel_smoother_1, 0.0_f64, &
         sim%degree_smoother, sim%efield_dofs, sim%bfield_dofs, &
         sim%rank, th_diag_id, sim%efield_dofs_n, rho_local, rho)

    if (sim%rank == 0 ) then
       call sll_s_set_time_mark(start_loop )
    end if

    
    ! Time loop
    do j=1, sim%n_time_steps
       !print*, 'TIME STEP', j
       ! Strang splitting
       call sim%propagator%strang_splitting(sim%delta_t,1)

       ! Diagnostics
       call solve_poisson( sim%particle_group, sim%kernel_smoother_0, sim%maxwell_solver, rho_local, rho, rho_local )
       call sll_s_time_history_diagnostics_pic_vm_1d2v( &
         sim%particle_group, sim%maxwell_solver, &
         sim%kernel_smoother_0, sim%kernel_smoother_1,  sim%delta_t*real(j,f64), &
         sim%degree_smoother, sim%efield_dofs, sim%bfield_dofs, &
         sim%rank, th_diag_id, sim%efield_dofs_n, rho_local, rho)

    end do

    if (sim%rank == 0 ) then
       call sll_s_set_time_mark( end_loop )
       write(*, "(A, F10.3)") "Main loop run time [s] = ", sll_f_time_elapsed_between( start_loop, end_loop)
    end if
    
    !!! Part for ctest
    ! Compute final rho
    rho_local = 0.0_f64
    do i_part = 1, sim%particle_group%n_particles
       xi = sim%particle_group%get_x(i_part)
       wi(1) = sim%particle_group%get_charge( i_part)
       call sim%kernel_smoother_0%add_charge(xi(1), wi(1), rho_local)
    end do
    ! MPI to sum up contributions from each processor
    rho = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, &
         rho_local, &
         sim%n_gcells, MPI_SUM, rho)

    if (sim%rank == 0) then

       call ctest( rho, rho_local, sim%ctest_passed )

    end if
    !!! Part for ctest end


  contains
    function beta_cos_k(x)
      sll_real64             :: beta_cos_k
      sll_real64, intent(in) :: x

      beta_cos_k = sim%beta * cos(2*sll_p_pi*x/sim%domain(3)) 
    end function beta_cos_k

    function beta_sin_k(x)
      sll_real64             :: beta_sin_k
      sll_real64, intent(in) :: x

      beta_sin_k = sim%beta * sin(2*sll_p_pi*x/sim%domain(3)) 
    end function beta_sin_k
    
  end subroutine run_pic_vm_1d2v

!------------------------------------------------------------------------------!
  ! local subroutine to handle ctest
  subroutine ctest(rho_simulated, rho_ref, passed)
    sll_real64, intent(in   ) :: rho_simulated(:)
    sll_real64, intent(inout) :: rho_ref(:)
    logical,    intent(  out) :: passed     

    ! For testing
    character(len=256) :: reffile
    sll_real64 :: error

    call sll_s_concatenate_filename_and_path( "reffile_pic_vm_1d2v_cart.dat", __FILE__,&
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
    class(sll_t_sim_pic_vm_1d2v_cart), intent(inout) :: sim
    SLL_ASSERT(storage_size(sim)>0)

    call sim%propagator%free()
    deallocate(sim%propagator)
    call sim%particle_group%free()
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

  end subroutine delete_pic_vm_1d2v

!------------------------------------------------------------------------------!

!> Diagnostics for PIC Vlasov-Maxwell 1d2v 
!> @todo (should be part of the library)
  subroutine sll_s_time_history_diagnostics_pic_vm_1d2v(&
       particle_group, &
       maxwell_solver, &
       kernel_smoother_0, &
       kernel_smoother_1, &
       time, &
       degree, &
       efield_dofs, &
       bfield_dofs, &
       mpi_rank, &
       file_id, &
       efield_dofs_n, &
       efield_poisson, scratch)
    class(sll_c_particle_group_base), intent(in) :: particle_group
    class(sll_c_maxwell_1d_base),     intent(in) :: maxwell_solver
    class(sll_c_particle_mesh_coupling),     intent(inout) :: kernel_smoother_0, kernel_smoother_1
    sll_real64,                       intent(in) :: time
    sll_real64,                       intent(in) :: efield_dofs(:,:)
    sll_real64,                       intent(in) :: efield_dofs_n(:,:)
    sll_real64,                       intent(in) :: efield_poisson(:)
    sll_real64,                       intent(out) :: scratch(:)
    sll_real64,                       intent(in) :: bfield_dofs(:)
    sll_int32,                        intent(in) :: degree
    sll_int32,                        intent(in) :: mpi_rank
    sll_int32,                        intent(in) :: file_id

    ! local variables
    sll_real64 :: diagnostics_local(3)
    sll_real64 :: diagnostics(3)
    sll_real64 :: potential_energy(3)
    sll_int32  :: i_part
    sll_real64 :: vi(3)
    sll_real64 :: wi(1)
    sll_real64 :: transfer(1), vvb(1), poynting

    diagnostics_local = 0.0_f64
    do i_part=1,particle_group%n_particles
       vi = particle_group%get_v(i_part)
       wi = particle_group%get_mass(i_part)
       ! Kinetic energy
       diagnostics_local(1) = diagnostics_local(1) + &
            (vi(1)**2+vi(2)**2)*wi(1)
       ! Momentum 1
       diagnostics_local(2) = diagnostics_local(2) + &
            vi(1)*wi(1)
       ! Momentum 2
       diagnostics_local(3) = diagnostics_local(3) + &
            vi(2)*wi(1)
    end do
    diagnostics = 0.0_f64
    call sll_s_collective_reduce_real64(sll_v_world_collective, diagnostics_local, 3,&
         MPI_SUM, 0, diagnostics)
    call sll_s_pic_diagnostics_transfer( particle_group, kernel_smoother_0, kernel_smoother_1, &
            efield_dofs, transfer )
    call sll_s_pic_diagnostics_vvb( particle_group, kernel_smoother_1, &
            bfield_dofs, vvb )
    call sll_s_pic_diagnostics_poynting( maxwell_solver, degree, efield_dofs(:,2), bfield_dofs, &
         scratch, poynting )
    
    
    
    if (mpi_rank == 0) then
       potential_energy(1) = maxwell_solver%inner_product&
            ( efield_dofs(:,1), efield_dofs_n(:,1), degree-1 )
       potential_energy(2) = maxwell_solver%inner_product&
            ( efield_dofs(:,2), efield_dofs_n(:,2), degree )
       potential_energy(3) = maxwell_solver%L2norm_squared&
            ( bfield_dofs, degree-1 )
       write(file_id,'(f12.5,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16)' ) &
            time,  potential_energy, diagnostics(1), &
            diagnostics(1) + sum(potential_energy), diagnostics(2:3), -transfer+vvb+poynting, &
            maxval(abs(efield_dofs(:,1)-efield_poisson))
    end if

  end subroutine sll_s_time_history_diagnostics_pic_vm_1d2v


  subroutine sll_s_diagnostics_fields( field_dofs,  field_grid, xi, n_dofs, kernel_smoother, file_id )
    sll_real64,                       intent(in) :: field_dofs(:)
    sll_real64,                       intent(inout) :: field_grid(:)
    sll_real64,                       intent(in) :: xi(:)
    sll_int32, intent(in) :: n_dofs
    class(sll_c_particle_mesh_coupling),     intent(inout) :: kernel_smoother
    sll_int32, intent(in) :: file_id
    

    sll_int32 :: j

    do j=1,n_dofs
       call kernel_smoother%evaluate( [xi(j)], field_dofs, field_grid(j))
    end do
    
    write(file_id, *) field_grid


  end subroutine sll_s_diagnostics_fields

  subroutine solve_poisson( particle_group, kernel_smoother_0, maxwell_solver, rho_local, rho, efield_dofs )
    class(sll_c_particle_group_base), intent(in) :: particle_group
    class(sll_c_maxwell_1d_base),     intent(in) :: maxwell_solver
    class(sll_c_particle_mesh_coupling),     intent(inout) :: kernel_smoother_0
    sll_real64,                       intent(inout) :: rho_local(:)
    sll_real64,                       intent(inout) :: rho(:)
    sll_real64,                       intent(inout) :: efield_dofs(:)
    
    
    sll_int32 :: i_part
    sll_real64 :: xi(3), wi(1)
    
    rho_local = 0.0_f64
    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x(i_part)
       wi(1) = particle_group%get_charge( i_part)
       call kernel_smoother_0%add_charge(xi(1), wi(1), rho_local)
    end do
    ! MPI to sum up contributions from each processor
    rho = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, &
         rho_local, &
         kernel_smoother_0%n_dofs, MPI_SUM, rho)
    ! Solve Poisson problem
    call maxwell_solver%compute_E_from_rho( efield_dofs,&
         rho )

  end subroutine solve_poisson
    

end module sll_m_sim_pic_vm_1d2v_cart
