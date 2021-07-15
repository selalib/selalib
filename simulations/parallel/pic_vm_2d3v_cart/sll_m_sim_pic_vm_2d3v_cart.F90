! Simulation of 1d2v Vlasov-Maxwell with simple PIC method, periodic boundary conditions, Weibel instability. FEM with splines, degree 3 for B and 2 for E

! author: Katharina Kormann, IPP

module sll_m_sim_pic_vm_2d3v_cart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_ascii_io, only: &
    sll_s_ascii_file_create

  use sll_m_collective, only: &
    sll_o_collective_allreduce, &
    sll_s_collective_reduce_real64, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_control_variate, only: &
    sll_t_control_variate, &
    sll_t_control_variates

  use sll_m_time_propagator_base, only: &
    sll_c_time_propagator_base

  use sll_m_time_propagator_pic_vm_2d3v_hs, only: &
    sll_s_new_time_propagator_pic_vm_2d3v_hs, &
    sll_t_time_propagator_pic_vm_2d3v_hs

  use sll_m_initial_distribution, only : &
       sll_c_distribution_params, &
       sll_s_initial_distribution_new
  
  use sll_m_io_utilities, only : &
    sll_s_read_data_real_array, &
    sll_s_concatenate_filename_and_path

  use sll_m_particle_mesh_coupling_spline_2d_feec, only: &
    sll_t_particle_mesh_coupling_spline_2d_feec

  use sll_m_maxwell_2d_fem_fft, only: &
       sll_t_maxwell_2d_fem_fft

  use sll_m_particle_group_2d3v, only: &
       sll_t_particle_group_2d3v, &
       sll_s_new_particle_group_2d3v_ptr

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base, &
       sll_t_particle_array

  use sll_m_particle_sampling, only: &
       sll_t_particle_sampling
  
  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_t_time_mark

  use sll_mpi, only: &
       mpi_sum

  implicit none

  public :: &
    sll_t_sim_pic_vm_2d3v_cart

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: sll_p_splitting_hs = 0
  sll_int32, parameter :: sll_p_splitting_boris = 1
  sll_int32, parameter :: sll_p_splitting_disgradE = 2
  sll_int32, parameter :: sll_p_splitting_nldisgradE = 3
  sll_int32, parameter :: sll_p_splitting_cef = 4
  sll_int32, parameter :: sll_p_splitting_disgradE_trunc = 5
  sll_int32, parameter :: sll_p_splitting_Hpsplit = 6


  sll_int32, parameter :: sll_p_onegaussian = 0
  sll_int32, parameter :: sll_p_twogaussian = 1

  sll_int32, parameter :: sll_p_bfield_cos = 0
  sll_int32, parameter :: sll_p_bfield_sin = 1
  sll_int32, parameter :: sll_p_bfield_constant = 2
  
    type, extends(sll_c_simulation_base_class) :: sll_t_sim_pic_vm_2d3v_cart

     ! Abstract particle group
     class(sll_t_particle_array), pointer :: particle_group

     ! Control variate
     type(sll_t_control_variates), pointer :: control_variate
     sll_int32 :: no_weights
     
     ! 
     sll_real64, pointer :: efield_dofs(:)
     sll_real64, pointer :: bfield_dofs(:)

     !sll_real64, allocatable :: x_array(:)
     !sll_real64, allocatable :: field_grid(:)

     ! Cartesian mesh
     sll_real64 :: delta_x(2)
     sll_real64 :: volume

     ! Maxwell solver 
     ! Abstract 
     type(sll_t_maxwell_2d_fem_fft), pointer :: maxwell_solver

     ! Abstract kernel smoothers
     type(sll_t_particle_mesh_coupling_spline_2d_feec), pointer :: particle_mesh_coupling


     ! Specific operator splitting
     class(sll_c_time_propagator_base), allocatable :: propagator
     sll_int32 :: splitting_case

     ! Fields on the grid
     sll_real64, allocatable :: fields_grid(:,:)
     
     ! Physical parameters
     class(sll_c_distribution_params), allocatable :: init_distrib_params
     sll_real64 :: beta
     sll_real64 :: domain(2,3) ! x_min, x_max, Lx
     type(sll_t_particle_sampling) :: sampler
     sll_real64 :: plasma_beta = 1.0_f64
     
     
     ! Simulation parameters
     sll_real64 :: delta_t
     sll_int32  :: n_time_steps
     sll_int32  :: n_particles
     sll_int32  :: n_total_particles
     sll_int32  :: degree_smoother
     sll_int32  :: n_gcells(2)
     sll_int32  :: n_totaldofs

     
     ! Parameters for MPI
     sll_int32  :: rank
     sll_int32  :: world_size
     
     ! Case definitions
     sll_int32  :: initial_bfield
     
     ! For ctest
     logical    :: ctest_passed
     logical    :: make_ctest = .false.
     character(len=256)   :: ctest_ref_file

     ! Output
     character(len=256)   :: file_prefix
     
     
     
   contains
     procedure :: init_from_file => init_pic_vm_2d3v
     procedure :: run => run_pic_vm_2d3v
     procedure :: delete => delete_pic_vm_2d3v

  end type sll_t_sim_pic_vm_2d3v_cart

  
contains
!------------------------------------------------------------------------------!
  ! Read in the simulation parameters from input file
  subroutine init_pic_vm_2d3v (sim, filename)
    class(sll_t_sim_pic_vm_2d3v_cart), intent(inout) :: sim
    character(len=*),                  intent(in)    :: filename

    sll_int32   :: io_stat
    sll_int32   :: n_time_steps
    sll_real64  :: delta_t,  beta
    sll_int32   :: ng_x(2)
    sll_real64  :: x_min(2), x_max(2)
    sll_int32   :: n_particles!, degree_smoother
    sll_real64  :: plasma_beta
    character(len=256)   :: sampling_case
    character(len=256)   :: splitting_case
    logical     :: with_control_variate
    character(len=256)   :: initial_distrib
    character(len=256)   :: initial_bfield
    sll_int32   :: spline_degree
    character(len=256) :: file_prefix
    logical :: make_ctest = .false.
    character(len=256) :: ctest_ref_file
    
    sll_int32   :: input_file
    sll_int32   :: ierr, j
    

    namelist /sim_params/         delta_t, n_time_steps, beta, initial_distrib, initial_bfield, plasma_beta
    namelist /output/  file_prefix
    
    namelist /grid_dims/          ng_x, x_min, x_max

    namelist /pic_params/         n_particles, sampling_case, splitting_case, spline_degree, with_control_variate
    namelist /ctest/ make_ctest, ctest_ref_file

    ! Read parameters from file
    open(newunit = input_file, file=trim(filename), IOStat=io_stat)
    if (io_stat /= 0) then
       print*, 'init_pic_2d3v() failed to open file ', filename
       STOP
    end if

    read(input_file, sim_params)

    call sll_s_initial_distribution_new( trim(initial_distrib), [2,3], input_file, sim%init_distrib_params )

    read(input_file, output)
    read(input_file, grid_dims)
    read(input_file, pic_params)
    read(input_file, ctest)
    close (input_file)

    ! Output and ctest
    sim%file_prefix = file_prefix
    sim%make_ctest = make_ctest
    sim%ctest_ref_file = ctest_ref_file
   
    ! Set MPI parameters
    sim%world_size = sll_f_get_collective_size(sll_v_world_collective)
    sim%rank = sll_f_get_collective_rank(sll_v_world_collective)

    ! Copy the read parameters into the simulation parameters
    sim%delta_t = delta_t
    sim%n_time_steps = n_time_steps

    sim%beta = beta

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

    sim%plasma_beta = plasma_beta
    
    sim%n_gcells = ng_x
    sim%n_totaldofs = product(ng_x)
    sim%domain(:,1) = x_min
    sim%domain(:,2) = x_max
    sim%domain(:,3) = x_max-x_min
    sim%delta_x = sim%domain(:,3)/sim%n_gcells
    sim%volume = product(sim%delta_x)
    
    sim%n_particles = n_particles/sim%world_size
    sim%degree_smoother = spline_degree

    call sim%sampler%init( trim(sampling_case), [2,3], sim%n_particles, sim%rank)
    sim%n_total_particles = sim%n_particles * sim%world_size
 
    ! Control variate
    if (with_control_variate .eqv. .true.) then
       allocate(sim%control_variate)
       allocate(sim%control_variate%cv(1))
       call sim%control_variate%cv(1)%init( control_variate_equi, &
            distribution_params=sim%init_distrib_params )
       sim%no_weights = 3
    else
       sim%no_weights = 1
    end if
    

    select case(splitting_case)
    case("splitting_hs")
       sim%splitting_case = sll_p_splitting_hs
    case("splitting_boris")
       sim%splitting_case = sll_p_splitting_boris
    case("splitting_disgradE")
       sim%splitting_case = sll_p_splitting_disgradE
    case("splitting_nldisgradE")
       sim%splitting_case = sll_p_splitting_nldisgradE
    case("splitting_cef")
       sim%splitting_case = sll_p_splitting_cef
    case("splitting_disgradE_trunc")
       sim%splitting_case = sll_p_splitting_disgradE_trunc
    case("splitting_Hpsplit")
       sim%splitting_case = sll_p_splitting_Hpsplit
    case default
       print*, '#splitting case ', splitting_case, ' not implemented.'
    end select

    ! Initialize the particles    (mass set to 1.0 and charge set to -1.0)
    allocate( sim%particle_group )
    sim%particle_group%n_species = 1
    allocate( sll_t_particle_group_2d3v :: sim%particle_group%group(sim%particle_group%n_species) )
    select type ( qp => sim%particle_group%group(1) )
    type is (  sll_t_particle_group_2d3v )
       call qp%init(sim%n_particles, &
            sim%n_total_particles, -1.0_f64, 1.0_f64, sim%no_weights )
    end select

    ! Initialize the field solver
    allocate( sim%maxwell_solver )
    call sim%maxwell_solver%init( sim%domain(:,1:2), sim%n_gcells, sim%degree_smoother )


     ! Initialize kernel smoother
    allocate( sim%particle_mesh_coupling )
    call sim%particle_mesh_coupling%init( sim%n_gcells, sim%domain(:,1:2), sim%n_particles, &
          sim%degree_smoother)

    ! Initialize the arrays for the spline coefficients of the fields
    SLL_ALLOCATE(sim%efield_dofs(sim%n_totaldofs*3), ierr)
    SLL_ALLOCATE(sim%bfield_dofs(sim%n_totaldofs*3), ierr)

    ! Initialize the time-splitting propagator
    if (sim%splitting_case == sll_p_splitting_hs) then
       if (sim%no_weights == 1) then
          call sll_s_new_time_propagator_pic_vm_2d3v_hs(&
               sim%propagator, sim%maxwell_solver, &
               sim%particle_mesh_coupling, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(:,1), sim%domain(:,3), &
               betar=1.0_f64/sim%plasma_beta)
       else
          call sll_s_new_time_propagator_pic_vm_2d3v_hs(&
               sim%propagator, sim%maxwell_solver, &
               sim%particle_mesh_coupling, sim%particle_group, &
               sim%efield_dofs, sim%bfield_dofs, &
               sim%domain(:,1), sim%domain(:,3), &
               control_variate=sim%control_variate, &
               betar=1.0_f64/sim%plasma_beta)
       end if
    else
       print*, 'Propagator not implemented.'
   
    end if
    

  end subroutine init_pic_vm_2d3v

!------------------------------------------------------------------------------!

  subroutine run_pic_vm_2d3v (sim)
    class(sll_t_sim_pic_vm_2d3v_cart), intent(inout) :: sim

    ! Local variables
    sll_int32 :: j, ierr, i_part
    sll_real64, allocatable :: rho(:), rho_local(:), efield_poisson(:)
    sll_int32 :: th_diag_id, dfield_id, efield_id, bfield_id
    sll_real64 :: rhof0
    character(len=256) :: diag_file_name 
    sll_real64 :: wi(1)
    sll_real64 :: xi(3)

    type(sll_t_time_mark) :: start_loop, end_loop

    
    ! Initialize file for diagnostics
    if (sim%rank == 0) then
       diag_file_name = trim(sim%file_prefix)//"_diag.dat"
       
       call sll_s_ascii_file_create(trim(diag_file_name), th_diag_id, ierr)
       call sll_s_ascii_file_create('dfield.dat', dfield_id, ierr)
       call sll_s_ascii_file_create('efield.dat', efield_id, ierr)
       call sll_s_ascii_file_create('bfield.dat', bfield_id, ierr)
    end if

    if (sim%no_weights == 1) then
       call sim%sampler%sample( sim%particle_group%group(1), sim%init_distrib_params, &
            sim%domain(:,1), sim%domain(:,3) )
       rhof0 = 0.0_f64
    else
       call sim%sampler%sample_cv( sim%particle_group%group(1), sim%init_distrib_params, &
            sim%domain(:,1), sim%domain(:,3), sim%control_variate%cv(1) )
       rhof0 = sim%volume  !!! TODO: Implement the general case (also in hamiltonian splitting for j, that way it is hard coded for the initial value)
    end if
       

    ! Set the initial fields
    SLL_ALLOCATE(rho_local(sim%particle_mesh_coupling%n_dofs), ierr)
    SLL_ALLOCATE(rho(sim%particle_mesh_coupling%n_dofs), ierr)
    SLL_ALLOCATE(efield_poisson(sim%particle_mesh_coupling%n_dofs*3), ierr)

    ! Efield  by Poisson
    call solve_poisson( sim%particle_group%group(1), sim%particle_mesh_coupling, sim%maxwell_solver, &
         sim%volume-rhof0, sim%no_weights, rho_local, rho, sim%efield_dofs )

    
    !call sim%particle_group%print('particle.dat')

!!$    ! Bfield = beta*cos(kx): Use b = M{-1}(N_i,beta*cos(kx))
!!$    if ( sim%initial_bfield == sll_p_bfield_cos ) then
!!$       call sim%maxwell_solver%L2projection( beta_cos_k, sim%degree_smoother-1, &
!!$            sim%bfield_dofs)
!!$    else
!!$       call sim%maxwell_solver%L2projection( beta_sin_k, sim%degree_smoother-1, &
!!$            sim%bfield_dofs)
!!$    end if
!!$    sim%bfield_dofs = - sim%bfield_dofs
    sim%bfield_dofs = 0.0_f64

    select case( sim%initial_bfield )
    case (  sll_p_bfield_cos )
       call sim%maxwell_solver%L2projection( beta_cos_k, 3, 2, sim%bfield_dofs(sim%particle_mesh_coupling%n_dofs*2+1:3*sim%particle_mesh_coupling%n_dofs) )
    case ( sll_p_bfield_sin )
       call sim%maxwell_solver%L2projection( beta_sin_k, 3, 2, sim%bfield_dofs(sim%particle_mesh_coupling%n_dofs*2+1:3*sim%particle_mesh_coupling%n_dofs) )
    case ( sll_p_bfield_constant )
       call sim%maxwell_solver%L2projection( beta_constant, 3, 2, sim%bfield_dofs(sim%particle_mesh_coupling%n_dofs*2+1:3*sim%particle_mesh_coupling%n_dofs) )
    end select
       
       
    ! Diagnostics
    call sll_s_time_history_diagnostics_pic_vm_2d3v( &
         sim%particle_group%group(1), sim%maxwell_solver, &
         sim%particle_mesh_coupling,  0.0_f64, &
         sim%degree_smoother, sim%volume, sim%efield_dofs, sim%bfield_dofs, &
         sim%rank, th_diag_id, rho_local, rho )
    
    if (sim%rank == 0 ) then
       call sll_s_set_time_mark( start_loop )
    end if

    
    ! Time loop
    do j=1, sim%n_time_steps
       !print*, 'TIME STEP', j
       ! Strang splitting
       call sim%propagator%strang_splitting( sim%delta_t, 1 )
       
       ! Diagnostics
       !call solve_poisson( sim%particle_group, sim%particle_mesh_coupling, sim%maxwell_solver, &
       !     sim%volume-rhof0, sim%no_weights, rho_local, rho, efield_poisson )
       call sll_s_time_history_diagnostics_pic_vm_2d3v( &
         sim%particle_group%group(1), sim%maxwell_solver, &
         sim%particle_mesh_coupling,  sim%delta_t*real(j,f64), &
         sim%degree_smoother, sim%volume, sim%efield_dofs, sim%bfield_dofs, &
         sim%rank, th_diag_id, rho_local, rho)

    end do

    if (sim%rank == 0 ) then
       call sll_s_set_time_mark( end_loop )
       write(*, "(A, F10.3)") "Main loop run time [s] = ", sll_f_time_elapsed_between( start_loop, end_loop)
       close(th_diag_id)

       if ( sim%make_ctest .eqv. .true. ) then
       ! Check for ctest
          call sll_s_check_diagnostics(trim(sim%ctest_ref_file),trim(diag_file_name) , 1E-13_f64, sim%ctest_passed)
       end if
    end if
    


  contains
    
    function beta_cos_k(x)
      sll_real64             :: beta_cos_k
      sll_real64, intent(in) :: x(2)
      sll_int32 :: r

      r=1

      beta_cos_k = sim%beta * cos(2*sll_p_pi*x(r)/sim%domain(r,3)) 
    end function beta_cos_k

    function beta_sin_k(x)
      sll_real64             :: beta_sin_k
      sll_real64, intent(in) :: x(2)
      sll_int32 :: r

      r=1

      beta_sin_k = sim%beta * sin(2*sll_p_pi*x(r)/sim%domain(r,3)) 
    end function beta_sin_k

    function beta_constant(x)
      sll_real64             :: beta_constant
      sll_real64, intent(in) :: x(2)

      beta_constant = sim%beta
      
    end function beta_constant
    
  end subroutine run_pic_vm_2d3v



!------------------------------------------------------------------------------!

  subroutine delete_pic_vm_2d3v (sim)
    class(sll_t_sim_pic_vm_2d3v_cart), intent(inout) :: sim
    SLL_ASSERT(storage_size(sim)>0)

    call sim%propagator%free()
    deallocate(sim%propagator)
    call sim%particle_group%group(1)%free()
    deallocate (sim%particle_group)
    call sim%maxwell_solver%free()
    deallocate(sim%maxwell_solver)
    call sim%particle_mesh_coupling%free()
    deallocate(sim%particle_mesh_coupling)

    call sim%init_distrib_params%free()
    deallocate(sim%init_distrib_params)
    call sim%sampler%free()

  end subroutine delete_pic_vm_2d3v


!------------------------------------------------------------------------------!
!Diagnostic functions and other helper functions
!> Diagnostics for PIC Vlasov-Maxwell 1d2v 
!> @todo (should be part of the library)
  subroutine sll_s_time_history_diagnostics_pic_vm_2d3v(&
       particle_group, &
       maxwell_solver, &
       particle_mesh_coupling, &
       time, &
       degree, &
       volume, &
       efield_dofs, &
       bfield_dofs, &
       mpi_rank, &
       file_id, &
       scratch1, scratch2)
    class(sll_c_particle_group_base), intent(in) :: particle_group
    type(sll_t_maxwell_2d_fem_fft),     intent(inout) :: maxwell_solver
    type(sll_t_particle_mesh_coupling_spline_2d_feec),     intent(inout) :: particle_mesh_coupling
    sll_real64,                       intent(in) :: time
    sll_real64,                       intent(in) :: efield_dofs(:)
    sll_real64,                       intent(out) :: scratch1(:)
    sll_real64,                       intent(out) :: scratch2(:)
    sll_real64,                       intent(in) :: bfield_dofs(:)
    sll_int32,                        intent(in) :: degree
    sll_int32,                        intent(in) :: mpi_rank
    sll_int32,                        intent(in) :: file_id
    sll_real64,                       intent(in) :: volume  ! product(deltax)

    ! local variables
    sll_real64 :: diagnostics_local(6)
    sll_real64 :: diagnostics(6)
    sll_real64 :: potential_energy(6)
    sll_int32  :: i_part
    sll_real64 :: vi(3),  xi(3)
    sll_real64 :: wi(1)
    sll_real64 :: transfer(1), vvb(1), poynting
    sll_real64 :: efield(2), bfield
    sll_int32  :: n_dofs
    sll_real64 :: error_gauss

    n_dofs = particle_mesh_coupling%n_dofs

    diagnostics_local = 0.0_f64
    do i_part=1,particle_group%n_particles
       vi = particle_group%get_v(i_part)
       xi = particle_group%get_x(i_part)
       wi = particle_group%get_mass(i_part)
       
       ! Kinetic energy
       diagnostics_local(1) = diagnostics_local(1) + &
            (vi(1)**2)*wi(1)
       diagnostics_local(2) = diagnostics_local(2) + &
            (vi(2)**2)*wi(1)
       diagnostics_local(3) = diagnostics_local(3) + &
            (vi(3)**2)*wi(1)
       ! Momentum 1
       diagnostics_local(4) = diagnostics_local(4) + &
            vi(1)*wi(1)
       ! Momentum 2
       diagnostics_local(5) = diagnostics_local(5) + &
            vi(2)*wi(1)
       ! Momentum 3
       diagnostics_local(6) = diagnostics_local(6) + &
            vi(3)*wi(1)
    end do
    diagnostics = 0.0_f64
    call sll_s_collective_reduce_real64(sll_v_world_collective, diagnostics_local, 6,&
         MPI_SUM, 0, diagnostics)
    ! Add ExB part ! TODO
    call compute_e_cross_b ( maxwell_solver, scratch1, efield_dofs, bfield_dofs, diagnostics_local(4:6) )
    !diagnostics(4:6) = diagnostics(4:6) + diagnostics_local(4:6)
    
    !diagnostics(2) = diagnostics(2) + maxwell_solver%inner_product( efield_dofs(:,2), bfield_dofs, degree, degree-1 )
    !diagnostics(3) = diagnostics(3) - maxwell_solver%inner_product( efield_dofs(:,1), bfield_dofs, degree-1 )
    
    ! Check error in Gauss law
    call check_gauss_law ( particle_group, particle_mesh_coupling, maxwell_solver, volume, efield_dofs, scratch1, scratch2, error_gauss )
    
    
    if (mpi_rank == 0) then
       potential_energy(1) = maxwell_solver%inner_product &
            ( efield_dofs(1:n_dofs), efield_dofs(1:n_dofs), 1, 1 )
       potential_energy(2) = maxwell_solver%inner_product &
            ( efield_dofs(n_dofs+1:2*n_dofs), efield_dofs(n_dofs+1:2*n_dofs), 2, 1 )
       potential_energy(3) = maxwell_solver%inner_product &
            ( efield_dofs(n_dofs*2+1:n_dofs*3), efield_dofs(n_dofs*2+1:n_dofs*3), 3, 1 )
       potential_energy(4) = maxwell_solver%inner_product &
            ( bfield_dofs(1:n_dofs), bfield_dofs(1:n_dofs), 1, 2 )
       potential_energy(5) = maxwell_solver%inner_product &
            ( bfield_dofs(n_dofs+1:2*n_dofs), bfield_dofs(n_dofs+1:2*n_dofs), 2, 2 )
       potential_energy(6) = maxwell_solver%inner_product &
            ( bfield_dofs(n_dofs*2+1:n_dofs*3), bfield_dofs(n_dofs*2+1:n_dofs*3), 3, 2 )
       write(file_id,'(f12.5,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16,2g24.16)' ) &
            time,  potential_energy, diagnostics(1:3), &
            sum(diagnostics(1:3)) + sum(potential_energy), diagnostics(4:6), &
            diagnostics_local(4:6), error_gauss
    end if

  end subroutine sll_s_time_history_diagnostics_pic_vm_2d3v

  

  
  !> Accumulate rho and solve Poisson
  subroutine solve_poisson( particle_group, particle_mesh_coupling, maxwell_solver, volume, i_weight, rho_local, rho, efield_dofs )
    class(sll_c_particle_group_base), intent(in) :: particle_group
    type(sll_t_maxwell_2d_fem_fft),     intent(in) :: maxwell_solver
    class(sll_t_particle_mesh_coupling_spline_2d_feec),     intent(inout) :: particle_mesh_coupling
    sll_real64,                       intent(in) :: volume  ! product(deltax)
    sll_int32,                        intent(in) :: i_weight
    sll_real64,                       intent(inout) :: rho_local(:)
    sll_real64,                       intent(inout) :: rho(:)
    sll_real64,                       intent(inout) :: efield_dofs(:)
    
    
    sll_int32 :: i_part
    sll_real64 :: xi(3), wi(1)

    
    rho_local = 0.0_f64
    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x(i_part)
       wi(1) = particle_group%get_charge( i_part , i_weight)
       call particle_mesh_coupling%add_charge(xi, wi(1), &
            [particle_mesh_coupling%spline_degree, &
            particle_mesh_coupling%spline_degree, &
            particle_mesh_coupling%spline_degree],rho_local)
    end do
    
    ! MPI to sum up contributions from each processor
    rho = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, &
         rho_local, &
         particle_mesh_coupling%n_dofs, MPI_SUM, rho)

    ! Add background distribution of neutralizing ions
    rho = rho+volume
    ! Solve Poisson problem
    efield_dofs = 0.0_f64
    call maxwell_solver%compute_E_from_rho( rho, efield_dofs)
    
  end subroutine solve_poisson

  subroutine check_gauss_law ( particle_group, particle_mesh_coupling, maxwell_solver, volume, efield_dofs, rho_gauss, rho, error )
    class(sll_c_particle_group_base), intent(in) :: particle_group
    type(sll_t_maxwell_2d_fem_fft),     intent(in) :: maxwell_solver
    class(sll_t_particle_mesh_coupling_spline_2d_feec),     intent(inout) :: particle_mesh_coupling
    sll_real64,                       intent(in) :: volume  ! product(deltax)
    sll_real64,                       intent(in) :: efield_dofs(:)
    sll_real64,                       intent(inout) :: rho_gauss(:)
    sll_real64,                       intent(inout) :: rho(:)
    sll_real64,                       intent(out) :: error
    
    
    sll_int32 :: i_part
    sll_real64 :: xi(3), wi(1)

    
    rho_gauss = 0.0_f64
    do i_part = 1, particle_group%n_particles
       xi = particle_group%get_x(i_part)
       wi(1) = particle_group%get_charge( i_part)
       call particle_mesh_coupling%add_charge(xi, wi(1), &
            [particle_mesh_coupling%spline_degree, &
            particle_mesh_coupling%spline_degree, &
            particle_mesh_coupling%spline_degree], rho_gauss)
    end do
    
    ! MPI to sum up contributions from each processor
    rho = 0.0_f64
    call sll_o_collective_allreduce( sll_v_world_collective, &
         rho_gauss, &
         particle_mesh_coupling%n_dofs, MPI_SUM, rho)

    ! Add background distribution of neutralizing ions
    rho = rho+volume

    ! Solve Gauss law
    call maxwell_solver%compute_rho_from_e( efield_dofs, rho_gauss )

    error = maxval(abs(rho - rho_gauss ))

  end subroutine check_gauss_law


  !------------------------------------------------------------------------------!

  !> As a control variate, we use the equilibrium (v part of the initial distribution)
  function control_variate_equi( this, xi, vi, time) result(sll_f_control_variate)
    class(sll_t_control_variate) :: this
    sll_real64, optional,  intent( in ) :: xi(:) !< particle position
    sll_real64, optional,  intent( in ) :: vi(:) !< particle velocity
    sll_real64, optional,  intent( in ) :: time  !< current time
    sll_real64               :: sll_f_control_variate


    sll_f_control_variate = &
         this%control_variate_distribution_params%eval_v_density( vi )


  end function control_variate_equi


  subroutine compute_e_cross_b ( maxwell, scratch, efield, bfield, e_cross_b )
    type(sll_t_maxwell_2d_fem_fft), intent(inout) :: maxwell
    sll_real64, intent( out ) :: scratch(:)
    sll_real64, intent( in ) :: efield(:)
    sll_real64, intent( in ) :: bfield(:)
    sll_real64, intent( out ) :: e_cross_b(3)


    
    call  maxwell%multiply_mass( [3, 2, 1], bfield(maxwell%n_total*2+1:maxwell%n_total*3), scratch )
    e_cross_b(1) = sum(efield(maxwell%n_total+1:maxwell%n_total*2)*scratch)
    call  maxwell%multiply_mass([3, 1, 2], bfield(maxwell%n_total+1:maxwell%n_total*2), scratch )
    e_cross_b(1) = e_cross_b(1) - sum(efield(maxwell%n_total*2+1:maxwell%n_total*3)*scratch)
    
    call maxwell%multiply_mass( [1, 3, 2], bfield(1:maxwell%n_total), scratch )
    e_cross_b(2) = sum(efield(maxwell%n_total*2+1:maxwell%n_total*3)*scratch)
    call maxwell%multiply_mass( [2, 3, 1], bfield(maxwell%n_total*2+1:maxwell%n_total*3), scratch )
    e_cross_b(2) = e_cross_b(2) - sum(efield(1:maxwell%n_total)*scratch)

    
    call  maxwell%multiply_mass( [2, 1, 3], bfield(maxwell%n_total+1:maxwell%n_total*2), scratch )
    e_cross_b(3) = sum(efield(1:maxwell%n_total)*scratch)
    call  maxwell%multiply_mass( [1, 2, 3], bfield(1:maxwell%n_total), scratch )
    e_cross_b(3) = e_cross_b(3) - sum(efield(maxwell%n_total+1:maxwell%n_total*2)*scratch)
    
    
  end subroutine compute_e_cross_b


  
  subroutine sll_s_check_diagnostics(reffile, simfile, tol_error, passed)
    character(*), intent(in) :: reffile !< Name of reference file (stored in same folder as source file)
    character(*), intent(in) :: simfile !< Name of file with simulation result
    sll_real64, intent(in)   :: tol_error
    logical, intent(out)     :: passed

    sll_real64 :: error
    sll_real64 :: data_sim(4,18)
    sll_real64 :: data_ref(4,18)
    sll_int32  :: file_id


    ! Read simulation result
    open(newunit=file_id, file=simfile, status='old', action='read')
    read(unit=file_id,fmt=*) data_sim
    close(file_id)

    ! Read reference
    open(newunit=file_id, file=reffile, status='old', action='read')
    read(unit=file_id,fmt=*) data_ref
    close(file_id)

    ! Compare
    data_sim = data_sim - data_ref
    error = maxval(abs(data_sim))

    print*, 'Max error in time history diagnostics: ', error
    if (error < tol_error) then
       passed = .true.
    else
       passed = .false.
    end if
    
  end subroutine sll_s_check_diagnostics
  

end module sll_m_sim_pic_vm_2d3v_cart
