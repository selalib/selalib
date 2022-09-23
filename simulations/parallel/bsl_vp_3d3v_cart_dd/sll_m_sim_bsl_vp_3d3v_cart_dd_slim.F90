module sll_m_sim_bsl_vp_3d3v_cart_dd_slim
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_utilities.h"

  use sll_m_advection_6d_lagrange_dd_slim , only : &
       sll_t_advection_6d_lagrange_dd_slim, &
       sll_s_advection_6d_lagrange_dd_slim_init, &
       sll_s_advection_6d_lagrange_dd_slim_free, &
       sll_s_advection_6d_lagrange_dd_slim_fadvect_eta1, &
       sll_s_advection_6d_lagrange_dd_slim_fadvect_eta2, &
       sll_s_advection_6d_lagrange_dd_slim_fadvect_eta3, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta1, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta2, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta3, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta4, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta5, &
       sll_s_advection_6d_lagrange_dd_slim_advect_eta6, &
       sll_s_advection_6d_lagrange_dd_slim_set_eta123

  use sll_m_advection_6d_spline_dd_slim, only : &
       sll_t_advection_6d_spline_dd_slim, &
       sll_s_advection_6d_spline_dd_slim_init, &
       sll_s_advection_6d_spline_dd_slim_free, &
       sll_s_advection_6d_spline_dd_slim_fadvect_eta1, &
       sll_s_advection_6d_spline_dd_slim_fadvect_eta2, &
       sll_s_advection_6d_spline_dd_slim_fadvect_eta3, &
       sll_s_advection_6d_spline_dd_slim_advect_eta4, &
       sll_s_advection_6d_spline_dd_slim_advect_eta5, &
       sll_s_advection_6d_spline_dd_slim_advect_eta6

  use sll_m_ascii_io, only: &
    sll_s_ascii_file_create

  use sll_m_collective, only : &
       sll_f_create_collective, &
       sll_f_get_collective_size, &
       sll_f_get_collective_rank, &
       sll_s_collective_barrier, &
       sll_s_collective_bcast_3d_real64, &
       sll_t_collective_t, &
       sll_v_world_collective

  use sll_m_sim_base, only : &
       sll_c_simulation_base_class
       
!  use sll_m_constants
  use sll_m_remapper, only : &
       sll_f_new_layout_3d, &
       sll_o_compute_local_sizes, &
       sll_o_initialize_layout_with_distributed_array, &
       sll_t_layout_3d

  use sll_m_poisson_3d_periodic_par, only : &
       sll_s_poisson_3d_periodic_par_compute_e_from_phi, &
       sll_s_poisson_3d_periodic_par_init, &
       sll_s_poisson_3d_periodic_par_solve, &
       sll_t_poisson_3d_periodic_par

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_6d

  use sll_m_decomposition, only: &
       sll_f_apply_halo_exchange, &
       sll_t_cartesian_topology_6d, &
       sll_f_new_cartesian_topology_6d, &
       sll_f_new_cartesian_topology_3d_from_6d, &
       sll_t_cartesian_topology_3d, &
       sll_t_decomposition_slim_6d, &
       sll_f_new_cartesian_domain_decomposition_slim_6d, &
       sll_f_select_dim, &
       sll_f_set_process_grid

  use sll_m_lagrange_interpolator_1d, only: &
       sll_p_lagrange_fixed, &
       sll_p_lagrange_centered

  use sll_m_sim_6d_utilities, only : &
       sll_s_compute_charge_density_6d_dd, &
       sll_s_compute_charge_density_6d_dd_slim, &
       sll_s_time_history_diagnostics, &
       sll_s_check_diagnostics, &
       sll_s_write_simulation_info, &
       sll_f_check_triggered_shutdown

  use sll_m_sim_6d_utilities, only : &
       sll_t_clocks, sll_t_stopwatch, &
       sll_s_init_clocks, sll_s_finalize_clocks, &
       sll_s_start_clock, sll_s_stop_clock

  use sll_m_distribution_function_initializer_6d, only: &
       sll_s_set_local_grid,  &
       sll_t_array, &
       sll_c_distribution_params_6d, &
       sll_s_distribution_params_6d_new, &
       sll_s_distribution_initializer_6d

  use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_f_time_elapsed_since, &
       sll_t_time_mark

  use sll_m_utilities, only: sll_f_query_environment

#ifdef _OPENMP
  use omp_lib
#endif

#ifdef USE_FMEMPOOL
  use fmempool
#endif

  implicit none

  public :: sll_t_sim_bsl_vp_3d3v_cart_dd_slim

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  sll_int32, parameter :: nd = 6
  sll_int32, parameter :: sll_p_splines = 1020

  type, extends(sll_c_simulation_base_class) :: sll_t_sim_bsl_vp_3d3v_cart_dd_slim
     ! Parallel environment parameters
     ! MPI
     sll_int32 :: mpi_world_size
     sll_int32 :: mpi_rank

     ! Physical parameters
     sll_real64 :: domain(6,2) ! Domain size: First column lower bound, second column upper bound.

     ! Simulation parameters
     sll_real64 :: delta_t ! Time step
     sll_real64 :: delta_t_max ! Maximum time step
     sll_real64 :: time_final ! Simulation time
     sll_int32  :: n_iterations !< Number of time steps
     sll_int32  :: n_cells(6) ! Number of cells per dimension
     sll_int32  :: n_pts(6) ! Number of points per dimension
     sll_int32  :: lagrange_width(2) ! Half of stencil width for lagrange interpolation
     sll_int32  :: n_disp(6) ! Maximum allowed number of cells that points are displaced along each dimension.
     logical    :: restart  ! Is it a restart? yes/no
     sll_int32  :: first_time_step !< Number of first time steps (1 usually, >1 for restart)
     character(len=256) :: restart_filename

     ! Mesh parameters
     type(sll_t_cartesian_mesh_6d) :: mesh6d
     type(sll_t_array ) :: etas(6)
     sll_real64, allocatable :: eta1(:)
     sll_real64, allocatable :: eta2(:)
     sll_real64, allocatable :: eta3(:)
     sll_real64, allocatable :: eta4(:)
     sll_real64, allocatable :: eta5(:)
     sll_real64, allocatable :: eta6(:)

     ! For domain decomposition
     sll_int32 :: procs_per_dimension(6)   ! Distribution of the dimensions
     logical :: periodic(6)    ! Periodic boundary conditions?
     type(sll_t_cartesian_topology_6d), pointer :: topology ! Processor topology

     type(sll_t_decomposition_slim_6d), pointer :: decomposition ! Decomposition object for 6d distribution function
     sll_int32 :: global_grid_points_per_dimension(6)
     sll_int32 :: halo_width_per_dimension(6)

     ! Additional part for 6d <-> 3d reduction and Poisson
     type(sll_t_cartesian_topology_3d), pointer :: topology_3d_velocity
     type(sll_t_cartesian_topology_3d), pointer :: topology_3d_spatial
     type(sll_t_collective_t), pointer :: collective_6d
     type(sll_t_collective_t), pointer :: collective_3d_velocity
     type(sll_t_collective_t), pointer :: collective_3d_spatial
     ! attempt to wire the domain decomposition with the existing poisson solver
     type(sll_t_layout_3d), pointer :: layout_3d
     type(sll_t_poisson_3d_periodic_par) :: poisson


     ! Distribution function
     sll_real64, pointer :: f6d(:,:,:,:,:,:)
     ! Phi, efields
     sll_real64, allocatable :: rho(:,:,:)
     sll_real64, allocatable :: phi(:,:,:)
     sll_real64, allocatable :: ex(:,:,:)
     sll_real64, allocatable :: ey(:,:,:)
     sll_real64, allocatable :: ez(:,:,:)

     type(sll_t_advection_6d_lagrange_dd_slim) :: ladvector
     type(sll_t_advection_6d_spline_dd_slim) :: sadvector

     ! For H_p part
     sll_real64, allocatable :: disp_eta1(:)
     sll_real64, allocatable :: disp_eta2(:)
     sll_real64, allocatable :: disp_eta3(:)

     ! Physical parameters
     class(sll_c_distribution_params_6d ), allocatable :: params

     !> For diagnostics, input/output
     sll_int32 :: thdiag_file_id
     logical :: ctest
     character(len=256) :: ctest_ref_file
     character(len=256) :: out_file_prefix
     character(len=256) :: add_file_prefix
     character(len=256) :: nml_filename

     !> Specify initial value
     sll_int32 :: test_case

     sll_int32 :: advector_type

     type(sll_t_clocks) :: clocks
     sll_int32 :: n_diagnostics

     logical :: time_in_phase = .false.

     contains
       procedure :: run => run_6d_vp_dd
       procedure :: init_from_file => init_6d_vp_dd
       procedure :: delete => delete_6d_vp_dd
       procedure :: get_distribution
       procedure :: set_distribution
       procedure :: get_local_size
       procedure :: advect_x
       procedure :: advect_v

    end type sll_t_sim_bsl_vp_3d3v_cart_dd_slim


  contains

    subroutine get_distribution( sim, distribution )
      class(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), intent(inout) :: sim
      sll_real64, pointer,  intent( out ) :: distribution(:,:,:,:,:,:)

      distribution => sim%f6d

    end subroutine get_distribution


    subroutine set_distribution( sim, distribution )
      class(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), intent(inout) :: sim
      sll_real64, pointer,  intent( in ) :: distribution(:,:,:,:,:,:)

      sim%f6d => distribution

    end subroutine set_distribution

    subroutine get_local_size(sim, local_size)
      class(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), intent(inout) :: sim
      sll_int32,  intent( out ) :: local_size(6)

      sll_int32 :: i

      do i=1,6
         local_size(i) = sim%decomposition%local%nw(i)
      end do

    end subroutine get_local_size


    subroutine init_6d_vp_dd(sim, filename)
      intrinsic :: trim
      class(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), intent(inout) :: sim
      character(len=*), intent(in) :: filename

      ! local variables
      sll_int32 :: input_file
      sll_int32 :: IO_stat
      sll_int32 :: aux_data_idx_mn(6)
      ! For namelists
      sll_real64 :: final_time
      sll_real64 :: delta_t
      sll_int32  :: stencil
      sll_int32  :: stencil_x
      character(len=256)   :: bc_type
      character(len=256)   :: interpolator_type
      character(len=256)   :: test_case
      character(len=256)   :: file_prefix
      character(len=256)   :: ctest_ref_file
      character(len=256)   :: restart_filename
      sll_int32  :: restart_itime = 0
      logical    :: restart = .false.
      sll_int32  :: max_disp
      sll_int32  :: num_cells_x1
      sll_int32  :: num_cells_x2
      sll_int32  :: num_cells_x3
      sll_int32  :: num_cells_x4
      sll_int32  :: num_cells_x5
      sll_int32  :: num_cells_x6
      sll_real64 :: v_max
      sll_real64 :: x1_max
      sll_real64 :: x2_max
      sll_real64 :: x3_max
      sll_int32 :: ierr
      logical :: keep_dim(6)
      sll_int32  :: loc_sizes(3)
      logical    :: ctest
      sll_int32 :: omp_world_size
      sll_int32 :: n_iterations
      sll_int32 :: n_diagnostics = 1
      sll_int32  :: n_blocks(nd), block_dim(nd), process_grid(nd)
      logical   :: time_in_phase = .true.

      ! Define namelist for input files
      namelist /sim_params/ final_time, delta_t, ctest, ctest_ref_file, test_case
      namelist /sim_params/ max_disp
      namelist /sim_params/ n_iterations
      namelist /restart_params/ restart, restart_filename, restart_itime, n_diagnostics
      namelist /grid_dims/ num_cells_x1, num_cells_x2, num_cells_x3
      namelist /grid_dims/ num_cells_x4, num_cells_x5, num_cells_x6
      namelist /domain_dims/ v_max, x1_max, x2_max, x3_max
      namelist /advect_params/ bc_type, stencil, interpolator_type, stencil_x
      namelist /output/ file_prefix, time_in_phase
      namelist /parallel_params/ n_blocks, block_dim, process_grid

      aux_data_idx_mn(:) = 1

      sim%nml_filename = filename
      n_iterations = -1

      ! Read the parameters from file
      open(newunit = input_file, file=trim(filename), iostat=IO_stat, status='old', action='read')
      if( IO_stat /= 0 ) then
         print *, 'init_6d_vpB_dd_slim() failed to open file ', filename
         STOP
      end if
      read(input_file, sim_params)
      read(input_file, grid_dims)
      read(input_file, domain_dims)
      read(input_file, advect_params)
      read(input_file, output)

      sim%restart = restart
      sim%restart_filename = restart_filename
      sim%first_time_step = restart_itime
      sim%n_diagnostics = n_diagnostics
      sim%time_in_phase = time_in_phase

      if ( sim%restart .EQV. .false.) then
            sim%first_time_step=1
      endif

      sim%out_file_prefix = file_prefix
      sim%add_file_prefix = trim(file_prefix) // 'add'

      select case(interpolator_type)
      case("fixed")
         sim%advector_type =  sll_p_lagrange_fixed
      case("centered")
         sim%advector_type = sll_p_lagrange_centered
      case("spline", "splines")
         sim%advector_type = sll_p_splines
      case default
         SLL_ERROR('init_bsl_vp_3d3v_dd_slim', 'Interpolator type not implemented.')
      end select

      !> parallelization/overlap parameters
      process_grid = [0,0,0,0,0,0]  ! "empty" default, use internal heuristics
      read(input_file, parallel_params)

      call sll_s_distribution_params_6d_new( sim%params, trim(test_case), input_file )
      
      close(input_file)

      ! set up parallel environment
      sim%mpi_world_size = sll_f_get_collective_size(sll_v_world_collective)
      sim%mpi_rank = sll_f_get_collective_rank(sll_v_world_collective)

      ! iff n_iterations is set it shall take preference over final_time
      if (n_iterations < 0) then
        sim%n_iterations = nint(final_time/delta_t, i32)
      else
        sim%n_iterations = n_iterations
        final_time = real(n_iterations,f64) * delta_t
      endif
      sim%time_final = final_time
      sim%delta_t_max = delta_t
      sim%delta_t = sim%delta_t_max
      sim%ctest = ctest
      sim%ctest_ref_file = ctest_ref_file

      sim%n_cells = [num_cells_x1, num_cells_x2, num_cells_x3, &
           num_cells_x4, num_cells_x5, num_cells_x6]
      sim%n_pts = sim%n_cells
      sim%domain(:,1) = [0.0_f64, 0.0_f64, 0.0_f64, -v_max, -v_max, -v_max]
      sim%domain(:,2) = [x1_max, x2_max, x3_max, v_max, v_max, v_max]

      call sim%mesh6d%init([num_cells_x1, num_cells_x2, num_cells_x3, &
           num_cells_x4, num_cells_x5, num_cells_x6], &
           sim%domain(:,1), sim%domain(:,2))

      sim%lagrange_width = [stencil_x, stencil]
      sim%n_disp = max_disp

      sim%procs_per_dimension = sll_f_set_process_grid(sim%mpi_world_size, process_grid)

      if (sim%mpi_rank == 0) then
        write(*,*) "Running 6D Vlasov simulation with MPI domain decomposition (slim) ..."
#ifdef SLL_GIT_VERSION
        write(*,'(A)') " git: " // SLL_GIT_VERSION
#endif
#ifdef SLL_COMP
        write(*,'(A)') " compiler: " // &
          SLL_COMP
#endif
#ifdef _OPENMP
        omp_world_size = omp_get_max_threads()
#else
        omp_world_size = 1
#endif
        write(*,"(A,I5,A,I3)") "parallelism :: MPI procs: ", sim%mpi_world_size, ", OMP threads: ", omp_world_size
        write(*,"(A,6I4)")     "               process grid: ", sim%procs_per_dimension(1:6)
      end if
      call sll_s_init_clocks(sim%clocks)

#ifdef USE_FMEMPOOL
      call mp_init()
#endif

      !> (1) Create a topology based on the distribution.  In addition, a nd boolean array
      !> is passed to indicate if there is a periodic BC associated to a certain dimension.
      sim%periodic(:) = .true.
      sim%topology => &
           sll_f_new_cartesian_topology_6d(sll_v_world_collective, sim%procs_per_dimension, sim%periodic)
      !> derive a sll_collective from the topologies' MPI communicator
      sim%collective_6d => sll_f_create_collective(sim%topology%comm)


      !> Create 3D sub-topologies from the 6D topology, keeping the velocity space 3D topology.
      keep_dim(1:3) = .false.
      keep_dim(4:6) = .true.
      sim%topology_3d_velocity => &
           sll_f_new_cartesian_topology_3d_from_6d(sim%topology, keep_dim)
      !> derive a sll_collective from the topologies' MPI communicator
      sim%collective_3d_velocity => sll_f_create_collective(sim%topology_3d_velocity%comm)


      !> Create 3D sub-topologies from the 6D topology, keeping the spatial 3D topology.
      keep_dim(1:3) = .true.
      keep_dim(4:6) = .false.
      sim%topology_3d_spatial => &
           sll_f_new_cartesian_topology_3d_from_6d(sim%topology, keep_dim)
      !> derive a sll_collective from the topologies' MPI communicator
      sim%collective_3d_spatial => sll_f_create_collective(sim%topology_3d_spatial%comm)


      !> (2) Create a domain decomposition, create a global array distributed over the
      !> MPI processes.  The size of the local box of the array is passed via an nd array
      !> as well as the widh of the halo to be exchanged between neighboring processes.
      !> The decomposition object contains all the information necessary to allocate and
      !> access the local array (ie it has all the indices), see below.
      sim%global_grid_points_per_dimension(:) = sim%n_cells

      ! Let us allow for a displacement smaller equal to sim%n_disp*dx. For a interpolator with we need lagrange%d points in each direction around the interval. Together this gives a required halo width of lagrange%d+sim%n_disp-1. (lagrange%d = sim%lagrange_width)
      !sim%halo_width_per_dimension(1:3) = sim%n_disp(1:3) - 1 + sim%lagrange_width(1)
      !sim%halo_width_per_dimension(4:6) = sim%n_disp(4:6) - 1 + sim%lagrange_width(2)
      ! We use a fixed stencil width of lagrange_width points around the displaced point, i.e. we need (lagrange_width-1)/2 points to the left and right respectively.
      sim%halo_width_per_dimension(1:3) = (sim%lagrange_width(1)-1)/2
      sim%halo_width_per_dimension(4:6) = (sim%lagrange_width(2)-1)/2
      sim%decomposition => &
           sll_f_new_cartesian_domain_decomposition_slim_6d(sim%topology, &
                    sim%global_grid_points_per_dimension)

      if (sim%mpi_rank == 0) then
        write(*,"(A,6I4)")     "         process-local grid: ", sim%decomposition%local%nw(1:6)
      end if

      !> Allocate the local array.
      allocate( sim%f6d(sim%decomposition%local%mn(1):sim%decomposition%local%mx(1), &
                        sim%decomposition%local%mn(2):sim%decomposition%local%mx(2), &
                        sim%decomposition%local%mn(3):sim%decomposition%local%mx(3), &
                        sim%decomposition%local%mn(4):sim%decomposition%local%mx(4), &
                        sim%decomposition%local%mn(5):sim%decomposition%local%mx(5), &
                        sim%decomposition%local%mn(6):sim%decomposition%local%mx(6)),&
                stat=ierr )
      SLL_ASSERT( ierr == 0 )

      ! Allocate the local array for rho
      allocate( sim%rho(sim%decomposition%local%mn(1):sim%decomposition%local%mx(1), &
                        sim%decomposition%local%mn(2):sim%decomposition%local%mx(2), &
                        sim%decomposition%local%mn(3):sim%decomposition%local%mx(3)),&
                stat=ierr )
      SLL_ASSERT( ierr == 0 )

      ! Allocate the local array for ex
      allocate( sim%ex(sim%decomposition%local%mn(1):sim%decomposition%local%mx(1), &
                       sim%decomposition%local%mn(2):sim%decomposition%local%mx(2), &
                       sim%decomposition%local%mn(3):sim%decomposition%local%mx(3)),&
                stat=ierr )
      SLL_ASSERT( ierr == 0 )

      ! Allocate the local array for ey
      allocate( sim%ey(sim%decomposition%local%mn(1):sim%decomposition%local%mx(1), &
                       sim%decomposition%local%mn(2):sim%decomposition%local%mx(2), &
                       sim%decomposition%local%mn(3):sim%decomposition%local%mx(3)),&
                stat=ierr )
      SLL_ASSERT( ierr == 0 )

      ! Allocate the local array for ez
      allocate( sim%ez(sim%decomposition%local%mn(1):sim%decomposition%local%mx(1), &
                       sim%decomposition%local%mn(2):sim%decomposition%local%mx(2), &
                       sim%decomposition%local%mn(3):sim%decomposition%local%mx(3)),&
                stat=ierr )
      SLL_ASSERT( ierr == 0 )

      ! select a subset of the processors to do the Poisson solve step
    !!  if (all(sim%topology_3d_velocity%coords == 0)) then
    !!     SLL_ASSERT_ALWAYS( sim%collective_3d_velocity%rank == 0 )
         sim%layout_3d => sll_f_new_layout_3d( sim%collective_3d_spatial )
         call sll_o_initialize_layout_with_distributed_array( &
              sim%mesh6d%num_cells(1:3), &
              sim%decomposition%local%mn(1:3), &
              sim%decomposition%local%mx(1:3), sim%layout_3d )
         call sll_o_compute_local_sizes( sim%layout_3d, &
              loc_sizes(1), loc_sizes(2), loc_sizes(3) )

         ! initialize the existing Poisson solver
         call sll_s_poisson_3d_periodic_par_init( &
                      sim%layout_3d, &
                      sim%global_grid_points_per_dimension(1), &
                      sim%global_grid_points_per_dimension(2), &
                      sim%global_grid_points_per_dimension(3), &
                      sim%mesh6d%eta_max(1)- sim%mesh6d%eta_min(1), &
                      sim%mesh6d%eta_max(2)- sim%mesh6d%eta_min(2), &
                      sim%mesh6d%eta_max(3)- sim%mesh6d%eta_min(3), &
                      sim%poisson )
         ! print*, poisson%loc_sizes(1,:)
         allocate( sim%phi(1:sim%poisson%loc_sizes(1,1), &
                           1:sim%poisson%loc_sizes(1,2), &
                           1:sim%poisson%loc_sizes(1,3)),&
                   stat=ierr )
         SLL_ASSERT( ierr == 0 )
         sim%phi = 0.0_f64
   !!   end if

      ! Initialize advector
      if (sim%advector_type ==  sll_p_lagrange_fixed .OR. sim%advector_type ==  sll_p_lagrange_centered) then
         call sll_s_advection_6d_lagrange_dd_slim_init( &
              sim%ladvector, sim%lagrange_width )
      end if

      allocate( sim%disp_eta1(sim%decomposition%local%nw(4)) )
      allocate( sim%disp_eta2(sim%decomposition%local%nw(5)) )
      allocate( sim%disp_eta3(sim%decomposition%local%nw(6)) )


      if(sll_f_get_collective_rank(sll_v_world_collective)==0) then
         call sll_s_ascii_file_create( trim(sim%out_file_prefix)//'.dat', sim%thdiag_file_id, ierr)
      endif



      call sll_s_set_local_grid( sim%decomposition%local%nw, &
           sim%decomposition%local%mn, &
           sim%mesh6d%eta_min, &
           sim%mesh6d%delta_eta, &
           sim%etas )

      call sll_s_distribution_initializer_6d (  &
           sim%decomposition%local%nw, &
           aux_data_idx_mn, &
           sim%params, sim%etas, sim%f6d )


! SLIM : disabled
!      call sll_f_apply_halo_exchange(sim%topology, sim%decomposition, sim%f6d)

      sim%disp_eta1 = -sim%etas(4)%vals*sim%delta_t/sim%mesh6d%delta_eta(1)
      sim%disp_eta2 = -sim%etas(5)%vals*sim%delta_t/sim%mesh6d%delta_eta(2)
      sim%disp_eta3 = -sim%etas(6)%vals*sim%delta_t/sim%mesh6d%delta_eta(3)

      select case ( sim%advector_type )
      case( sll_p_lagrange_centered )
         ! Note: n_disp(1:3) need to be identical
         sim%n_disp(1:3) = sim%lagrange_width(1)/2
         call sll_s_advection_6d_lagrange_dd_slim_set_eta123(sim%ladvector, sim%decomposition, &
              sim%disp_eta1, sim%disp_eta2, sim%disp_eta3 )
      case( sll_p_splines )
         call sll_s_advection_6d_spline_dd_slim_init (sim%sadvector, sim%decomposition, &
              sim%disp_eta1, sim%disp_eta2, sim%disp_eta3 )
      end select

      call sll_s_compute_charge_density_6d_dd_slim(sim%f6d, &
                 sim%decomposition, &
                 sim%rho, &
                 sim%collective_3d_velocity, &
                 sim%mesh6d%volume_eta456)

      ! Compute electric fields
      !if (all(sim%topology_3d_velocity%coords == 0)) then
          ! call the existing Poisson solver
         call sll_s_poisson_3d_periodic_par_solve( sim%poisson, sim%rho, sim%phi )
         call sll_s_poisson_3d_periodic_par_compute_e_from_phi( sim%poisson, sim%phi, &
              sim%ex, sim%ey, sim%ez)
      !end if

      !call sll_s_collective_bcast_3d_real64(sim%collective_3d_velocity, sim%ex, 0)
      !call sll_s_collective_bcast_3d_real64(sim%collective_3d_velocity, sim%ey, 0)
      !call sll_s_collective_bcast_3d_real64(sim%collective_3d_velocity, sim%ez, 0)

      call sll_s_time_history_diagnostics(&
                 sim%decomposition%local%mn, &
                 sim%decomposition%local%mx, &
                 sim%decomposition%local%mn, &  ! previously: `sim%decomposition%local%lo`
                 0.0_f64, &
                 sim%mesh6d%volume, &
                 sim%mesh6d%volume_eta123, &
                 sim%etas,&
                 sim%f6d, &
                 sim%rho, &
                 sim%phi, &
                 sim%ex, &
                 sim%ey, &
                 sim%ez, &
                 sim%thdiag_file_id, &
                 sim%topology_3d_velocity%coords)

    end subroutine init_6d_vp_dd


    subroutine run_6d_vp_dd(sim)
      class(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), intent(inout) :: sim

      sll_int32 :: itime
      type(sll_t_time_mark) :: t0, t1, t0_step, t1_step
      sll_real64 :: wall_time(1)

      ! load initial value


!      call sll_f_apply_halo_exchange(sim%topology, sim%decomposition, sim%f6d, sll_f_select_dim(5))
!      call sll_s_advection_6d_lagrange_dd_slim_advect_eta5(sim%advector, -sim%ey*sim%delta_t*0.5_f64/sim%mesh6d%delta_eta(5), sim%f6d)
!
!      call sll_f_apply_halo_exchange(sim%topology, sim%decomposition, sim%f6d, sll_f_select_dim(6))
      !      call sll_s_advection_6d_lagrange_dd_slim_advect_eta6(sim%advector, -sim%ez*sim%delta_t*0.5_f64/sim%mesh6d%delta_eta(6), sim%f6d)


      ! ------------------------------------------------------------------------
      !
      !                                MAIN LOOP
      !
      ! ------------------------------------------------------------------------


      call sim%advect_v(0.5_f64*sim%delta_t)
      
      if (sim%mpi_rank == 0 ) then
         write(*,*) "Entering main loop ... "
         call sll_s_set_time_mark(t0)
      end if

      do itime = sim%first_time_step, (sim%first_time_step + sim%n_iterations - 1)
        if (sim%mpi_rank == 0 ) &
          call sll_s_set_time_mark(t0_step)

        call sim%advect_x()


        if (sll_f_query_environment("SLL_PRE_POISSON_BARRIER", .false.)) then
          call sll_s_collective_barrier(sll_v_world_collective)
        endif
        call sll_s_start_clock(sim%clocks, 'P')  ! total time for Poisson

        call sll_s_start_clock(sim%clocks, 'PC')  ! Poisson, charge part
        call sll_s_compute_charge_density_6d_dd_slim(sim%f6d, &
                   sim%decomposition, &
                   sim%rho, &
                   sim%collective_3d_velocity, &
                   sim%mesh6d%volume_eta456)
        call sll_s_stop_clock(sim%clocks, 'PC')

        call sll_s_start_clock(sim%clocks, 'PF')  ! Poisson, fields part
        ! Compute electric fields
        !if (all(sim%topology_3d_velocity%coords == 1)) then
           ! call the existing Poisson solver
           call sll_s_poisson_3d_periodic_par_solve( sim%poisson, sim%rho, sim%phi )

           call sll_s_poisson_3d_periodic_par_compute_e_from_phi( sim%poisson, sim%phi, &
                sim%ex, sim%ey, sim%ez)
        !end if
        call sll_s_stop_clock(sim%clocks, 'PF')


        call sll_s_stop_clock(sim%clocks, 'P')


        call sll_s_start_clock(sim%clocks, 'D')

        if ( mod(itime, int( sim%n_diagnostics)) == 0) then

          call sll_s_time_history_diagnostics(&
                    sim%decomposition%local%mn, &
                    sim%decomposition%local%mx, &
                    sim%decomposition%local%mn, &  ! previously: `sim%decomposition%local%lo`
                    real(itime, f64)*sim%delta_t, &
                    sim%mesh6d%volume, &
                    sim%mesh6d%volume_eta123, &
                    sim%etas,&
                    sim%f6d, &
                    sim%rho, &
                    sim%phi, &
                    sim%ex, &
                    sim%ey, &
                    sim%ez, &
                    sim%thdiag_file_id, &
                    sim%topology_3d_velocity%coords)

        endif

        call sll_s_stop_clock(sim%clocks, 'D')

        if ((sim%time_in_phase .eqv. .true.) .and. (itime == (sim%first_time_step + sim%n_iterations - 1)) ) then
           call sim%advect_v(0.5_f64*sim%delta_t)
        else
           call sim%advect_v(sim%delta_t)
        end if
           

        if (sim%mpi_rank == 0) then
          call sll_s_set_time_mark(t1_step)
          write (*, "(A, F7.3, A, F7.3, A, F7.3)") " Time ", real(itime, f64)*sim%delta_t, " of ", sim%time_final, &
                " :: step run time [s] = ", sll_f_time_elapsed_between(t0_step, t1_step)
#ifdef USE_FMEMPOOL
          call mp_statistics()
#endif
        end if

        ! Check if a file named "./stop" is present and perform a safe shutdown.
        if (sll_f_check_triggered_shutdown()) then
          exit
        endif
      end do  ! time loop
      sim%first_time_step = itime

      if (sim%mpi_rank == 0) then
         call sll_s_set_time_mark(t1)

         call sll_s_finalize_clocks(sim%clocks)

         wall_time = sll_f_time_elapsed_between(t0, t1)
         write(*,"(A, F10.3)") " Leaving main loop.  Main loop run time [s] = ", wall_time
         call sll_s_write_simulation_info( sim%out_file_prefix, &
              sim%nml_filename, &
              sim%mesh6d%num_cells, &
              [sim%delta_t], &
              [sim%mpi_world_size], wall_time )
      end if

      if (sim%mpi_rank == 0) then
         close(sim%thdiag_file_id)
      end if

      if (sim%mpi_rank == 0 .and. sim%ctest .eqv. .true.) then
         call sll_s_check_diagnostics(trim(sim%ctest_ref_file), trim(sim%out_file_prefix)//'.dat')
      end if

    end subroutine run_6d_vp_dd


    subroutine delete_6d_vp_dd (sim)
      class(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), intent(inout) :: sim

      sll_int32 :: ierr

#ifdef USE_FMEMPOOL
      call mp_finalize()
#endif
      SLL_DEALLOCATE(sim%topology, ierr)
      SLL_ASSERT(ierr == 0 )
      SLL_DEALLOCATE(sim%decomposition, ierr)
      SLL_ASSERT(ierr == 0 )
    end subroutine delete_6d_vp_dd

    
    subroutine advect_x( sim )
      class(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), intent(inout) :: sim

       select case( sim%advector_type )
        case ( sll_p_lagrange_centered )

           call sll_s_start_clock(sim%clocks, 'X')
           call sll_s_advection_6d_lagrange_dd_slim_fadvect_eta1(sim%ladvector, &
                sim%topology, &
                sim%decomposition, &
                sim%f6d)
           call sll_s_advection_6d_lagrange_dd_slim_fadvect_eta2(sim%ladvector, &
                sim%topology, &
                sim%decomposition, &
                sim%f6d)
           call sll_s_advection_6d_lagrange_dd_slim_fadvect_eta3(sim%ladvector, &
                sim%topology, &
                sim%decomposition, &
                sim%f6d)
           call sll_s_stop_clock(sim%clocks, 'X')

        case ( sll_p_lagrange_fixed )

           call sll_s_start_clock(sim%clocks, 'H1')
           call sll_f_apply_halo_exchange(sim%topology, &
                sim%decomposition, &
                sim%f6d, &
                1, &
                sim%halo_width_per_dimension(1), &
                sim%halo_width_per_dimension(1))
           call sll_s_stop_clock(sim%clocks, 'H1')

           call sll_s_start_clock(sim%clocks, 'X1')
           call sll_s_advection_6d_lagrange_dd_slim_advect_eta1(sim%ladvector, &
                sim%decomposition, &
                sim%disp_eta1, &
                sim%f6d)
           call sll_s_stop_clock(sim%clocks, 'X1')

           call sll_s_start_clock(sim%clocks, 'H2')
           call sll_f_apply_halo_exchange(sim%topology, &
                sim%decomposition, &
                sim%f6d, &
                2, &
                sim%halo_width_per_dimension(2), &
                sim%halo_width_per_dimension(2))
           call sll_s_stop_clock(sim%clocks, 'H2')

           call sll_s_start_clock(sim%clocks, 'X2')
           call sll_s_advection_6d_lagrange_dd_slim_advect_eta2(sim%ladvector, &
                sim%decomposition, &
                sim%disp_eta2, &
                sim%f6d)
           call sll_s_stop_clock(sim%clocks, 'X2')

           call sll_s_start_clock(sim%clocks, 'H3')
           call sll_f_apply_halo_exchange(sim%topology, &
                sim%decomposition, &
                sim%f6d, &
                3, &
                sim%halo_width_per_dimension(3), &
                sim%halo_width_per_dimension(3))
           call sll_s_stop_clock(sim%clocks, 'H3')

           call sll_s_start_clock(sim%clocks, 'X3')
           call sll_s_advection_6d_lagrange_dd_slim_advect_eta3(sim%ladvector, &
                sim%decomposition, &
                sim%disp_eta3, &
                sim%f6d)
           call sll_s_stop_clock(sim%clocks, 'X3')

        case ( sll_p_splines )
           call sll_s_start_clock(sim%clocks, 'X')
           call sll_s_advection_6d_spline_dd_slim_fadvect_eta1(sim%sadvector, &
                sim%topology, &
                sim%decomposition, &
                sim%f6d)
           call sll_s_advection_6d_spline_dd_slim_fadvect_eta2(sim%sadvector, &
                sim%topology, &
                sim%decomposition, &
                sim%f6d)
           call sll_s_advection_6d_spline_dd_slim_fadvect_eta3(sim%sadvector, &
                sim%topology, &
                sim%decomposition, &
                sim%f6d)
           call sll_s_stop_clock(sim%clocks, 'X')
        end select

      

      end subroutine advect_x


      subroutine advect_v( sim, delta_t )
      class(sll_t_sim_bsl_vp_3d3v_cart_dd_slim), intent(inout) :: sim
      sll_real64, intent(in) :: delta_t
      

        select case( sim%advector_type )
        case ( sll_p_splines )
           call sll_s_start_clock(sim%clocks, 'V')
           call sll_s_advection_6d_spline_dd_slim_advect_eta4(sim%sadvector, &
                sim%topology, sim%decomposition, &
                sim%ex*delta_t/sim%mesh6d%delta_eta(4), &
                sim%f6d)
           call sll_s_advection_6d_spline_dd_slim_advect_eta5(sim%sadvector, &
                sim%topology, sim%decomposition, &
                sim%ey*delta_t/sim%mesh6d%delta_eta(5), &
                sim%f6d)
           call sll_s_advection_6d_spline_dd_slim_advect_eta6(sim%sadvector, &
                sim%topology, sim%decomposition, &
                sim%ez*delta_t/sim%mesh6d%delta_eta(6), &
                sim%f6d)
           call sll_s_stop_clock(sim%clocks, 'V')
        case default
           call sll_s_start_clock(sim%clocks, 'H4')
           call sll_f_apply_halo_exchange(sim%topology, &
                sim%decomposition, &
                sim%f6d, &
                4, &
                sim%halo_width_per_dimension(4), &
                sim%halo_width_per_dimension(4))
           call sll_s_stop_clock(sim%clocks, 'H4')

           call sll_s_start_clock(sim%clocks, 'X4')
           call sll_s_advection_6d_lagrange_dd_slim_advect_eta4(sim%ladvector, &
                sim%decomposition, &
                sim%ex*delta_t/sim%mesh6d%delta_eta(4), &
                sim%f6d)
           call sll_s_stop_clock(sim%clocks, 'X4')

           call sll_s_start_clock(sim%clocks, 'H5')
           call sll_f_apply_halo_exchange(sim%topology, &
                sim%decomposition, &
                sim%f6d, &
                5, &
                sim%halo_width_per_dimension(5), &
                sim%halo_width_per_dimension(5))
           call sll_s_stop_clock(sim%clocks, 'H5')

           call sll_s_start_clock(sim%clocks, 'X5')
           call sll_s_advection_6d_lagrange_dd_slim_advect_eta5(sim%ladvector, &
                sim%decomposition, &
                sim%ey*delta_t/sim%mesh6d%delta_eta(5), &
                sim%f6d)
           call sll_s_stop_clock(sim%clocks, 'X5')

           call sll_s_start_clock(sim%clocks, 'H6')
           call sll_f_apply_halo_exchange(sim%topology, &
                sim%decomposition, &
                sim%f6d, &
                6, &
                sim%halo_width_per_dimension(6), &
                sim%halo_width_per_dimension(6))
           call sll_s_stop_clock(sim%clocks, 'H6')

           call sll_s_start_clock(sim%clocks, 'X6')
           call sll_s_advection_6d_lagrange_dd_slim_advect_eta6(sim%ladvector, &
                sim%decomposition, &
                sim%ez*delta_t/sim%mesh6d%delta_eta(6), &
                sim%f6d)
           call sll_s_stop_clock(sim%clocks, 'X6')
        end select

      end subroutine advect_v
      

  end module sll_m_sim_bsl_vp_3d3v_cart_dd_slim
