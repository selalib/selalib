!> @ingroup particle_methods
!> @author MCP ALH
!> @brief Generic simulation algorithms for PIC particles

!> @details Generic simulation algorithm for simple PIC particles of type ::sll_m_simple_pic_4d_group <!--
!> [[selalib:src/particle_methods/pic_remapped/simple_pic/sll_m_simple_pic_4d_group.F90::sll_m_simple_pic_4d_group]]
!> --> and LTPIC particles based on the generic class ::sll_m_remapped_pic_base <!--
!> [[selalib:src/particle_methods/pic_remapped/sll_m_remapped_pic_base.F90::sll_m_remapped_pic_base]] -->.

! (The doxygen page for this simulation is
! [[selalib:doc/build/html/doxygen/html/namespacesll__simulation__4d__vp__generic__pic__cartesian__module.html]]
! produced by [[elisp:(compile "cd ${SELALIB}/build && make doc")]])

module sll_m_sim_pic_vp_2d2v_cart_remapped

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_accumulators.h" 
#include "particle_representation.h"


  use sll_m_constants
  use sll_m_sim_base
  use sll_m_cartesian_meshes
  use sll_m_timer
  use sll_m_poisson_2d_fft
  use sll_m_poisson_2d_base
  use sll_m_gnuplot
  use sll_m_collective
  use sll_m_ascii_io

  use sll_m_remapped_pic_base
  use sll_m_simple_pic_4d_group
  use sll_m_bsl_lt_pic_4d_group
  ! use sll_m_particle_initializers_4d
  ! use sll_m_particle_sort
  use sll_m_charge_to_density
  use sll_m_pic_utilities
  ! use sll_m_representation_conversion

  implicit none

  ! uses [[selalib:src/simulations/simulation_base_class.F90::sll_simulation_base_class]]. We have chosen to store all
  ! simulation parameters (both ltpic and simple) in this sim object.

  ! doxygen page
  ! [[selalib:doc/build/html/doxygen/html/structsll__simulation__4d__vp__generic__pic__cartesian__module_1_1sll__simulation__4d__vp__generic__pic__cartesian.html]]
  ! produced by [[elisp:(compile "cd ${SELALIB}/build && make doc")]])

  !> @ingroup particle_methods

  !> @brief Test simulation applied to PIC particles of type
  !> sll_m_simple_pic_4d_particle::sll_simple_pic_4d_particle

  ! [[selalib:src/particle_methods/pic_remapped/simple_pic/sll_m_simple_pic_4d_particle.F90::sll_simple_pic_4d_particle]]
  
  type, extends(sll_simulation_base_class) :: sll_simulation_4d_vp_generic_pic_cartesian

     !> the abstract particle group
     
     ! <<particle_group>>
     
     class(sll_c_remapped_particle_group),  pointer     :: particle_group

     !> @name Physics/numerical parameters
     !> @{
     
     !> Time step
     sll_real64 :: dt

     !> Total number of iteration
     sll_int32  :: num_iterations

     !> Plot period in terms of number of iterations
     sll_int32  :: plot_period

     !> unused at the moment
     sll_int32  :: remap_period

     !> cf sll_c_remapped_particle_group::set_landau_parameters

     ! [[selalib:src/particle_methods/pic_remapped/bsl_lt_pic/sll_m_bsl_lt_pic_4d_group.F90::bsl_lt_pic_4d_set_landau_parameters]]

     sll_real64 :: thermal_speed_ions
     
     sll_int32  :: number_particles
     sll_int32  :: number_deposition_particles
     logical :: domain_is_x_periodic
     logical :: domain_is_y_periodic
     sll_real64, dimension(1:6) :: elec_params
     !> @}

     !> Underlying 2D cartesian sll_m_cartesian_meshes::sll_cartesian_mesh_2d

     type ( sll_cartesian_mesh_2d ),    pointer :: mesh_2d

     !> called q_accumulator in Sever simulation
     type(sll_charge_accumulator_2d_ptr), dimension(:), pointer     :: q_accumulator_ptr

     !> uses sll_m_accumulators::sll_charge_accumulator_2d
     ! [[selalib:src/particle_methods/pic_2d_standard/pic_accumulators/sll_m_accumulators.F90::sll_charge_accumulator_2d]]
     type(sll_charge_accumulator_2d),     dimension(:), pointer     :: charge_accumulator
     type(electric_field_accumulator),                  pointer     :: E_accumulator
     logical :: use_lt_pic_scheme        ! if false then use pic scheme
     type(sll_charge_accumulator_2d_CS_ptr), dimension(:), pointer  :: q_accumulator_CS
     type(electric_field_accumulator_CS), pointer :: E_accumulator_CS
     sll_real64, dimension(:,:), pointer :: rho
     type(poisson_2d_fft_solver), pointer :: poisson
     sll_real64 :: total_density

     sll_real64, dimension(:,:), pointer :: E1, E2
     sll_int32 :: my_rank
     sll_int32 :: world_size
     sll_int32 :: n_threads
   contains
     procedure, pass(sim) :: init_from_file => init_4d_generic_pic_cartesian
     procedure, pass(sim) :: run => run_4d_generic_pic_cartesian
  end type sll_simulation_4d_vp_generic_pic_cartesian

  !> Standard Selalib deallocation
  interface sll_delete
     module procedure delete_4d_generic_pic_cartesian
  end interface sll_delete
  
contains

  ! <<particles_snapshot>> gnuplot-compatible output of all particles with their position and speed at a given time ALH_BUG Develop the
  ! corresponding gnuplot visualisation script <<ALH>>

  subroutine particles_snapshot(time,sim)

    class(sll_simulation_4d_vp_generic_pic_cartesian), intent(in) :: sim
    
    ! Actual time
    
    sll_real64::time

    ! Local data
    
    sll_int32          :: fileid
    character(len=100) :: filename
    sll_int32          :: error
    sll_int32          :: k
    sll_real64         :: x(3),v(3)
    
    ! creating unique file name [[http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html]]

    write(filename,'(A19,F0.3,A4)') 'particles_snapshot_',time,'.gnu'

    ! Inspired from [[selalib:src/io/file_io/sll_m_gnuplot.F90::subroutine sll_gnuplot_corect_2d]]. Calls
    ! [[selalib:src/io/file_io/sll_m_ascii_io.F90::sll_ascii_file_create]]

    call sll_ascii_file_create(filename,fileid,error)

    ! placing data inside the gnuplot script
    ! [[http://stackoverflow.com/questions/12722048/gnuplot-data-points-in-script]]
    
    write(fileid,*) 'set xrange [',sim%mesh_2d%eta1_min,':',sim%mesh_2d%eta1_max,']'
    write(fileid,*) 'set yrange [',sim%mesh_2d%eta2_min,':',sim%mesh_2d%eta2_max,']'
    write(fileid,*) 'set size ratio -1'
    write(fileid,*) 'plot "-" using 1:2:3:4 with vectors'

    ! Loop over all particles. Do not display all particles for a clearer picture. AAA_ALH_TODO try with less particles.
    
    do k = 1,sim%number_particles,1000

       ! reference to [[particle_group]]
       
       x = sim%particle_group%get_x(k)
       v = sim%particle_group%get_v(k)
       
       write(fileid,*) x(1),x(2),v(1),v(2)
    enddo
    close(fileid)
  end subroutine particles_snapshot

  !> redefines sll_simulation_base_class::init_from_file <!--
  !> [[selalib:src/simulations/simulation_base_class.F90::init_from_file]] -->

  subroutine init_4d_generic_pic_cartesian   ( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_vp_generic_pic_cartesian), intent(inout) :: sim
    character(len=*), intent(in)                          :: filename
    sll_int32   :: IO_stat
    sll_int32   :: ierr
    sll_int32   :: j
    sll_real64  :: dt
    sll_int32   :: number_iterations, plot_period, remap_period
    sll_int32   :: NUM_PARTICLES, GUARD_SIZE, PARTICLE_ARRAY_SIZE
    sll_real64  :: THERM_SPEED
    sll_real64  :: SPECIES_CHARGE, SPECIES_MASS, ALPHA
    logical     :: UseLtPicScheme
    logical     :: DOMAIN_IS_X_PERIODIC, DOMAIN_IS_Y_PERIODIC
    sll_int32   :: NC_X,  NC_Y
    sll_real64  :: XMIN, KX_LANDAU, XMAX, YMIN, YMAX
    sll_real64  :: target_total_charge
    logical     :: enforce_total_charge
    sll_int32   :: particle_group_id
    sll_int32   :: remap_f_type
    sll_int32   :: remap_degree
    sll_real64  :: remapping_grid_vx_min
    sll_real64  :: remapping_grid_vx_max
    sll_real64  :: remapping_grid_vy_min
    sll_real64  :: remapping_grid_vy_max
    sll_int32   :: remapping_cart_grid_number_cells_x
    sll_int32   :: remapping_cart_grid_number_cells_y
    sll_int32   :: remapping_cart_grid_number_cells_vx
    sll_int32   :: remapping_cart_grid_number_cells_vy
    sll_int32   :: remapping_sparse_grid_max_level_x
    sll_int32   :: remapping_sparse_grid_max_level_y
    sll_int32   :: remapping_sparse_grid_max_level_vx
    sll_int32   :: remapping_sparse_grid_max_level_vy
    sll_int32   :: number_deposition_particles
    sll_int32   :: number_markers_x
    sll_int32   :: number_markers_y
    sll_int32   :: number_markers_vx
    sll_int32   :: number_markers_vy
    sll_int32   :: flow_grid_number_cells_x
    sll_int32   :: flow_grid_number_cells_y
    sll_int32   :: flow_grid_number_cells_vx
    sll_int32   :: flow_grid_number_cells_vy

    sll_int32   :: initial_density_identifier

    sll_int32, dimension(4)   :: remapping_sparse_grid_max_levels

    ! <<simple_pic_particle_group>>

    type(sll_simple_pic_4d_group),      pointer     :: simple_pic_particle_group

    ! <<bsl_lt_pic_particle_group>>
    type(sll_bsl_lt_pic_4d_group),      pointer     :: bsl_lt_pic_particle_group

    type(sll_charge_accumulator_2d),    pointer :: charge_accumulator

    sll_real64  :: er, psi, omega_i, omega_r
    sll_int32, parameter  :: input_file = 99
    sll_int32, dimension(:), allocatable  :: rand_seed
    sll_int32   :: rand_seed_size
    sll_int32  :: thread_id

    !> The generic simulation reads all parameters for both simple and ltpic options from the same file
    namelist /sim_params/           THERM_SPEED,            &
                                    dt,                     &
                                    number_iterations,      &
                                    plot_period,            &
                                    SPECIES_CHARGE,         &
                                    SPECIES_MASS,           &
                                    ALPHA,                  &
                                    UseLtPicScheme

    namelist /elec_params/          er, psi, omega_r, omega_i

    ! cartesian grid for the Poisson solver      [previous name was grid_dims]
    namelist /poisson_grid_params/  NC_X, NC_Y, XMIN, KX_LANDAU, YMIN, YMAX, &
                                    DOMAIN_IS_X_PERIODIC, DOMAIN_IS_Y_PERIODIC

    ! MCP put new parameters here on monday nov 9 -- todo: remove this line once the new parameters are used consistently

    ! there are 3 sets of parameters:
    !   - one for the discretization of the remapped f,
    !   - one for the discretization of the deposited f
    !   - one for the discretization of the flow

    ! discretization of the remapped f
    !   -> remapping period, type of interpolation structure and polynomial degree
    !   -> size of the remapping grid / sparse grid
    namelist /lt_pic_remap_params/  remap_period,                       &
                                    remap_f_type,                       &
                                    remap_degree,                       &
                                    remapping_grid_vx_min,              &   ! note: the x-y domain is given by the Poisson grid
                                    remapping_grid_vx_max,              &
                                    remapping_grid_vy_min,              &
                                    remapping_grid_vy_max,              &
                                    remapping_cart_grid_number_cells_x,          &   ! for splines
                                    remapping_cart_grid_number_cells_y,          &   ! for splines
                                    remapping_cart_grid_number_cells_vx,         &   ! for splines
                                    remapping_cart_grid_number_cells_vy,         &   ! for splines
                                    remapping_sparse_grid_max_level_x,           &   ! for the sparse grid
                                    remapping_sparse_grid_max_level_y,           &   ! for the sparse grid
                                    remapping_sparse_grid_max_level_vx,          &   ! for the sparse grid
                                    remapping_sparse_grid_max_level_vy               ! for the sparse grid


    ! discretization of the deposited f
    !   -> number of deposition particles
    namelist /lt_pic_deposition_params/ number_deposition_particles

    ! discretization of the flow:
    !   -> number of markers to be pushed forward
    !   -> size of the cells where the flow is linearized (the flow grid)
    !      note: the bounds of the flow grid are the same as those of the remap grid
    namelist /lt_pic_markers_params/number_markers_x,             &
                                    number_markers_y,             &
                                    number_markers_vx,            &
                                    number_markers_vy,            &
                                    flow_grid_number_cells_x,     &
                                    flow_grid_number_cells_y,     &
                                    flow_grid_number_cells_vx,    &
                                    flow_grid_number_cells_vy

    namelist /simple_pic_params/    NUM_PARTICLES

    print *, "AA0"

    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_vp4d_par_cart() failed to open file ', filename
       STOP
    end if

    print *, "AA1"
    read(input_file, sim_params)
    read(input_file, elec_params)
    read(input_file, poisson_grid_params)
    if( UseLtPicScheme )then
        read(input_file, lt_pic_remap_params)
        read(input_file, lt_pic_deposition_params)
        read(input_file, lt_pic_markers_params)
    else
        read(input_file, simple_pic_params)
    end if
    close(input_file)

    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)
    
    print*, 'sim%world_size=',sim%world_size, ' sim%my_rank=', sim%my_rank

    XMAX = (2._f64*sll_pi/KX_LANDAU)
    sim%use_lt_pic_scheme = UseLtPicScheme
    print *, "AAA"

    sim%thermal_speed_ions = THERM_SPEED
    sim%total_density = 1._f64 * (XMAX - XMIN) * (YMAX - YMIN)

    sim%dt = dt
    sim%num_iterations = number_iterations
    sim%plot_period = plot_period
    sim%remap_period = remap_period
    sim%elec_params = (/KX_LANDAU, ALPHA, er, psi, omega_r, omega_i /)
    
    sim%mesh_2d =>  new_cartesian_mesh_2d( NC_X, NC_Y, &
                                           XMIN, XMAX, YMIN, YMAX )

    sim%poisson => new_poisson_2d_fft_solver( sim%mesh_2d%eta1_min,    &
                                              sim%mesh_2d%eta1_max,    &
                                              sim%mesh_2d%num_cells1,  &
                                              sim%mesh_2d%eta2_min,    &
                                              sim%mesh_2d%eta2_max,    &
                                              sim%mesh_2d%num_cells2   )

    print *, "[simulation] constructing the particle group..."

    ! construct [[particle_group]]
    if( sim%use_lt_pic_scheme )then


      ! construct [[bsl_lt_pic_particle_group]]
      particle_group_id = 1
      remapping_sparse_grid_max_levels(1) = remapping_sparse_grid_max_level_x
      remapping_sparse_grid_max_levels(2) = remapping_sparse_grid_max_level_y
      remapping_sparse_grid_max_levels(3) = remapping_sparse_grid_max_level_vx
      remapping_sparse_grid_max_levels(4) = remapping_sparse_grid_max_level_vy

      bsl_lt_pic_particle_group => sll_bsl_lt_pic_4d_group_new( &
        SPECIES_CHARGE,                             &
        SPECIES_MASS,                               &
        particle_group_id,                          &
        DOMAIN_IS_X_PERIODIC,                       &
        DOMAIN_IS_Y_PERIODIC,                       &
        remap_f_type,                               &
        remap_degree,                               &
        remapping_grid_vx_min,                          &
        remapping_grid_vx_max,                          &
        remapping_grid_vy_min,                          &
        remapping_grid_vy_max,                          &
        remapping_cart_grid_number_cells_x,                  &   ! for splines
        remapping_cart_grid_number_cells_y,                  &   ! for splines
        remapping_cart_grid_number_cells_vx,                 &   ! for splines
        remapping_cart_grid_number_cells_vy,                 &   ! for splines
        remapping_sparse_grid_max_levels,                     &   ! for the sparse grid: for now, same level in each dimension
        number_deposition_particles,                &
        number_markers_x,                           &
        number_markers_y,                           &
        number_markers_vx,                          &
        number_markers_vy,                          &
        flow_grid_number_cells_x,                   &
        flow_grid_number_cells_y,                   &
        flow_grid_number_cells_vx,                  &
        flow_grid_number_cells_vy,                  &
        sim%mesh_2d )

      ! here, [[particle_group]] will contain a reference to a bsl_lt_pic group
      sim%particle_group => bsl_lt_pic_particle_group
      sim%number_deposition_particles = number_deposition_particles


    else

      ! construct [[simple_pic_particle_group]]
      particle_group_id = 0

      simple_pic_particle_group => sll_simple_pic_4d_group_new(     &
            SPECIES_CHARGE,                                             &
            SPECIES_MASS,                                               &
            particle_group_id,                                          &
            DOMAIN_IS_X_PERIODIC,                                       &
            DOMAIN_IS_Y_PERIODIC,                                       &
            NUM_PARTICLES,                                              &
            sim%mesh_2d )

      ! here [[particle_group]] will contain a reference to a simple pic group
      sim%particle_group => simple_pic_particle_group
      sim%number_deposition_particles = NUM_PARTICLES

    end if

    sim%number_particles = sim%particle_group%number_particles    ! with bsl_lt_pic this is actually the number of markers

    sim%domain_is_x_periodic = DOMAIN_IS_X_PERIODIC
    sim%domain_is_y_periodic = DOMAIN_IS_Y_PERIODIC

    ! initialization of the particle group
    ! todo: should we write a structure for the initialization?
    call sim%particle_group%set_landau_parameters( sim%thermal_speed_ions, ALPHA, KX_LANDAU )

    rand_seed_size = 10
    print *, "rand_seed_size = ", rand_seed_size
    print *, "sim%my_rank = ", sim%my_rank
    call random_seed (SIZE=rand_seed_size)
    print *, "AA BB"

    SLL_ALLOCATE( rand_seed(1:rand_seed_size), ierr )
    do j=1, rand_seed_size
      rand_seed(j) = (-1)**j*(100 + 15*j)*(2*sim%my_rank + 1)
    end do

    print *, "BBB"

    initial_density_identifier = 0     ! for the moment we only use one density (landau)
    call sim%particle_group%initializer( initial_density_identifier, rand_seed, sim%my_rank, sim%world_size )
    SLL_DEALLOCATE_ARRAY(rand_seed, ierr)

    print *, "CCC"

    sim%n_threads = 1
    print*, 'number of threads is ', sim%n_threads

   SLL_ALLOCATE(sim%q_accumulator_ptr(1:sim%n_threads), ierr)

   thread_id = 0
   sim%q_accumulator_ptr(thread_id+1)%q => new_charge_accumulator_2d( sim%mesh_2d )
   sim%E_accumulator => new_field_accumulator_2d( sim%mesh_2d )

   !! -- --  First charge deposition [begin]  -- --

   charge_accumulator => sim%q_accumulator_ptr(thread_id+1)%q

   target_total_charge = SPECIES_CHARGE * 1._f64 * (XMAX - XMIN) * (YMAX - YMIN)
   enforce_total_charge = .true.
   call sim%particle_group%deposit_charge_2d( charge_accumulator, target_total_charge, enforce_total_charge )

   !! -- --  First charge deposition [end]  -- --


  end subroutine init_4d_generic_pic_cartesian


  ! <<run_4d_generic_pic_cartesian>>
  
  !> Run the Vlasov-Poisson simulation. Redefines sll_simulation_base_class::run <!--
  !> [[selalib:src/simulations/simulation_base_class.F90::run]] -->
  !!
  !! \note 1: this is a skeleton-in-progress: some routines are not implemented, some variables are not needed
  !!
  !! \note 2: use of cubic spline particles (routines with _CS) is disabled for now
  !!
  !! \todo 1: write and check the remapping step with adequate frequency
  !!
  !!
  !! @author this version written by MCP, ALH

  subroutine run_4d_generic_pic_cartesian( sim )
    class(sll_simulation_4d_vp_generic_pic_cartesian), intent(inout)  :: sim
    sll_int32  :: ierr, it, jj, counter
    sll_int32  :: i, j, k
    sll_real64 :: tmp1, tmp2, tmp3, tmp4
    sll_real64 :: tmp5, tmp6, temp
    sll_real32 :: ttmp(1:4,1:2), ttmp1(1:4,1:2), ttmp2(1:4,1:2)
    sll_real64 :: valeur, val2
    sll_real64, dimension(:,:), pointer :: phi
    sll_int32  :: ncx, ncy, ic_x,ic_y
    sll_int32  :: ic_x1,ic_y1
    sll_real32 :: off_x, off_y,off_x1,off_y1
    sll_real64 :: xmin, ymin
    sll_real64 :: rdx, rdy
    ! sll_int32  :: gi ! counter index for guard list
    sll_real64 :: Ex, Ey, Ex1, Ey1, Ex_CS, Ey_CS
    sll_real64 :: dt_q_over_m ! dt * qoverm
    sll_real64 :: x, x1  ! for global position
    sll_real64 :: y, y1  ! for global position
    sll_real64 :: dt
    sll_real64 :: pp_x,pp_y,pp_vx, pp_vy
    sll_real32 :: pp_dx, pp_dy
    sll_int32  :: pp_icell, pp_icell_x, pp_icell_y
    type(field_accumulator_cell), dimension(:), pointer :: accumE
    ! type(sll_generic_pic_4d_particle_guard), dimension(:), pointer :: p_guard
    sll_real64, dimension(:,:), allocatable  ::  diag_energy
    sll_real64, dimension(:),   allocatable  ::  diag_TOTmoment
    sll_real64, dimension(:),   allocatable  ::  diag_TOTenergy
    sll_real64, dimension(:,:), allocatable :: diag_AccMem! a memory buffer
    sll_real64 :: t2, t3
    sll_real64, dimension(:), allocatable :: rho1d_send
    sll_real64, dimension(:), allocatable :: rho1d_receive
    sll_real64   :: t_init, t_fin, time
    sll_int32 :: save_nb
    sll_int32 :: thread_id
    sll_int32 :: n_threads
    sll_int32 :: plot_np_x    !< nb of points in the x  plotting grid for a (x,vx) plot
    sll_int32 :: plot_np_y    !< nb of points in the y  plotting grid for a (x,vx) plot
    sll_int32 :: plot_np_vx   !< nb of points in the vx plotting grid for a (x,vx) plot
    sll_int32 :: plot_np_vy   !< nb of points in the vy plotting grid for a (x,vx) plot

    type(sll_charge_accumulator_2d),    pointer :: charge_accumulator
    sll_int32 :: sort_nb
    sll_real64 :: some_val, une_cst
    sll_real64 :: val_lee, exval_ee
    sll_real64 :: tot_ee, val_ee
    sll_real64 :: omega_i, omega_r, psi
    sll_real64 :: bors
    sll_real64 :: coords(3)
    sll_real64 :: target_total_charge
    logical    :: enforce_total_charge


    ! Timings and statistics
    sll_real64 :: deposit_time
    sll_real64 :: loop_time
    type(sll_time_mark) :: deposit_time_mark, loop_time_mark
    
    ! ------------------------------
    ncx = sim%mesh_2d%num_cells1
    ncy = sim%mesh_2d%num_cells2
    n_threads = sim%n_threads
    thread_id = 0
    save_nb = sim%num_iterations/10

    SLL_ALLOCATE( rho1d_send(1:(ncx+1)*(ncy+1)),    ierr)
    SLL_ALLOCATE( rho1d_receive(1:(ncx+1)*(ncy+1)), ierr)

    SLL_ALLOCATE(sim%rho(ncx+1,ncy+1),ierr)
    SLL_ALLOCATE( sim%E1(1:ncx+1,1:ncy+1), ierr )
    SLL_ALLOCATE( sim%E2(1:ncx+1,1:ncy+1), ierr )
    SLL_ALLOCATE(phi(1:ncx+1, 1:ncy+1), ierr)
    SLL_ALLOCATE(diag_energy(1:save_nb, 1:5), ierr)
    SLL_ALLOCATE(diag_TOTmoment(1:save_nb), ierr)
    SLL_ALLOCATE(diag_TOTenergy(0:sim%num_iterations-1), ierr)
    SLL_ALLOCATE(diag_AccMem(0:sim%num_iterations-1, 1:2), ierr)

    sort_nb = 10
    dt = sim%dt

    ! <<dt_q_over_m>>
    dt_q_over_m = dt * sim%particle_group%species%q_over_m()
    print*,  "dt_q_over_m = ", dt_q_over_m
    
    xmin = sim%mesh_2d%eta1_min
    ymin = sim%mesh_2d%eta2_min
    rdx = 1._f64/sim%mesh_2d%delta_eta1
    rdy = 1._f64/sim%mesh_2d%delta_eta2

    print*,  "aaaaaa"

    !  ----------------------------------------------------------------------------------------------------
    !> ## Time loop initialisation
    !>
    !> The time loop is prepared by computing the \f$E\f$ field
    !>
    !> - starting from
    !>   * \f$(x,y)^0_k\f$, \f$(v_x, v_y)^0_k\f$ stored in the particle group
    !>     sll_simulation_4d_vp_generic_pic_cartesian::particle_group
    !>
    !> - we end with
    !>   * \f$E^0\f$ stored in sim\%E1 (\f$E^0_x\f$), sim\%E2 (\f$E^0_y\f$)
    !>   * \f$(x,y)^0_k, (v_x, v_y)^{-1/2}_k\f$ stored in the the particle group
    !>     sll_simulation_4d_vp_generic_pic_cartesian::particle_group
    !  ----------------------------------------------------------------------------------------------------

    !! -- --  ?? [begin]  -- --

    accumE => sim%E_accumulator%e_acc
    call sum_accumulators( sim%q_accumulator_ptr, n_threads, ncx*ncy )

    call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator_ptr(1)%q, sim%rho )     ! this name not clear enough

    print*,  "aaa bb"

    !! -- --  ?? [end]  -- --


    !! -- --  MPI communications of rho [begin]  -- --

    do j = 1, ncy+1
       do i = 1, ncx+1
          rho1d_send(i+(j-1)*(ncx+1)) = sim%rho(i, j)
          rho1d_receive(i+(j-1)*(ncx+1)) = 0._f64
       enddo
    enddo
    
    call sll_collective_allreduce( sll_world_collective, rho1d_send, (ncx+1)*(ncy+1), &
         MPI_SUM, rho1d_receive   )

    print*,  "aaaa cc"

    do j = 1, ncy+1
       do i = 1, ncx+1
          sim%rho(i, j) = rho1d_receive(i+(j-1)*(ncx+1))
       enddo
    enddo

    !! -- --  MPI communications of rho [end]  -- --

    print*,  "aa  dd"

    print*,  "aa  dd 1 ", xmin
    print*,  "aa  dd 2 ", sim%mesh_2d%eta1_max
    print*,  "aa  dd 3 ", ncx+1
    print*,  "aa  dd 4 ", ymin
    print*,  "aa  dd 5 ", sim%mesh_2d%eta2_max
    print*, sim%my_rank
    print*, sim%rho
    print*, size(sim%rho)


    if (sim%my_rank == 0) then
       it = 0

       ! [[selalib:src/io/file_io/sll_m_gnuplot.F90::sll_gnuplot_corect_2d]] uses a plot number starting from 1

       ! <<rho_init_standPUSH>> This will also generate the corresponding gnuplot script
       call sll_gnuplot_2d(xmin, sim%mesh_2d%eta1_max, ncx+1, ymin, &
            sim%mesh_2d%eta2_max, ncy+1,                            &
            sim%rho, 'rho_init_standPUSH', it+1, ierr )
    endif

    !> The initial field \f$E^0\f$ is obtained with a call to the Poisson solver. Note that here sim\%rho has the proper
    !> sign (hence there is no need to multiply it by an additional physical constant). The resulting field \f$E^0_x\f$
    !> is stored in sim\%E1, and \f$E^0_y\f$ in sim\%E2.

    print*,  "aa  ee"

    call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, sim%rho )

    print*,  "aa  ff"

    ! <<Ex_Ey_output>> using the [[selalib:src/io/file_io/sll_m_gnuplot.F90::sll_gnuplot_2d]] interface and most probably
    ! the [[selalib:src/io/file_io/sll_m_gnuplot.F90::sll_gnuplot_corect_2d]] implementation.
    
    if (sim%my_rank == 0) then
       it = 0
       
       ! [[selalib:src/io/file_io/sll_m_gnuplot.F90::sll_gnuplot_corect_2d]] uses a plot number starting from 1

       print *, "writing Ex, Ey in gnuplot format for iteration # it = ", it," # plot = ", it+1
       call sll_gnuplot_2d(xmin, sim%mesh_2d%eta1_max, ncx+1, ymin, &
            sim%mesh_2d%eta2_max, ncy+1,                            &
            sim%E1, 'Ex', it+1, ierr )

       call sll_gnuplot_2d(xmin, sim%mesh_2d%eta1_max, ncx+1, ymin, &
            sim%mesh_2d%eta2_max, ncy+1,                            &
            sim%E2, 'Ey', it+1, ierr )

    endif


    !! -- --  diagnostics: compute energy [begin]  -- --

    call electric_energy( tot_ee, sim%E1, sim%E2, ncx, &
         ncy, sim%mesh_2d%eta1_max, sim%mesh_2d%eta2_max )
    bors = 0.0_f64

    do k = 1, sim%number_particles

       !> This simulation does not have access to the particles (because they may be of different incompatible types
       !> like "ltpic" or "simple") so we use the standard interface defined in
       !> sll_m_remapped_pic_base::sll_c_remapped_particle_group

       coords = sim%particle_group%get_v(k)
       bors = bors + coords(1)**2 + coords(2)**2
    enddo
    diag_TOTenergy(0) = bors * 0.5_f64*(sim%mesh_2d%eta1_max - sim%mesh_2d%eta1_min)  &
         * (sim%mesh_2d%eta2_max - sim%mesh_2d%eta2_min)/( sim%world_size*sim%number_particles)  &
         + tot_ee * 0.5_f64

    !! -- --  diagnostics [end]  -- --


    !! -- --  half v-push (computing v^0 -> v^{-1/2})  [begin]  -- --

    call reset_field_accumulator_to_zero( sim%E_accumulator )
    call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )

    do k = 1, sim%number_particles

      ! particle position
      coords=sim%particle_group%get_x(k)
      pp_x = coords(1)
      pp_y = coords(2)

      ! particle speed
      coords=sim%particle_group%get_v(k)
      pp_vx = coords(1)
      pp_vy = coords(2)

      !> compute cell dx and dy again locally (although we could extract this info for some types of particles)
      !> because the base PIC class does not give access to this type of info (because some particles may not have
      !> it)

      call global_to_cell_offset_extended(pp_x, pp_y, sim%mesh_2d, pp_icell_x, pp_icell_y, pp_dx, pp_dy)
      call get_poisson_cell_index(sim%mesh_2d, pp_icell_x, pp_icell_y, pp_icell)

      SLL_INTERPOLATE_FIELD_IN_CELL(Ex,Ey, accumE, pp_dx, pp_dy, tmp5,tmp6, pp_icell)


      ! Set particle speed [[dt_q_over_m]]
      coords(1) = pp_vx - 0.5_f64 * dt_q_over_m * Ex
      coords(2) = pp_vy - 0.5_f64 * dt_q_over_m * Ey
      coords(3) = 0.0_f64
      call sim%particle_group%set_v(k, coords)
    enddo

    ! -- --  half v-push  [end]  -- --

    ! this means (pi / KX_LANDAU) * (4._f64 * ALPHA * er)**2
    une_cst = (sll_pi / sim%elec_params(1)) * (4._f64 * sim%elec_params(2) * sim%elec_params(3))**2
    omega_i = sim%elec_params(6)
    omega_r = sim%elec_params(5)
    psi = sim%elec_params(4)

    ! <<logE_standPush>>
    if (sim%my_rank ==0) open(65,file='logE_standPush.dat')
    !#ifdef _OPENMP
    !    t2 = omp_get_wtime() !   call sll_set_time_mark(t2)
    !#endif

    !  ----------------------------------------------------------------------------------------------------
    !> ## Content of the time loop
    !>          - starting from
    !>              - \f$E^n\f$ stored in sim\%E1, sim\%E2
    !>              - \f$(x,y)^n_k\f$, \f$(v_x, v_y)^{n-1/2}_k\f$ stored in the particle group (with \f$n = it\f$)
    !>          - we end with
    !>              - \f$E^{n+1}\f$ stored in sim\%E1, sim\%E2
    !>              - \f$(x,y)^{n+1}_k\f$, \f$(v_x, v_y)^{n+1/2}_k\f$ stored in the particle group
    !  ----------------------------------------------------------------------------------------------------

    ! Time statistics
    call sll_set_time_mark(loop_time_mark)
    deposit_time = 0._f64

    ! First snapshot at time 0 [[particles_snapshot]]
    ! AAA_ALH_TODO call particles_snapshot(0.0_8,sim)
    
    do it = 0, sim%num_iterations-1

       print *, "BEGIN one loop in time, it+1 = ", it+1, " / ", sim%num_iterations
       !! -- --  <<diagnostics>> (computing energy) [begin]  -- --

       if (sim%my_rank == 0) then
          exval_ee = une_cst * exp(2._f64 * omega_i * real(it,f64) * sim%dt)                          &
                             * ( 0.5_f64 + 0.5_f64 * cos(2._f64 * (omega_r * real(it,f64) * sim%dt - psi)) )
          exval_ee = max(exval_ee, 1e-30_f64)    ! to avoid taking the log of 0

          ! [[normL2_field_Ex]]
          call normL2_field_Ex ( val_lee, val_ee, ncx, ncy,                         &
                                 sim%E1,                                            &
                                 sim%mesh_2d%delta_eta1, sim%mesh_2d%delta_eta2 )
          counter = 1 + mod(it,save_nb)
          diag_energy(counter,:) = (/ it*sim%dt, val_lee, log(sqrt(exval_ee)), &
               val_ee, exval_ee /)

          if ( mod(it+1,save_nb)==0 ) then
             do jj=1,save_nb
                write(65,*) diag_energy(jj,:)!, diag_TOTmoment(jj)!, diag_TOTenergy(jj)
             enddo
          endif
       endif

       ! -- --  diagnostics [end]  -- --

       !  ------------------------------------------------------------------
       !  ------
       !  ------  PUSH PARTICLES [begin]
       !  ------
       !  ------------------------------------------------------------------

       !> ### Particles are pushed as follows

      some_val = 0.0_f64

      thread_id = 0
      call reset_charge_accumulator_to_zero ( sim%q_accumulator_ptr(thread_id+1)%q )
      charge_accumulator => sim%q_accumulator_ptr(thread_id+1)%q

      !! -- --  particle loop -- --

      do k = 1, sim%number_particles

         ! -- --  v-push [begin]  -- --
         !> #### v-push
         !> using \f$E^n\f$  we compute  \f$v^{n-1/2} \rightarrow v^{n+1/2}\f$

         ! particle position
         coords=sim%particle_group%get_x(k)
         pp_x = coords(1)
         pp_y = coords(2)

         ! particle speed
         coords=sim%particle_group%get_v(k)
         pp_vx = coords(1)
         pp_vy = coords(2)

         !> compute cell dx and dy again locally (although we could extract this info for some types of particles)
         !> because the base PIC class does not give access to this type of info (because some particles may not
         !> have it).

         call global_to_cell_offset_extended(pp_x, pp_y, sim%mesh_2d, pp_icell_x, pp_icell_y, pp_dx, pp_dy)
         call get_poisson_cell_index(sim%mesh_2d, pp_icell_x, pp_icell_y, pp_icell)

         SLL_INTERPOLATE_FIELD_IN_CELL(Ex,Ey, accumE, pp_dx, pp_dy, tmp5,tmp6, pp_icell)

         ! Set particle speed [[dt_q_over_m]]
         coords(1) = pp_vx + dt_q_over_m * Ex
         coords(2) = pp_vy + dt_q_over_m * Ey
         coords(3) = 0.0_f64
         call sim%particle_group%set_v(k, coords)

         !> remember new speed for x-push
         pp_vx=coords(1)
         pp_vy=coords(2)

         !! -- --  v-push [end]  -- --

         !! -- --  x-push [begin]  -- --
         !> #### x-push
         !> using \f$v^{n+1/2}\f$  we compute \f$x^n \rightarrow x^{n+1}\f$

         coords(1) = pp_x + dt * pp_vx
         coords(2) = pp_y + dt * pp_vy
         coords(3) = 0.0_f64

         !! -- --  x-push [end]  -- --

         if( .not. x_is_in_domain_2d( coords(1), coords(2), sim%mesh_2d,            &
                                       sim%domain_is_x_periodic,                    &
                                       sim%domain_is_y_periodic )                   &
               ) then
            !! -- -- put outside particles back in domain
            ! [[selalib:src/particle_methods/pic_remapped/remapped_pic_utilities/sll_m_remapped_pic_utilities.F90::apply_periodic_bc_on_cartesian_mesh_2d]]
            call apply_periodic_bc_on_cartesian_mesh_2d( sim%mesh_2d, coords(1), coords(2))
         end if

         call sim%particle_group%set_x(k, coords)

      end do

      !> ### Charge deposit
      
      print *, "deposit charge begin"

      call sll_set_time_mark(deposit_time_mark)

      ! [[file:~/selalib/src/particle_methods/sll_pic_base.F90::deposit_charge_2d]]
      target_total_charge = 0._f64
      enforce_total_charge = .false.   ! todo: try with true
      charge_accumulator => sim%q_accumulator_ptr(thread_id+1)%q
      call sim%particle_group%deposit_charge_2d( charge_accumulator, target_total_charge, enforce_total_charge )

      deposit_time=deposit_time+sll_time_elapsed_since(deposit_time_mark)

      print *, "deposit charge end"

      !  ------------------------------------------------------------------
      !  ------
      !  ------  PUSH PARTICLES [end]
      !  ------
      !  ------------------------------------------------------------------

      !! -- --  parallel communications for rho ?? [begin]  -- --

      call sum_accumulators( sim%q_accumulator_ptr, n_threads, ncx*ncy )
      call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator_ptr(1)%q, sim%rho )     ! this name not clear enough
      do j = 1, ncy+1
         do i = 1, ncx+1
            rho1d_send(i+(j-1)*(ncx+1)) = sim%rho(i, j)
            rho1d_receive(i+(j-1)*(ncx+1)) = 0._f64
         end do
      end do

      call sll_collective_allreduce( sll_world_collective, rho1d_send, (ncx+1)*(ncy+1), &
            MPI_SUM, rho1d_receive   )

      do j = 1, ncy+1
         do i = 1, ncx+1
            sim%rho(i, j) = rho1d_receive(i+(j-1)*(ncx+1))
         end do
      end do

      !! -- --  parallel communications for rho ?? [end]  -- --

      !> ### Poisson solver (computing \f$E^{n+1}\f$)

      !> In the time loop, the field \f$E^{n+1}\f$ is obtained with a call to the Poisson solver.  Again, sim\%rho has
      !> the proper sign so that we do not need to multiply it by an additional physical constant.

      call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, sim%rho )

      !! -- --  diagnostics (plotting) [begin]  -- --

      if (sim%my_rank == 0 .and. mod(it, sim%plot_period)==0 ) then

        print *, "plotting f slice in gnuplot format for iteration # it = ", it, " / ", sim%num_iterations
        plot_np_x  = 100
        plot_np_y  = 6
        plot_np_vx = 30
        plot_np_vy = 5

        ! <<f_slice>> base class definition of visualize_f_slice_x_vx:
        !   [[selalib:src/particle_methods/pic_remapped/sll_m_remapped_pic_base.F90::visualize_f_slice_x_vx]]
        ! specialized in:
        ! - [[selalib:src/particle_methods/pic_remapped/bsl_lt_pic/sll_m_bsl_lt_pic_4d_group.F90::bsl_lt_pic_4d_visualize_f_slice_x_vx]]
        ! - [[selalib:src/particle_methods/pic_remapped/simple_pic/sll_m_simple_pic_4d_group.F90::simple_pic_4d_visualize_f_slice_x_vx]]

        call sim%particle_group%visualize_f_slice_x_vx("f_slice", plot_np_x, plot_np_y, plot_np_vx, plot_np_vy, it+1)

      end if

      if (sim%use_lt_pic_scheme .and. sim%my_rank == 0 .and. mod(it+1, sim%remap_period)==0 ) then
        ! note: condition on rank == 0 needs to be revised for the actual parallel version...

        print *, "remapping f..."
        call sim%particle_group%remap()
      end if

       if (sim%my_rank == 0 .and. mod(it+1, sim%plot_period)==0 ) then

          print *, "writing Ex, Ey  in gnuplot format for iteration # it = ", it, " / ", sim%num_iterations, &
               " # plot = ",it+1
            call sll_gnuplot_2d(xmin, sim%mesh_2d%eta1_max, ncx+1, ymin,            &
                                sim%mesh_2d%eta2_max, ncy+1,                        &
                                sim%E1, 'Ex', it+1, ierr )

            call sll_gnuplot_2d(xmin, sim%mesh_2d%eta1_max, ncx+1, ymin,            &
                                sim%mesh_2d%eta2_max, ncy+1,                        &
                                sim%E2, 'Ey', it+1, ierr )
            print *, "done."

       end if

       !! -- --  diagnostics (plotting E ) [end]  -- --


       !! -- --  diagnostics (computing energy) [begin]  -- --
        print *, "diag energy"

        call reset_field_accumulator_to_zero( sim%E_accumulator )
        call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
        some_val = 0.0_f64

        do k = 1, sim%number_particles

             ! particle position
             coords=sim%particle_group%get_x(k)
             pp_x = coords(1)
             pp_y = coords(2)

             !> compute cell dx and dy again locally (although we could extract this info for some types of particles)
             !> because the base PIC class does not give access to this type of info (because some particles may not
             !> have it).  todo: we may want to optimize this out for speed

             call global_to_cell_offset_extended(pp_x, pp_y, sim%mesh_2d, pp_icell_x, pp_icell_y, pp_dx, pp_dy)
             call get_poisson_cell_index(sim%mesh_2d, pp_icell_x, pp_icell_y, pp_icell)

             SLL_INTERPOLATE_FIELD_IN_CELL(Ex,Ey, accumE, pp_dx, pp_dy, tmp5,tmp6, pp_icell)

             ! call global_to_cell_offset(pp_x,pp_y,sim%mesh_2d,pp_icell,pp_dx,pp_dy)

             ! SLL_INTERPOLATE_FIELD_IN_CELL(Ex,Ey,accumE,pp_dx,pp_dy,tmp5,tmp6,pp_icell)

             ! particle speed
             coords=sim%particle_group%get_v(k)
             pp_vx = coords(1)
             pp_vy = coords(2)

            ! todo: check the energy diagnostics ([[dt_q_over_m]]??)
             some_val = some_val &
                    + (pp_vx + 0.5_f64 * dt_q_over_m * Ex)**2 &
                    + (pp_vy + 0.5_f64 * dt_q_over_m * Ey)**2
        enddo

        call electric_energy( tot_ee, sim%E1, sim%E2, ncx, &
           ncy, sim%mesh_2d%eta1_max, sim%mesh_2d%eta2_max )
        diag_TOTenergy(it) = some_val* 0.5_f64 *(sim%mesh_2d%eta1_max - sim%mesh_2d%eta1_min) * &
           (sim%mesh_2d%eta2_max - sim%mesh_2d%eta2_min)/( sim%world_size*sim%number_particles)  &
           + tot_ee * 0.5_f64

       !! -- --  diagnostics [end]  -- --

        ! Another snapshot after each iteration [[particles_snapshot]]
        ! AAA_ALH_TODO call particles_snapshot(it*dt,sim)
            
    print *, "end one loop in time"
    enddo

    !  ----------------------------------------------------------------------------------------------------
    !  ------
    !  ------  TIME LOOP [end]
    !  ------
    !  ----------------------------------------------------------------------------------------------------

    !#ifdef _OPENMP
    !    time = omp_get_wtime()!! time = sll_time_elapsed_since(t2)
    !
    !    if (sim%my_rank ==0) then
    !       close(65)
    !       open(93,file='time_parts_sec_omp.dat',position='append')
    !       write(93,*) '# Nb of threads  ||  time (sec)  ||  average pushes/sec'
    !       write(93,*) sim%n_threads, time-t2, int(sim%num_iterations,i64)*int(sim%number_particles,i64)/(time-t2)
    !       close(93)
    !       open(65,file='Energie_totale_standPush.dat')
    !       do it=0,sim%num_iterations-1
    !          write(65,*) it*dt,  diag_TOTenergy(it), (diag_TOTenergy(it)-diag_TOTenergy(0))/diag_TOTenergy(0)
    !       enddo
    !       close(65)
    !    endif
    !#endif

    ! ALH - Time statistics. The number of deposits is computed separately to take the number of LTP BSL virtual
    ! particles into account.
    
    loop_time=sll_time_elapsed_since(loop_time_mark)
    write(*,'(A,ES8.2,A)') 'sim stats: ',1 / loop_time * sim%num_iterations * sim%number_particles,' pushes/sec '
    if(sim%use_lt_pic_scheme)then
       write(*,'(A,ES8.2,A)') 'lt_pic stats: ',                                     &
            1 / deposit_time * sim%num_iterations * sim%number_deposition_particles,    &
            ' deposits/sec'
    end if

    SLL_DEALLOCATE(sim%rho,   ierr)
    SLL_DEALLOCATE(sim%E1,    ierr)
    SLL_DEALLOCATE(sim%E2,    ierr)
    SLL_DEALLOCATE(phi, ierr)
    SLL_DEALLOCATE_ARRAY(rho1d_send, ierr)
    SLL_DEALLOCATE_ARRAY(rho1d_receive, ierr)
    SLL_DEALLOCATE_ARRAY(diag_energy, ierr)
    SLL_DEALLOCATE_ARRAY(diag_TOTenergy, ierr)
    SLL_DEALLOCATE_ARRAY(diag_TOTmoment, ierr)
    SLL_DEALLOCATE_ARRAY(diag_AccMem, ierr)

  end subroutine run_4d_generic_pic_cartesian


  subroutine delete_4d_generic_pic_cartesian( sim )
    type(sll_simulation_4d_vp_generic_pic_cartesian) :: sim
  end subroutine delete_4d_generic_pic_cartesian


!  function in_bounds( x, y, mesh ) result(res)
!    logical :: res
!    sll_real64, intent(in) :: x
!    sll_real64, intent(in) :: y
!    type(sll_cartesian_mesh_2d), pointer :: mesh
!
!    res = (x >= mesh%eta1_min) .and. (x <= mesh%eta1_max) .and. &
!          (y >= mesh%eta2_min) .and. (y <= mesh%eta2_max)
!  end function in_bounds

!  subroutine apply_periodic_bc( mesh, x, y )
!    type(sll_cartesian_mesh_2d), pointer :: mesh
!    sll_real64, intent(inout) :: x
!    sll_real64, intent(inout) :: y
!
!    x = modulo(x,mesh%eta1_max - mesh%eta1_min)! Ca marche quand xmin=ymin=0
!    y = modulo(y,mesh%eta2_max - mesh%eta2_min)! Sinon, il faut modulo(...) + xmin/mesh%delta_x
!  end subroutine apply_periodic_bc

  ! <<normL2_field_Ex>>
  subroutine normL2_field_Ex (lee,ee,nx,ny,e,dx,dy)
    sll_real64, intent(out) :: lee, ee
    sll_real64, intent(in) :: dx,dy
    sll_int32, intent(in) :: nx,ny
    sll_real64, dimension(1:nx+1,1:ny+1),intent(in) :: e
    sll_int32 :: i,j
    
    lee = 0._f64
    do j=1,ny
       do i=1,nx
          lee = lee + e(i,j)**2 !*e(i,j)
       enddo
    enddo
    lee = lee*dx*dy
    ee = lee
    lee = max(lee, 1e-30_f64)    ! to avoid taking the log of 0
    lee = log(lee)*0.5_f64
  end subroutine normL2_field_Ex
  
    subroutine electric_energy(ee,ex,ey,nx,ny,dx,dy)
    ! ee  = the electric energy of (Ex,Ey)
    sll_real64, intent(out) :: ee
    sll_real64, intent(in)  :: dx,dy
    sll_int32, intent(in)   :: nx,ny
    sll_real64, dimension(1:nx+1,1:ny+1),intent(in) :: ex, ey
    sll_int32 :: i, j

    ee = 0._f64
    do j=1,ny
       do i=1,nx
          ee = ee + ex(i,j)**2 + ey(i,j)**2
       enddo
    enddo
    ee = ee * dx * dy    
  end subroutine electric_energy

end module sll_m_sim_pic_vp_2d2v_cart_remapped
