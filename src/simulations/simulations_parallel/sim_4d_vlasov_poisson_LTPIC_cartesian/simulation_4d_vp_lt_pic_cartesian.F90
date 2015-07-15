!> @details
!> The base for this development was \ref sll_pic_simulation_4d_cartesian_module

! [[selalib:src/simulations/simulations_parallel/sim_4d_vlasov_poisson_PIC_cartesian/simulation_4d_vp_pic_cartesian.F90::sll_pic_simulation_4d_cartesian_module]]

module sll_simulation_4d_vp_lt_pic_cartesian_module

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_utilities.h"
#include "sll_accumulators.h" 
#include "particle_representation.h"

  use sll_constants
  use sll_simulation_base
  use sll_cartesian_meshes
  use sll_timer
  use sll_particle_group_4d_module

  use sll_lt_pic_4d_group_module
  use sll_lt_pic_4d_init

  use sll_particle_initializers_4d
  use sll_particle_sort_module
  use sll_charge_to_density_module
  use sll_pic_utilities
  use sll_module_poisson_2d_fft
  use sll_module_poisson_2d_base
  use sll_representation_conversion_module
  use sll_gnuplot
  use sll_collective 
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none

  type, extends(sll_simulation_base_class) :: sll_simulation_4d_vp_lt_pic_cartesian
    ! Physics/numerical parameters
    sll_real64 :: dt
    sll_int32  :: num_iterations
    sll_int32  :: plot_period
    sll_int32  :: remap_period
    sll_real64 :: thermal_speed_ions
    sll_int32  :: ions_number
    sll_int32  :: virtual_particle_number
!    sll_int32  :: guard_size
!    sll_int32  :: array_size
    sll_real64, dimension(1:6) :: elec_params
    !!!!   note (march 26): to start with we use only lt_pic particles. When the simulation runs fine, try using a common type
    type(sll_lt_pic_4d_group),  pointer :: part_group
    !     type(sll_particle_group_4d),  pointer :: part_group           !! pic case
    type(sll_cartesian_mesh_2d),    pointer :: mesh_2d ! [[file:~/mcp/selalib/src/meshes/sll_cartesian_meshes.F90::sll_cartesian_mesh_2d]]
    type(sll_particle_sorter_2d), pointer :: sorter
    type(sll_charge_accumulator_2d_ptr), dimension(:), pointer  :: q_accumulator_ptr      ! called q_accumulator in Sever simulation
    type(electric_field_accumulator), pointer :: E_accumulator
    logical :: use_cubic_splines        ! disabled for now (will raise an error if true)
    logical :: use_lt_pic_scheme        ! if false then use pic scheme
    type(sll_charge_accumulator_2d_CS_ptr), dimension(:), pointer  :: q_accumulator_CS
    type(electric_field_accumulator_CS), pointer :: E_accumulator_CS
    sll_real64, dimension(:,:), pointer :: rho
    sll_int32 :: n_virtual_x_for_deposition, n_virtual_y_for_deposition
    sll_int32 :: n_virtual_vx_for_deposition, n_virtual_vy_for_deposition
    type(poisson_2d_fft_solver), pointer :: poisson
    sll_real64 :: total_density

    sll_real64, dimension(:,:), pointer :: E1, E2
    sll_int32 :: my_rank
    sll_int32 :: world_size
    sll_int32 :: n_threads
  contains
    procedure, pass(sim) :: init_from_file => init_4d_lt_pic_cartesian
    procedure, pass(sim) :: run => run_4d_lt_pic_cartesian
  end type sll_simulation_4d_vp_lt_pic_cartesian

  interface sll_delete
     module procedure delete_4d_lt_pic_cartesian
  end interface sll_delete

!!$  interface initialize
!!$     module procedure initialize_4d_qns_general
!!$  end interface initialize

contains

  subroutine init_4d_lt_pic_cartesian( sim, filename )
    intrinsic :: trim
    class(sll_simulation_4d_vp_lt_pic_cartesian), intent(inout) :: sim
    character(len=*), intent(in)                          :: filename
    sll_int32   :: IO_stat
    sll_int32   :: ierr
    sll_int32   :: j
    sll_real64  :: dt
    sll_int32   :: number_iterations, plot_period, remap_period
    sll_int32   :: NUM_PARTICLES, GUARD_SIZE, PARTICLE_ARRAY_SIZE
    sll_real64  :: THERM_SPEED
    sll_real64  :: QoverM, ALPHA
    logical     :: UseCubicSplines
    logical     :: UseLtPicScheme
    logical     :: DOMAIN_IS_X_PERIODIC, DOMAIN_IS_Y_PERIODIC
    sll_int32   :: NC_X,  NC_Y
    sll_real64  :: XMIN, KX_LANDAU, XMAX, YMIN, YMAX
    sll_int32   :: NUM_PARTS_X, NUM_PARTS_Y, NUM_PARTS_VX, NUM_PARTS_VY
    sll_real64  :: REMAP_GRID_VX_MIN, REMAP_GRID_VX_MAX
    sll_real64  :: REMAP_GRID_VY_MIN, REMAP_GRID_VY_MAX
    sll_int32   :: SPLINE_DEGREE
    sll_int32   :: NVirtual_X_ForDeposition
    sll_int32   :: NVirtual_Y_ForDeposition
    sll_int32   :: NVirtual_VX_ForDeposition
    sll_int32   :: NVirtual_VY_ForDeposition

    logical     :: UseExactF0
    sll_real64  :: er, psi, omega_i, omega_r
    sll_int32, parameter  :: input_file = 99
    sll_int32, dimension(:), allocatable  :: rand_seed
    sll_int32   :: rand_seed_size
    sll_int32  :: thread_id
    type(sll_particle_group_4d), pointer  :: pa_gr

    namelist /sim_params/   NUM_PARTICLES, GUARD_SIZE, &
                            PARTICLE_ARRAY_SIZE, &
                            THERM_SPEED, dt, number_iterations, plot_period, remap_period, &
                            QoverM, ALPHA, UseCubicSplines, UseLtPicScheme
    namelist /grid_dims/    NC_X, NC_Y, XMIN, KX_LANDAU, YMIN, YMAX, &
                            DOMAIN_IS_X_PERIODIC, DOMAIN_IS_Y_PERIODIC
    namelist /lt_pic_params/NUM_PARTS_X,            &
                            NUM_PARTS_Y,            &
                            NUM_PARTS_VX,           &
                            NUM_PARTS_VY,           &
                            REMAP_GRID_VX_MIN,      &
                            REMAP_GRID_VX_MAX,      &
                            REMAP_GRID_VY_MIN,      &
                            REMAP_GRID_VY_MAX,      &
                            NVirtual_X_ForDeposition,   &
                            NVirtual_Y_ForDeposition,   &
                            NVirtual_VX_ForDeposition,  &
                            NVirtual_VY_ForDeposition,  &
                            UseExactF0
    namelist /elec_params/  er, psi, omega_r, omega_i
    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_vp4d_par_cart() failed to open file ', filename
       STOP
    end if
    read(input_file, sim_params)
    read(input_file, grid_dims)
    read(input_file, lt_pic_params)
    read(input_file, elec_params)
    close(input_file)

    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)
    
    print*, 'sim%world_size=',sim%world_size, ' sim%my_rank=', sim%my_rank

    XMAX = (2._f64*sll_pi/KX_LANDAU)
    sim%use_cubic_splines = UseCubicSplines
    sim%use_lt_pic_scheme = UseLtPicScheme
    sim%n_virtual_x_for_deposition = NVirtual_X_ForDeposition
    sim%n_virtual_y_for_deposition = NVirtual_Y_ForDeposition
    sim%n_virtual_vx_for_deposition = NVirtual_VX_ForDeposition
    sim%n_virtual_vy_for_deposition = NVirtual_VY_ForDeposition
    sim%thermal_speed_ions = THERM_SPEED
!    sim%guard_size = GUARD_SIZE
!    sim%array_size = PARTICLE_ARRAY_SIZE
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

    sim%sorter => sll_new_particle_sorter_2d( sim%mesh_2d )

    if( sim%use_lt_pic_scheme )then
        !! initialize the group of PIC particles

        if( sim%use_cubic_splines )then
            SPLINE_DEGREE = 3
            print *, "Error (7634876568576) -- cubic spline particles not implemented yet"
        else
            SPLINE_DEGREE = 1
        end if

!        print *, 'tmp debug 65373654 -- REMAP_GRID_VX_MIN, REMAP_GRID_VX_MAX = ', REMAP_GRID_VX_MIN, REMAP_GRID_VX_MAX
!        print *, 'tmp debug 98536536 -- REMAP_GRID_VY_MIN, REMAP_GRID_VY_MAX = ', REMAP_GRID_VY_MIN, REMAP_GRID_VY_MAX

        sim%part_group => sll_lt_pic_4d_group_new(    &
            SPLINE_DEGREE,                            &
            NUM_PARTS_X,                              &
            NUM_PARTS_Y,                              &
            NUM_PARTS_VX,                             &
            NUM_PARTS_VY,                             &
            REMAP_GRID_VX_MIN,                        &
            REMAP_GRID_VX_MAX,                        &
            REMAP_GRID_VY_MIN,                        &
            REMAP_GRID_VY_MAX,                        &
            QoverM,                                   &
            DOMAIN_IS_X_PERIODIC,                     &
            DOMAIN_IS_Y_PERIODIC,                     &
            sim%mesh_2d )

            !            PARTICLE_ARRAY_SIZE,                      &       ! MCP: problems if this value too small ?
            !            GUARD_SIZE,                               &       ! MCP: problems if this value too small ?

        sim%part_group%use_exact_f0      = UseExactF0

        print *, "WARNING 65373654 -- writing landau parameters in the particle group -- this is temporary..."
        sim%part_group%thermal_speed = sim%thermal_speed_ions   !  temporary or not?
        sim%part_group%alpha_landau = ALPHA    !  temporary or not?
        sim%part_group%k_landau = KX_LANDAU    !  temporary or not?

        call sll_lt_pic_4d_init_landau (                &
            sim%thermal_speed_ions,                     &
            ALPHA, KX_LANDAU,                           &
            sim%part_group )

        sim%ions_number = sim%part_group%number_particles
        SLL_ASSERT( sim%ions_number == NUM_PARTS_X * NUM_PARTS_Y * NUM_PARTS_VX * NUM_PARTS_VY)

        print *, "sim%ions_number (pushed markers) = ", sim%ions_number

        sim%virtual_particle_number =                                               &
              sim%n_virtual_x_for_deposition * sim%mesh_2d%num_cells1               &
            * sim%n_virtual_y_for_deposition * sim%mesh_2d%num_cells2               &
            * sim%n_virtual_vx_for_deposition * sim%part_group%number_parts_vx      &
            * sim%n_virtual_vy_for_deposition * sim%part_group%number_parts_vy

        print *, "sim%virtual_particle_number (deposited particles) = ", sim%virtual_particle_number

    else
        !! initialize the group of PIC particles

        !! -- PIC_VERSION begin
        !        sim%part_group => new_particle_4d_group(    &
        !            NUM_PARTICLES,                          &
        !            PARTICLE_ARRAY_SIZE,                    &
        !            GUARD_SIZE,                             &
        !            QoverM,                                 &
        !            sim%mesh_2d )
        !
        !        call random_seed (SIZE=rand_seed_size)
        !
        !        SLL_ALLOCATE( rand_seed(1:rand_seed_size), ierr )
        !        do j=1, rand_seed_size
        !           rand_seed(j) = (-1)**j*(100 + 15*j)*(2*sim%my_rank + 1)
        !        enddo
        !
        !        pa_gr => sim%part_group
        !        call sll_initial_particles_4d( sim%thermal_speed_ions,          &
        !                                       ALPHA, KX_LANDAU, sim%mesh_2d,   &
        !                                       sim%ions_number,                 &
        !                                       pa_gr,                           &
        !                                       rand_seed, sim%my_rank,          &
        !                                       sim%world_size  )
        !        SLL_DEALLOCATE_ARRAY( rand_seed, ierr )
        !
        !        !$omp parallel
        !        !$omp do
        !        do j=1,sim%ions_number
        !           sim%part_group%p_list(j) = pa_gr%p_list(j)
        !        enddo
        !        !$omp end do
        !
        !        !$omp single
        !        call sll_sort_particles_2d( sim%sorter, sim%part_group )
        !        sim%n_threads =  1
        !        !$omp end single
        !        !$omp end parallel

        !! -- PIC_VERSION end

    end if


    !$omp parallel
#ifdef _OPENMP
    if (OMP_GET_THREAD_NUM() == 0) then
       sim%n_threads =  OMP_GET_NUM_THREADS()
    endif
#else
    sim%n_threads =  1
#endif
    !$omp end parallel

    print*, 'number of threads is ', sim%n_threads

    if (sim%use_cubic_splines) then
      print*, "error (0976765) cubic splines not implemented yet"
      stop

    else
       
       SLL_ALLOCATE(sim%q_accumulator_ptr(1:sim%n_threads), ierr)
       thread_id = 0
       !$omp parallel PRIVATE(thread_id)
#ifdef _OPENMP
       thread_id = OMP_GET_THREAD_NUM()
#endif 
       sim%q_accumulator_ptr(thread_id+1)%q => new_charge_accumulator_2d( sim%mesh_2d )
       !$omp end parallel
       sim%E_accumulator => new_field_accumulator_2d( sim%mesh_2d )

       !! -- --  First charge deposition [begin]  -- --

       if( sim%use_lt_pic_scheme )then

          SLL_ASSERT(thread_id == 0)

          ! [[selalib:src/pic_utilities/lt_pic_4d_utilities.F90::sll_lt_pic_4d_deposit_charge_on_2d_mesh]]
           call sll_lt_pic_4d_deposit_charge_on_2d_mesh( sim%part_group,                    &
                                                         sim%q_accumulator_ptr(1)%q,        &
                                                         sim%n_virtual_x_for_deposition,    &
                                                         sim%n_virtual_y_for_deposition,    &
                                                         sim%n_virtual_vx_for_deposition,   &
                                                         sim%n_virtual_vy_for_deposition,   &
                                                         sim%total_density )
       else
           !! -- PIC_VERSION
           !                   call sll_first_charge_accumulation_2d( sim%part_group, sim%q_accumulator_ptr)!(1)%q )

       end if

       !! -- --  First charge deposition [end]  -- --

    endif
    
  end subroutine init_4d_lt_pic_cartesian






  !> run_4d_lt_pic_cartesian: run the Vlasov-Poisson simulation
  !!
  !! note 1: this is a skeleton-in-progress: some routines are not implemented, some variables are not needed
  !!
  !! note 2: use of cubic spline particles (routines with _CS) is disabled for now
  !!
  !! todo 1: run the code and later write a (non-virtual) remapping step with adequate frequency
  !!
  !! todo 2: use a common type for PIC and LT_PIC (commented calls to PIC structures are labelled with "PIC_VERSION"
  !!
  !! this version written by MCP, ALH

  subroutine run_4d_lt_pic_cartesian( sim )
    class(sll_simulation_4d_vp_lt_pic_cartesian), intent(inout)  :: sim
    sll_int32  :: ierr, it, jj, counter
    sll_int32  :: i, j, k
    sll_real64 :: tmp1, tmp2, tmp3, tmp4
    sll_real32 :: tmp5, tmp6, temp
    sll_real32 :: ttmp(1:4,1:2), ttmp1(1:4,1:2), ttmp2(1:4,1:2)
    sll_real64 :: valeur, val2
    sll_real64, dimension(:,:), pointer :: phi
    sll_int32  :: ncx, ncy, ic_x,ic_y
    sll_int32  :: ic_x1,ic_y1
    sll_real32 :: off_x, off_y,off_x1,off_y1
    sll_real64 :: xmin, ymin
    sll_real64 :: rdx, rdy
    sll_int32  :: icell
    sll_int32  :: gi ! counter index for guard list
    sll_real64 :: Ex, Ey, Ex1, Ey1, Ex_CS, Ey_CS
    sll_real64 :: dt_q_over_m ! dt * qoverm
    sll_real64 :: x, x1  ! for global position
    sll_real64 :: y, y1  ! for global position
    sll_real64 :: dt
    sll_real64 :: pp_vx, pp_vy
    type(sll_lt_pic_4d_particle), dimension(:), pointer :: particles        !! todo: use a common class later
    !    type(sll_particle_4d), dimension(:), pointer :: p
    type(field_accumulator_cell), dimension(:), pointer :: accumE
    type(field_accumulator_CS), dimension(:), pointer :: accumE_CS
    type(sll_lt_pic_4d_particle_guard), dimension(:), pointer :: p_guard     !! todo: use a common type when you feel ready :-)
    !  type(sll_particle_4d_guard), dimension(:), pointer :: p_guard
    sll_real64, dimension(:,:), allocatable  ::  diag_energy
    sll_real64, dimension(:),   allocatable  ::  diag_TOTmoment
    sll_real64, dimension(:),   allocatable  ::  diag_TOTenergy
    sll_real64, dimension(:,:), allocatable :: diag_AccMem! a memory buffer
!    type(sll_time_mark)  :: t2, t3
    sll_real64 :: t2, t3
    sll_real64, dimension(:), allocatable :: rho1d_send
    sll_real64, dimension(:), allocatable :: rho1d_receive
    sll_real64   :: t_init, t_fin, time
    sll_int32 :: save_nb
    sll_int32 :: thread_id
    sll_int32 :: n_threads
    type(sll_charge_accumulator_2d),    pointer :: q_accum
    type(sll_charge_accumulator_2d_CS), pointer :: q_accum_CS
    sll_int32 :: sort_nb
    sll_real64 :: some_val, une_cst
    sll_real64 :: val_lee, exval_ee
    sll_real64 :: tot_ee, val_ee
    sll_real64 :: omega_i, omega_r, psi
    sll_real64 :: bors

    ! ALH - Timing statistics
    sll_real64 :: deposit_time,loop_time,intermediate_deposit_time
    type(sll_time_mark) :: deposit_time_mark,loop_time_mark,intermediate_deposit_time_mark
    sll_int32  :: speed_sample_period = 1 ! <<speed_sample_period>>
    
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
    particles => sim%part_group%p_list      ! previous name was 'p'
    !!!    p_guard => sim%part_group%p_guard
!    dt = sim%dt
    dt = sim%dt
    dt_q_over_m = dt * sim%part_group%qoverm
    print*,  "dt_q_over_m = ", dt_q_over_m
    xmin = sim%mesh_2d%eta1_min
    ymin = sim%mesh_2d%eta2_min
    rdx = 1._f64/sim%mesh_2d%delta_eta1
    rdy = 1._f64/sim%mesh_2d%delta_eta2

    !  ----------------------------------------------------------------------------------------------------
    !  ------
    !  ------  PREPARING THE TIME LOOP:
    !  ------
    !  ------  - begins with:
    !  ------      * (x,y)^0_k, (vx, vy)^0_k  stored in particles => sim%part_group%p_list
    !
    !  ------  - ends with:
    !  ------      * E^0 stored in sim%E1, sim%E2
    !  ------      * (x,y)^0_k, (vx, vy)^{-1/2}_k  stored in particles => sim%part_group%p_list
    !  ------
    !  ----------------------------------------------------------------------------------------------------

    !! -- --  ?? [begin]  -- --

    if (sim%use_cubic_splines) then
      print*, "error (0976765) cubic splines not implemented yet"
      stop
    else
       accumE => sim%E_accumulator%e_acc
       call sum_accumulators( sim%q_accumulator_ptr, n_threads, ncx*ncy )
       call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator_ptr(1)%q, sim%rho )     ! this name not clear enough
    endif

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
    
    do j = 1, ncy+1
       do i = 1, ncx+1
          sim%rho(i, j) = rho1d_receive(i+(j-1)*(ncx+1))
       enddo
    enddo

    !! -- --  MPI communications of rho [end]  -- --

    if (sim%my_rank == 0) then
    it = 0
    call sll_gnuplot_2d(xmin, sim%mesh_2d%eta1_max, ncx+1, ymin,            &
                        sim%mesh_2d%eta2_max, ncy+1,                        &
                        sim%rho, 'rho_init_standPUSH', it, ierr )
    endif

    !! -- --  Poisson solver (computing E^0) -- --

    call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, -sim%rho )         !  I don't like this - sign...

    if (sim%my_rank == 0) then
    it = 0
    print *, "writing Ex, Ey in gnuplot format for iteration # it = ", it
    call sll_gnuplot_2d(xmin, sim%mesh_2d%eta1_max, ncx+1, ymin,            &
                        sim%mesh_2d%eta2_max, ncy+1,                        &
                        sim%E1, 'Ex', it, ierr )

    call sll_gnuplot_2d(xmin, sim%mesh_2d%eta1_max, ncx+1, ymin,            &
                        sim%mesh_2d%eta2_max, ncy+1,                        &
                        sim%E2, 'Ey', it, ierr )
    endif


    !! -- --  diagnostics: compute energy [begin]  -- --

    call electric_energy( tot_ee, sim%E1, sim%E2, ncx, &
         ncy, sim%mesh_2d%eta1_max, sim%mesh_2d%eta2_max )
    bors = 0.0_f64
    !$omp parallel
    !$omp do reduction(+:bors)
    do k = 1, sim%ions_number
       bors = bors + particles(k)%vx**2 + particles(k)%vy**2
    enddo
    !$omp end do
    !$omp end parallel
    diag_TOTenergy(0) = bors * 0.5_f64*(sim%mesh_2d%eta1_max - sim%mesh_2d%eta1_min)  &
         * (sim%mesh_2d%eta2_max - sim%mesh_2d%eta2_min)/( sim%world_size*sim%ions_number)  &
         + tot_ee * 0.5_f64

    !! -- --  diagnostics [end]  -- --


    !! -- --  half v-push (computing v^0 -> v^{-1/2})  [begin]  -- --

    if (sim%use_cubic_splines) then
      print*, "error (0976765) cubic splines not implemented yet"
      stop

    else

       call reset_field_accumulator_to_zero( sim%E_accumulator )
       call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
       !$omp parallel do PRIVATE (pp_vx, pp_vy, Ex, Ey, tmp5, tmp6)
       !$&omp FIRSTPRIVATE(dt_q_over_m)
       do k = 1, sim%ions_number
          pp_vx = particles(k)%vx
          pp_vy = particles(k)%vy
          call get_poisson_cell_index(sim%mesh_2d, particles(k)%ic_x, particles(k)%ic_y, icell)
          SLL_INTERPOLATE_FIELD_EXTENDED(Ex,Ey,accumE,particles(k),tmp5,tmp6,icell)
          particles(k)%vx = pp_vx - 0.5_f64 * dt_q_over_m * Ex
          particles(k)%vy = pp_vy - 0.5_f64 * dt_q_over_m * Ey
       enddo
       !$omp end parallel do
    endif

    !! -- --  half v-push  [end]  -- --

    ! une_cst = (pi / KX_LANDAU) * (4._f64 * ALPHA * er)**2
    une_cst = (sll_pi / sim%elec_params(1)) * (4._f64 * sim%elec_params(2) * sim%elec_params(3))**2
    omega_i = sim%elec_params(6)
    omega_r = sim%elec_params(5)
    psi = sim%elec_params(4)

    if (sim%my_rank ==0) open(65,file='logE_standPush.dat')
#ifdef _OPENMP
    t2 = omp_get_wtime()!!   call sll_set_time_mark(t2)
#endif

    !  ----------------------------------------------------------------------------------------------------
    !  ------
    !  ------  TIME LOOP
    !  ------
    !  ------  - begins with:
    !  ------      * E^n stored in sim%E1, sim%E2
    !  ------      * (x,y)^n_k, (vx, vy)^{n-1/2}_k  stored in particles => sim%part_group%p_list
    !  ------    (where n = it)
    !  ------
    !  ------  - ends with:
    !  ------      * E^{n+1} stored in sim%E1, sim%E2
    !  ------      * (x,y)^{n+1}_k, (vx, vy)^{n+1/2}_k  stored in particles => sim%part_group%p_list
    !  ------
    !  ----------------------------------------------------------------------------------------------------

    ! Time statistics
    call sll_set_time_mark(loop_time_mark)
    deposit_time=0
    
    do it = 0, sim%num_iterations-1

       print *, "BEGIN one loop in time, it = ", it
       !! -- --  diagnostics (computing energy) [begin]  -- --

       if (sim%my_rank == 0) then
          exval_ee = une_cst * exp(2._f64 * omega_i * real(it,f64) * sim%dt)                          &
                             * ( 0.5_f64 + 0.5_f64 * cos(2._f64 * (omega_r * real(it,f64) * sim%dt - psi)) )
          exval_ee = max(exval_ee, 1e-30_f64)    ! to avoid taking the log of 0
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

       !! -- --  diagnostics [end]  -- --

       !! -- --  sort particles  -- --
       if( .not. sim%use_lt_pic_scheme )then
           if (mod(it+1,sort_nb)==0) then
            !  PIC_VERSION            call sll_sort_particles_2d( sim%sorter, sim%part_group )
           endif
       end if

       !  ------------------------------------------------------------------
       !  ------
       !  ------  PUSH PARTICLES [begin]
       !  ------
       !  ------------------------------------------------------------------

       some_val = 0.0_f64

      if (sim%use_cubic_splines) then

        print*, "error (0976765) cubic splines not implemented yet"
        stop

      else
       
          !$omp parallel PRIVATE(x,y,x1,y1,Ex,Ey,gi,tmp1,tmp2,tmp5,tmp6,off_x,off_y,ic_x,ic_y,thread_id,p_guard,q_accum)
          !$&omp FIRSTPRIVATE(dt,dt_q_over_m,ncx,xmin,ymin,rdx,rdy)
#ifdef _OPENMP
          thread_id = OMP_GET_THREAD_NUM()
#endif
          call reset_charge_accumulator_to_zero ( sim%q_accumulator_ptr(thread_id+1)%q )
          q_accum => sim%q_accumulator_ptr(thread_id+1)%q
          p_guard => sim%part_group%p_guard(thread_id+1)%g_list
          gi = 0
          !$omp do!! reduction(+:some_val)

          !! -- --  particle loop: treat the particles by pair for faster treatment ?  -- --

          do k = 1, sim%ions_number

             !! -- --  v-push (v^{n-1/2} -> v^{n+1/2}  using  E^n)  [begin]  -- --

             call get_poisson_cell_index(sim%mesh_2d, particles(k)%ic_x, particles(k)%ic_y, icell)
             SLL_INTERPOLATE_FIELD_EXTENDED(Ex,Ey,accumE,particles(k),tmp5,tmp6, icell)
            ! SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,particles(k),tmp5,tmp6)
             particles(k)%vx = particles(k)%vx + dt_q_over_m * Ex
             particles(k)%vy = particles(k)%vy + dt_q_over_m * Ey

             !! -- --  v-push [end]  -- --

             !! -- --  x-push (x^n -> x^{n+1} using v^{n+1/2})  [begin]  -- --

             GET_PARTICLE_POSITION_EXTENDED(particles(k),sim%mesh_2d,x,y)
             x = x + dt * particles(k)%vx
             y = y + dt * particles(k)%vy

             !! -- --  x-push [end]  -- --


             !! note: with the ltp_bsl charge deposition we will not be able to deposit the charge of the "just-pushed" particles
             !!       since we need to first compute the weights of the virtual particles (with the remapping algorithm)


             !! -- --  LTPIC: put outside particles back in domain                              [begin]  -- --
             !! -- --  PIC: deposit charge (if particle is inside, otherwise reserve it)        [begin]  -- --

             if( sim%part_group%track_markers_outside_domain                            &
                .or. ( in_bounds_periodic( x, y, sim%mesh_2d,                           &
                                           sim%part_group%domain_is_x_periodic,         &
                                           sim%part_group%domain_is_y_periodic )        &
                ) ) then ! finish push
                SET_PARTICLE_POSITION_EXTENDED(particles(k),xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
                if( .not. sim%use_lt_pic_scheme )then
                    print*,  "WARNING: charge deposition discarded for pic case, please update"
                    ! SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,particles(k),tmp5,tmp6)
                end if
             else
                ! particle outside domain
                if( sim%use_lt_pic_scheme )then
                     call apply_periodic_bc( sim%mesh_2d, x, y)
                     SET_PARTICLE_POSITION_EXTENDED(particles(k),xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
                else
                  ! store reference for later processing
                  gi = gi + 1
                  p_guard(gi)%p => particles(k)
                end if
             end if

             !! -- --  LTPIC: put outside particles back in domain                              [end]  -- --
             !! -- --  PIC: deposit charge (if particle is inside, otherwise reserve it)        [end]  -- --

          enddo


          !! -- --  deposit charge LTPIC [begin]  -- --

            print *, "deposit charge begin"

          if( sim%use_lt_pic_scheme )then
              SLL_ASSERT(thread_id == 0)
             
              ! [[file:~/mcp/selalib/src/pic_utilities/lt_pic_4d_utilities.F90::sll_lt_pic_4d_deposit_charge_on_2d_mesh]]
              call sll_set_time_mark(deposit_time_mark)
              call sll_set_time_mark(intermediate_deposit_time_mark)
              call sll_lt_pic_4d_deposit_charge_on_2d_mesh( sim%part_group,                     &
                                                            sim%q_accumulator_ptr(1)%q,         &
                                                            sim%n_virtual_x_for_deposition,     &
                                                            sim%n_virtual_y_for_deposition,     &
                                                            sim%n_virtual_vx_for_deposition,    &
                                                            sim%n_virtual_vy_for_deposition,    &
                                                            sim%total_density )
              deposit_time=deposit_time+sll_time_elapsed_since(deposit_time_mark)
              intermediate_deposit_time=sll_time_elapsed_since(intermediate_deposit_time_mark)
          else
              ! nothing to do, charge already deposited in the push loop
          end if

            print *, "deposit charge end"

          !! -- --  deposit charge LTPIC [end]  -- --

          !! -- --  [PIC ONLY] process the reserved particles (in the guard list) [begin]  -- --

          if( .not. sim%use_lt_pic_scheme )then
              !$omp end do
              sim%part_group%num_postprocess_particles(thread_id+1) = gi
              !$omp end parallel

              ! Process the particles in the guard list. In the periodic case, no
              ! destruction of particles is needed, so this is simple.

              !$omp parallel PRIVATE(x,y,ic_x,ic_y,off_x,off_y,tmp1,tmp2,tmp5,tmp6,p_guard,q_accum,p,thread_id)
              !$&omp FIRSTPRIVATE(dt,ncx,xmin,ymin,rdx,rdy)
#ifdef _OPENMP
              thread_id = OMP_GET_THREAD_NUM()
#endif
              q_accum => sim%q_accumulator_ptr(thread_id+1)%q
              p_guard => sim%part_group%p_guard(thread_id+1)%g_list
                ! !  !          p => sim%part_group%p_list
              do k = 1, sim%part_group%num_postprocess_particles(thread_id+1)

                print*,  "WARNING (895848764): macros discarded for pic case, please update"

                !                 GET_PARTICLE_POSITION(p_guard(k)%p,sim%mesh_2d,x,y)
                !                 x = x + dt * p_guard(k)%p%vx
                !                 y = y + dt * p_guard(k)%p%vy
                !                 call apply_periodic_bc( sim%mesh_2d, x, y)
                !
                !                 SET_PARTICLE_POSITION(p_guard(k)%p,xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
                !                 SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p_guard(k)%p,tmp5,tmp6)
              end do
              !$omp end parallel
              !       ! reset any counters
              gi = 0
          else
              ! nothing to do, particles have been put back in the periodic domain already
          end if
          !! -- --  [PIC ONLY] process the reserved particles [end]  -- --

       endif

       !  ------------------------------------------------------------------
       !  ------
       !  ------  PUSH PARTICLES [end]
       !  ------
       !  ------------------------------------------------------------------

       !! -- --  parallel communications for rho ?? [begin]  -- --

       if (sim%use_cubic_splines) then
         print*, "error (0976765) cubic splines not implemented yet"
         stop

       else
          call sum_accumulators( sim%q_accumulator_ptr, n_threads, ncx*ncy )
          call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator_ptr(1)%q, sim%rho )     ! this name not clear enough
       endif
       do j = 1, ncy+1
          do i = 1, ncx+1
             rho1d_send(i+(j-1)*(ncx+1)) = sim%rho(i, j)
             rho1d_receive(i+(j-1)*(ncx+1)) = 0._f64
          enddo
       enddo

       call sll_collective_allreduce( sll_world_collective, rho1d_send, (ncx+1)*(ncy+1), &
            MPI_SUM, rho1d_receive   )

       do j = 1, ncy+1
          do i = 1, ncx+1
             sim%rho(i, j) = rho1d_receive(i+(j-1)*(ncx+1))
          enddo
       enddo

       !! -- --  parallel communications for rho ?? [end]  -- --

       !! -- --  Poisson solver (computing E^{n+1}) -- --

       call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, -sim%rho )      !  I don't like this - sign...

       !! -- --  diagnostics (plotting) [begin]  -- --

        if (sim%my_rank == 0 .and. mod(it, sim%plot_period)==0 ) then

            print *, "writing f slice in gnuplot format for iteration # it = ", it, " / ", sim%num_iterations
!            call plot_f_slice(sim%part_group, "f_slice", it)

            call plot_f_slice_x_vx(sim%part_group,           &
                                   sim%part_group%remapping_grid%eta1_min,   &
                                   sim%part_group%remapping_grid%eta1_max,   &
                                   sim%part_group%remapping_grid%eta2_min,   &
                                   sim%part_group%remapping_grid%eta2_max,   &
                                   sim%part_group%remapping_grid%eta3_min,   &
                                   sim%part_group%remapping_grid%eta3_max,   &
                                   sim%part_group%remapping_grid%eta4_min,   &
                                   sim%part_group%remapping_grid%eta4_max,   &
                                   sim%part_group%remapping_grid%num_cells1, &
                                   sim%part_group%remapping_grid%num_cells2, &
                                   sim%part_group%remapping_grid%num_cells3, &
                                   sim%part_group%remapping_grid%num_cells4, &
                                   sim%n_virtual_x_for_deposition,    &
                                   sim%n_virtual_y_for_deposition,    &
                                   sim%n_virtual_vx_for_deposition,   &
                                   sim%n_virtual_vy_for_deposition,   &
                                   "f_slice", it)

        end if
        if (sim%my_rank == 0 .and. mod(it+1, sim%remap_period)==0 ) then

            print *, "remapping f..."
            call sll_lt_pic_4d_remap(sim%part_group)

            print *, "writing (remapped) f slice in gnuplot format for iteration # it = ", it, " / ", sim%num_iterations

            call plot_f_slice_x_vx(sim%part_group,           &
                                   sim%part_group%remapping_grid%eta1_min,   &
                                   sim%part_group%remapping_grid%eta1_max,   &
                                   sim%part_group%remapping_grid%eta2_min,   &
                                   sim%part_group%remapping_grid%eta2_max,   &
                                   sim%part_group%remapping_grid%eta3_min,   &
                                   sim%part_group%remapping_grid%eta3_max,   &
                                   sim%part_group%remapping_grid%eta4_min,   &
                                   sim%part_group%remapping_grid%eta4_max,   &
                                   sim%part_group%remapping_grid%num_cells1, &
                                   sim%part_group%remapping_grid%num_cells2, &
                                   sim%part_group%remapping_grid%num_cells3, &
                                   sim%part_group%remapping_grid%num_cells4, &
                                   sim%n_virtual_x_for_deposition,    &
                                   sim%n_virtual_y_for_deposition,    &
                                   sim%n_virtual_vx_for_deposition,   &
                                   sim%n_virtual_vy_for_deposition,   &
                                   "f_slice_remapped", it)
        end if

        if (sim%my_rank == 0 .and. mod(it+1, sim%plot_period)==0 ) then

            print *, "writing Ex, Ey  in gnuplot format for iteration # it = ", it+1, " / ", sim%num_iterations
            call sll_gnuplot_2d(xmin, sim%mesh_2d%eta1_max, ncx+1, ymin,            &
                                sim%mesh_2d%eta2_max, ncy+1,                        &
                                sim%E1, 'Ex', it+1, ierr )

            call sll_gnuplot_2d(xmin, sim%mesh_2d%eta1_max, ncx+1, ymin,            &
                                sim%mesh_2d%eta2_max, ncy+1,                        &
                                sim%E2, 'Ey', it+1, ierr )
            print *, "done."

        endif



        ! ALH - display intermediate times (cf [[speed_sample_period]]) to determine how much the speed of the method
        ! varies with the age of the last remap

        if (sim%my_rank == 0 .and. mod(it,speed_sample_period)==0 ) then
           write(*,'(A,I5.1,A,ES8.2,A)') ' intermediate speed at iteration ',it,' : ', &
                1/intermediate_deposit_time * sim%virtual_particle_number,' deposits/sec'
        endif
            
       !! -- --  diagnostics (plotting E ) [end]  -- --

       !! -- --  diagnostics (computing energy) [begin]  -- --
        print *, "diag energy"

       if (sim%use_cubic_splines) then
         print*, "error (0976765) cubic splines not implemented yet"
         stop

       else
          call reset_field_accumulator_to_zero( sim%E_accumulator )
          call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
          some_val = 0.0_f64
          !$omp parallel PRIVATE(Ex,Ey,Ex1,Ey1,tmp5,tmp6)
          !$omp do reduction(+:some_val)
          do k = 1, sim%ions_number
             call get_poisson_cell_index(sim%mesh_2d, particles(k)%ic_x, particles(k)%ic_y, icell)
             SLL_INTERPOLATE_FIELD_EXTENDED(Ex,Ey,accumE,particles(k),tmp5,tmp6, icell)
             ! SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,particles(k),tmp5,tmp6)
             some_val = some_val &
                    + (particles(k)%vx + 0.5_f64 * dt_q_over_m * Ex)**2 &
                    + (particles(k)%vy + 0.5_f64 * dt_q_over_m * Ey)**2
          enddo
          !$omp end do
          !$omp end parallel
          call electric_energy( tot_ee, sim%E1, sim%E2, ncx, &
               ncy, sim%mesh_2d%eta1_max, sim%mesh_2d%eta2_max )
          diag_TOTenergy(it) = some_val* 0.5_f64 *(sim%mesh_2d%eta1_max - sim%mesh_2d%eta1_min) * &
               (sim%mesh_2d%eta2_max - sim%mesh_2d%eta2_min)/( sim%world_size*sim%ions_number)  &
               + tot_ee * 0.5_f64
       endif

       !! -- --  diagnostics [end]  -- --

    print *, "end one loop in time"
    enddo

    !  ----------------------------------------------------------------------------------------------------
    !  ------
    !  ------  TIME LOOP [end]
    !  ------
    !  ----------------------------------------------------------------------------------------------------

#ifdef _OPENMP
    time = omp_get_wtime()!! time = sll_time_elapsed_since(t2)

    if (sim%my_rank ==0) then 
       close(65)
       open(93,file='time_parts_sec_omp.dat',position='append')
       write(93,*) '# Nb of threads  ||  time (sec)  ||  average pushes/sec'
       write(93,*) sim%n_threads, time-t2, int(sim%num_iterations,i64)*int(sim%ions_number,i64)/(time-t2)
       close(93)
       open(65,file='Energie_totale_standPush.dat')
       do it=0,sim%num_iterations-1
          write(65,*) it*dt,  diag_TOTenergy(it), (diag_TOTenergy(it)-diag_TOTenergy(0))/diag_TOTenergy(0)
       enddo
       close(65)
    endif
#endif

    ! ALH - Time statistics. The number of deposits is computed separately to take the number of LTP BSL virtual
    ! particles into account.
    
    loop_time=sll_time_elapsed_since(loop_time_mark)
    if(sim%use_lt_pic_scheme)then

       ! Different label for loop_time in PIC and LTPIC to make sure that the two definitions are not confused
       
       write(*,'(A,ES8.2,A,ES8.2,A)') 'LTPIC stats: ',                           &
            1 / loop_time * sim%num_iterations * sim%ions_number,                &
            ' markers/sec ',                                                     &
            1 / deposit_time * sim%num_iterations * sim%virtual_particle_number, &
            ' deposits/sec'
    else
       write(*,'(A,ES8.2,A)') 'PIC stats: ',1 / loop_time * sim%num_iterations * sim%ions_number,' pushes/sec '
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

  end subroutine run_4d_lt_pic_cartesian


  subroutine delete_4d_lt_pic_cartesian( sim )
    type(sll_simulation_4d_vp_lt_pic_cartesian) :: sim
  end subroutine delete_4d_lt_pic_cartesian


  function in_bounds( x, y, mesh ) result(res)
    logical :: res
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    type(sll_cartesian_mesh_2d), pointer :: mesh

    res = (x >= mesh%eta1_min) .and. (x <= mesh%eta1_max) .and. &
          (y >= mesh%eta2_min) .and. (y <= mesh%eta2_max)
  end function in_bounds

!  subroutine apply_periodic_bc( mesh, x, y )
!    type(sll_cartesian_mesh_2d), pointer :: mesh
!    sll_real64, intent(inout) :: x
!    sll_real64, intent(inout) :: y
!
!    x = modulo(x,mesh%eta1_max - mesh%eta1_min)! Ca marche quand xmin=ymin=0
!    y = modulo(y,mesh%eta2_max - mesh%eta2_min)! Sinon, il faut modulo(...) + xmin/mesh%delta_x
!  end subroutine apply_periodic_bc

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

    
end module sll_simulation_4d_vp_lt_pic_cartesian_module
