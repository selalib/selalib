module sll_pic_simulation_4d_cartesian_module

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
  use sll_particle_initializers_4d
  use sll_particle_sort_module
  use sll_charge_to_density_module
  use sll_pic_utilities
  use sll_module_poisson_2d_fft
  use sll_module_poisson_2d_base
  use sll_representation_conversion_module
  use sll_gnuplot
  use sll_timer
  use sll_collective 
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none

  type, extends(sll_simulation_base_class) :: sll_pic_simulation_4d_cartesian
     ! Physics/numerical parameters
     sll_real64 :: dt
     sll_int32  :: num_iterations
     sll_real64 :: thermal_speed_ions
     sll_int32  :: ions_number
     sll_int32  :: guard_size
     sll_int32  :: array_size
     sll_real64, dimension(1:6) :: elec_params
     type(sll_particle_group_4d),  pointer :: part_group
     type(sll_cartesian_mesh_2d),    pointer :: m2d
     type(sll_particle_sorter_2d), pointer :: sorter
     type(sll_charge_accumulator_2d_ptr), dimension(:), pointer  :: q_accumulator
     type(electric_field_accumulator), pointer :: E_accumulator
     logical :: use_cubic_splines
     type(sll_charge_accumulator_2d_CS_ptr), dimension(:), pointer  :: q_accumulator_CS
     type(electric_field_accumulator_CS), pointer :: E_accumulator_CS
     sll_real64, dimension(:,:), pointer :: rho
     type(poisson_2d_fft_solver), pointer :: poisson
     sll_real64, dimension(:,:), pointer :: E1, E2
     sll_int32 :: my_rank
     sll_int32 :: world_size
     sll_int32 :: n_threads
   contains
     procedure, pass(sim) :: init_from_file => init_4d_pic_cartesian
     procedure, pass(sim) :: run => run_4d_pic_cartesian
  end type sll_pic_simulation_4d_cartesian

  interface sll_delete
     module procedure delete_4d_pic_cart
  end interface sll_delete

!!$  interface initialize
!!$     module procedure initialize_4d_qns_general
!!$  end interface initialize

contains

  subroutine init_4d_pic_cartesian( sim, filename )
    intrinsic :: trim
    class(sll_pic_simulation_4d_cartesian), intent(inout) :: sim
    character(len=*), intent(in)                          :: filename
    sll_int32   :: IO_stat
    sll_int32   :: ierr
    sll_int32   :: j
    sll_real64  :: dt
    sll_int32   :: number_iterations
    sll_int32   :: NUM_PARTICLES, GUARD_SIZE, PARTICLE_ARRAY_SIZE
    sll_real64  :: THERM_SPEED
    sll_real64  :: QoverM, ALPHA
    logical     :: UseCubicSplines
    sll_int32   :: NC_X,  NC_Y
    sll_real64  :: XMIN, KX, XMAX, YMIN, YMAX
    sll_real64  :: er, psi, omega_i, omega_r
    sll_int32, parameter  :: input_file = 99
    sll_int32, dimension(:), allocatable  :: rand_seed
    sll_int32   :: rand_seed_size
    sll_int32  :: thread_id
    type(sll_particle_group_4d), pointer  :: pa_gr

    namelist /sim_params/ NUM_PARTICLES, GUARD_SIZE, &
                          PARTICLE_ARRAY_SIZE, &
                          THERM_SPEED, dt, number_iterations, &
                          QoverM, ALPHA, UseCubicSplines
    namelist /grid_dims/  NC_X, NC_Y, XMIN, KX, YMIN, YMAX
    namelist /elec_params/ er, psi, omega_r, omega_i
    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_vp4d_par_cart() failed to open file ', filename
       STOP
    end if
    read(input_file, sim_params)
    read(input_file, grid_dims)
    read(input_file, elec_params)
    close(input_file)

    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)
    
    print*, 'sim%world_size=',sim%world_size, ' sim%my_rank=', sim%my_rank

    XMAX = (2._f64*sll_pi/KX)
    sim%use_cubic_splines = UseCubicSplines
    sim%thermal_speed_ions = THERM_SPEED
    sim%ions_number = NUM_PARTICLES
    sim%guard_size = GUARD_SIZE  
    sim%array_size = PARTICLE_ARRAY_SIZE
    sim%dt = dt
    sim%num_iterations = number_iterations

    sim%elec_params = (/KX, ALPHA, er, psi, omega_r, omega_i /)
    
    sim%m2d =>  new_cartesian_mesh_2d( NC_X, NC_Y, &
                XMIN, XMAX, YMIN, YMAX )

    sim%part_group => new_particle_4d_group( &
         NUM_PARTICLES, &
         PARTICLE_ARRAY_SIZE, &
         GUARD_SIZE, &
         QoverM,     &
         sim%m2d )

    sim%sorter => sll_new_particle_sorter_2d( sim%m2d )

    sim%poisson => new_poisson_2d_fft_solver( sim%m2d%eta1_min,    &
                                              sim%m2d%eta1_max,    & 
                                              sim%m2d%num_cells1,  &
                                              sim%m2d%eta2_min,    &
                                              sim%m2d%eta2_max,    &
                                              sim%m2d%num_cells2   )

    call random_seed (SIZE=rand_seed_size)

    SLL_ALLOCATE( rand_seed(1:rand_seed_size), ierr )
    do j=1, rand_seed_size
       rand_seed(j) = (-1)**j*(100 + 15*j)*(2*sim%my_rank + 1)
    enddo

    pa_gr => sim%part_group
    call sll_initial_particles_4d( sim%thermal_speed_ions, & 
                                   ALPHA, KX, sim%m2d,     &
                                   sim%ions_number,        &
                                   pa_gr, &
                                   rand_seed, sim%my_rank, &
                                   sim%world_size  )
    SLL_DEALLOCATE_ARRAY( rand_seed, ierr )
    !$omp parallel
    !$omp do
    do j=1,sim%ions_number
       sim%part_group%p_list(j) = pa_gr%p_list(j)
    enddo
    !$omp end do

    !$omp single
    call sll_sort_particles_2d( sim%sorter, sim%part_group )
    sim%n_threads =  1
    !$omp end single

#ifdef _OPENMP
    if (OMP_GET_THREAD_NUM() == 0) then
       sim%n_threads =  OMP_GET_NUM_THREADS()
    endif
#endif
    !$omp end parallel
    print*, 'number of threads is ', sim%n_threads

    if (sim%use_cubic_splines) then
      print*, "error (0976765) cubic splines not implemented yet"
      stop

    else
       
       SLL_ALLOCATE(sim%q_accumulator(1:sim%n_threads), ierr)
       thread_id = 0
       !$omp parallel PRIVATE(thread_id)
#ifdef _OPENMP
       thread_id = OMP_GET_THREAD_NUM()
#endif 
       sim%q_accumulator(thread_id+1)%q => new_charge_accumulator_2d( sim%m2d )
       !$omp end parallel
       sim%E_accumulator => new_field_accumulator_2d( sim%m2d )
       call sll_first_charge_accumulation_2d( sim%part_group, sim%q_accumulator)!(1)%q )
    endif
    
  end subroutine init_4d_pic_cartesian






  !> run_4d_lt_pic_cartesian: run the Vlasov-Poisson simulation
  !!
  !! note 1: this is a skeleton-in-progress: some routines are not implemented, some variables are not needed
  !!
  !! note 2: use of cubic spline particles (routines with _CS) is disabled for now
  !!
  !! todo 1: write the ltp_bsl charge deposition algorithm in place of the simple "pic" charge deposition algorithm
  !!
  !! todo 2: check that the remaining procedures are correct with ltp particles
  !!
  !! todo 3: run the code and later write a (non-virtual) remapping step with adequate frequency
  !!
  !! this version written by MCP, ALH

  subroutine run_4d_lt_pic_cartesian( sim )
    class(sll_pic_simulation_4d_cartesian), intent(inout)  :: sim
    sll_int32  :: ierr, it, jj, counter
    sll_int32  :: i, j
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
    sll_int32  :: gi ! counter index for guard list
    sll_real64 :: Ex, Ey, Ex1, Ey1, Ex_CS, Ey_CS
    sll_real64 :: dt_q_over_m ! dt * qoverm
    sll_real64 :: x, x1  ! for global position
    sll_real64 :: y, y1  ! for global position
    sll_real64 :: dt, ttime
    sll_real64 :: pp_vx, pp_vy
    type(sll_particle_4d), dimension(:), pointer :: p
    type(field_accumulator_cell), dimension(:), pointer :: accumE
    type(field_accumulator_CS), dimension(:), pointer :: accumE_CS
    type(sll_particle_4d_guard), dimension(:), pointer :: p_guard
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
    
    ncx = sim%m2d%num_cells1
    ncy = sim%m2d%num_cells2
    n_threads = sim%n_threads
    thread_id = 0
    save_nb = sim%num_iterations/2

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
    dt_q_over_m = dt * sim%part_group%qoverm
    !!!    p_guard => sim%part_group%p_guard
    dt = sim%dt
    xmin = sim%m2d%eta1_min
    ymin = sim%m2d%eta2_min
    rdx = 1._f64/sim%m2d%delta_eta1
    rdy = 1._f64/sim%m2d%delta_eta2

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
       call sum_accumulators( sim%q_accumulator, n_threads, ncx*ncy )
       call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator(1)%q, sim%rho )     ! this name not clear enough
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
    call sll_gnuplot_2d(xmin, sim%m2d%eta1_max, ncx+1, ymin,            &
                        sim%m2d%eta2_max, ncy+1,                        &
                        sim%rho, 'rho_init_standPUSH', it, ierr )
    endif

    !! -- --  Poisson solver (computing E^0) -- --

    call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, -sim%rho )         !  I don't like this - sign...


    !! -- --  diagnostics: compute energy [begin]  -- --

    call electric_energy( tot_ee, sim%E1, sim%E2, sim%m2d%num_cells1, &
         sim%m2d%num_cells2, sim%m2d%eta1_max, sim%m2d%eta2_max )
    bors = 0.0_f64
    !$omp parallel
    !$omp do reduction(+:bors)
    do i =1, sim%ions_number,2
       bors = bors + particles(i)%vx**2 + particles(i)%vy**2 &
              + particles(i+1)%vx**2 + particles(i+1)%vy**2
    enddo
    !$omp end do
    !$omp end parallel
    diag_TOTenergy(0) = bors * 0.5_f64*(sim%m2d%eta1_max - sim%m2d%eta1_min)  &
         * (sim%m2d%eta2_max - sim%m2d%eta2_min)/( sim%world_size*sim%ions_number)  &
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
       do i = 1, sim%ions_number
          pp_vx = particles(i)%vx
          pp_vy = particles(i)%vy
          SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,particles(i),tmp5,tmp6)
          particles(i)%vx = pp_vx - 0.5_f64 * dt_q_over_m * Ex
          particles(i)%vy = pp_vy - 0.5_f64 * dt_q_over_m * Ey
       enddo
       !$omp end parallel do   
    endif

    !! -- --  half v-push  [end]  -- --

    une_cst = 16._f64*(sim%elec_params(2))**2 *(sim%elec_params(3))**2 &
         * sll_pi/sim%elec_params(1) !pi/KX*(4._f64*ALPHA*er)**2
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

    do it = 0, sim%num_iterations-1

       !! -- --  diagnostics (computing energy) [begin]  -- --

       if (sim%my_rank == 0) then
          exval_ee = une_cst * exp(2._f64*omega_i*real(it,f64)*sim%dt) * &
             ( 0.5_f64 + 0.5_f64*cos(2._f64*(omega_r*real(it,f64)*sim%dt-psi)) )
          call normL2_field_Ex ( val_lee, val_ee, ncx, ncy, &
                                 sim%E1,  &
                                 sim%m2d%delta_eta1, sim%m2d%delta_eta2 )
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

       if (mod(it+1,sort_nb)==0) then 
          call sll_sort_particles_2d( sim%sorter, sim%part_group )
       endif

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
          call reset_charge_accumulator_to_zero ( sim%q_accumulator(thread_id+1)%q )
          q_accum => sim%q_accumulator(thread_id+1)%q
          p_guard => sim%part_group%p_guard(thread_id+1)%g_list
          gi = 0
          !$omp do!! reduction(+:some_val)

          !! -- --  particle loop: treat the particles by pair for faster treatment ?  -- --
          do i = 1, sim%ions_number,2

             !! -- --  v-push (v^{n-1/2} -> v^{n+1/2}  using  E^n)  [begin]  -- --

             SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,particles(i),tmp5,tmp6)
             particles(i)%vx = particles(i)%vx + dt_q_over_m * Ex
             particles(i)%vy = particles(i)%vy + dt_q_over_m * Ey
             SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,particles(i+1),tmp5,tmp6)
             particles(i+1)%vx = particles(i+1)%vx + dt_q_over_m * Ex
             particles(i+1)%vy = particles(i+1)%vy + dt_q_over_m * Ey

             !! -- --  v-push [end]  -- --

             !! -- --  x-push (x^n -> x^{n+1} using v^{n+1/2})  [begin]  -- --

             GET_PARTICLE_POSITION(particles(i),sim%m2d,x,y)
             x = x + dt * particles(i)%vx
             y = y + dt * particles(i)%vy
             GET_PARTICLE_POSITION(particles(i+1),sim%m2d,x1,y1)
             x1 = x1 + dt * particles(i+1)%vx
             y1 = y1 + dt * particles(i+1)%vy

             !! -- --  x-push [end]  -- --


             !! )) -- here with the ltp_bsl charge deposition we will not be able to deposit the charge of the "just-pushed" particles
             !! ))    since we need to first compute the weights of the virtual particles (with the remapping algorithm)

             !! -- --  deposit charge (if particle is inside, otherwise reserve it) [begin]  -- --

             if (in_bounds( x, y, sim%m2d )) then ! finish push
                SET_PARTICLE_POSITION(particles(i),xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
                SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,particles(i),tmp5,tmp6)
             else ! store reference for later processing
                gi = gi + 1
                p_guard(gi)%p => particles(i)
             end if
             if (in_bounds( x1, y1, sim%m2d )) then ! finish push
                SET_PARTICLE_POSITION(particles(i+1),xmin,ymin,ncx,x1,y1,ic_x1,ic_y1,off_x1,off_y1,rdx,rdy,tmp1,tmp2)
                SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,particles(i+1),tmp5,tmp6)
             else ! store reference for later processing
                gi = gi + 1
                p_guard(gi)%p => particles(i+1)
             end if

             !! -- --  deposit charge [end]  -- --

          enddo

          !! -- --  process the reserved particles (in the guard list) [begin]  -- --

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
          q_accum => sim%q_accumulator(thread_id+1)%q
          p_guard => sim%part_group%p_guard(thread_id+1)%g_list
! !  !          p => sim%part_group%p_list
          do i=1, sim%part_group%num_postprocess_particles(thread_id+1)
             GET_PARTICLE_POSITION(p_guard(i)%p,sim%m2d,x,y)
             x = x + dt * p_guard(i)%p%vx
             y = y + dt * p_guard(i)%p%vy
             call apply_periodic_bc( sim%m2d, x, y)
             SET_PARTICLE_POSITION(p_guard(i)%p,xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
             SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p_guard(i)%p,tmp5,tmp6)
          end do
          !$omp end parallel
          !       ! reset any counters
          gi = 0

          !! -- --  process the reserved particles [end]  -- --

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
          call sum_accumulators( sim%q_accumulator, n_threads, ncx*ncy )
          call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator(1)%q, sim%rho )     ! this name not clear enough
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


       !! -- --  diagnostics (computing energy) [begin]  -- --

       if (sim%use_cubic_splines) then
         print*, "error (0976765) cubic splines not implemented yet"
         stop

       else
          call reset_field_accumulator_to_zero( sim%E_accumulator )
          call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
          some_val = 0.0_f64
          !$omp parallel PRIVATE(Ex,Ey,Ex1,Ey1,tmp5,tmp6)
          !$omp do reduction(+:some_val)
          do i =1,sim%ions_number,2
             SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,particles(i),tmp5,tmp6)
             SLL_INTERPOLATE_FIELD(Ex1,Ey1,accumE,particles(i+1),tmp5,tmp6)
             some_val = some_val + (particles(i)%vx + 0.5_f64 * dt_q_over_m * Ex)**2 &
                  + (particles(i)%vy + 0.5_f64 * dt_q_over_m * Ey)**2 &
                  + (particles(i+1)%vx + 0.5_f64 * dt_q_over_m * Ex1)**2  &
                  + (particles(i+1)%vy + 0.5_f64 * dt_q_over_m * Ey1)**2
          enddo
          !$omp end do
          !$omp end parallel
          call electric_energy( tot_ee, sim%E1, sim%E2, sim%m2d%num_cells1, &
               sim%m2d%num_cells2, sim%m2d%eta1_max, sim%m2d%eta2_max )
          diag_TOTenergy(it) = some_val* 0.5_f64 *(sim%m2d%eta1_max - sim%m2d%eta1_min) * &
               (sim%m2d%eta2_max - sim%m2d%eta2_min)/( sim%world_size*sim%ions_number)  &
               + tot_ee * 0.5_f64
       endif

       !! -- --  diagnostics [end]  -- --

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





!!$  ! Note that the following function has no local variables, which is silly...
!!$  ! This just happened since the guts of the unit test were transplanted here
!!$  ! directly, but this should be cleaned up.
!!$  subroutine run_4d_pic_cart(sim)
!!$    class(sll_pic_simulation_4d_cartesian), intent(inout) :: sim
!!$
!!$
!!$
!!$  end subroutine run_4d_pic_cart


  subroutine delete_4d_pic_cart( sim )
    type(sll_pic_simulation_4d_cartesian) :: sim
  end subroutine delete_4d_pic_cart


  function in_bounds( x, y, mesh ) result(res)
    logical :: res
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    type(sll_cartesian_mesh_2d), pointer :: mesh

    res = (x >= mesh%eta1_min) .and. (x <= mesh%eta1_max) .and. &
          (y >= mesh%eta2_min) .and. (y <= mesh%eta2_max)
  end function in_bounds

  subroutine apply_periodic_bc( mesh, x, y )
    type(sll_cartesian_mesh_2d), pointer :: mesh
    sll_real64, intent(inout) :: x
    sll_real64, intent(inout) :: y

    x = modulo(x,mesh%eta1_max - mesh%eta1_min)! Ca marche quand xmin=ymin=0
    y = modulo(y,mesh%eta2_max - mesh%eta2_min)! Sinon, il faut modulo(...) + xmin/mesh%delta_x
  end subroutine apply_periodic_bc

  subroutine normL2_field_Ex (lee,ee,nx,ny,e,dx,dy)
    sll_real64, intent(out) :: lee, ee
    sll_real64, intent(in) :: dx,dy
    sll_int32, intent(in) :: nx,ny
    sll_real64, dimension(1:nx+1,1:ny+1),intent(in) :: e
    sll_int32 :: i,j
    
    lee = 0._f64
    do j=1,ny
       do i=1,nx
          lee = lee + e(i,j)**2!*e(i,j)
       enddo
    enddo
    lee = lee*dx*dy
    ee = lee
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

    
end module sll_pic_simulation_4d_cartesian_module
