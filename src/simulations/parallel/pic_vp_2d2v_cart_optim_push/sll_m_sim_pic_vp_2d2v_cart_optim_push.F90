module sll_m_sim_pic_vp_2d2v_cart_optim_push

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_accumulators.h"
#include "particle_representation.h"

  use sll_m_accumulators, only: &
    sll_t_electric_field_accumulator, &
    sll_t_electric_field_accumulator_cs, &
    sll_t_field_accumulator_cell, &
    sll_t_field_accumulator_cs, &
    sll_f_new_charge_accumulator_2d, &
    sll_f_new_charge_accumulator_2d_cs, &
    sll_f_new_field_accumulator_2d, &
    sll_f_new_field_accumulator_cs_2d, &
    sll_s_reset_charge_accumulator_to_zero, &
    sll_s_reset_charge_accumulator_to_zero_cs, &
    sll_s_reset_field_accumulator_cs_to_zero, &
    sll_s_reset_field_accumulator_to_zero, &
    sll_t_charge_accumulator_2d, &
    sll_t_charge_accumulator_2d_cs, &
    sll_t_charge_accumulator_2d_cs_ptr, &
    sll_t_charge_accumulator_2d_ptr, &
    sll_s_sum_accumulators, &
    sll_s_sum_accumulators_cs

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d

  use sll_m_charge_to_density, only: &
    sll_s_accumulate_field, &
    sll_s_accumulate_field_cs, &
    sll_s_convert_charge_to_rho_2d_per_per, &
    sll_s_convert_charge_to_rho_2d_per_per_cs

  use sll_m_collective, only: &
    sll_o_collective_allreduce, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_v_world_collective

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_particle_group_4d, only: &
    sll_f_new_particle_4d_group, &
    sll_t_particle_group_4d

  use sll_m_particle_initializers_4d, only: &
    sll_s_initial_particles_4d

  use sll_m_particle_representations, only: &
    sll_t_particle_4d, &
    sll_t_particle_4d_guard

  use sll_m_particle_sort, only: &
    sll_f_new_particle_sorter_2d, &
    sll_t_particle_sorter_2d, &
    sll_s_sort_particles_2d

  use sll_m_pic_utilities, only: &
    sll_s_first_charge_accumulation_2d, &
    sll_s_first_charge_accumulation_2d_cs

  use sll_m_poisson_2d_fft, only: &
    sll_f_new_poisson_2d_fft_solver, &
    sll_t_poisson_2d_fft_solver

  use sll_m_sim_base, only: &
    sll_c_simulation_base_class

  use sll_mpi, only: &
    mpi_sum

#ifdef _OPENMP
  use omp_lib, only: &
    omp_get_num_threads, &
    omp_get_thread_num, &
    omp_get_wtime

#endif
  implicit none

  public :: &
    sll_o_delete, &
    sll_t_pic_simulation_4d_cartesian

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_c_simulation_base_class) :: sll_t_pic_simulation_4d_cartesian
     ! Physics/numerical parameters
     sll_real64 :: dt
     sll_int32  :: num_iterations
     sll_real64 :: thermal_speed_parts
     sll_int32  :: parts_number
     sll_int32  :: guard_size
     sll_int32  :: array_size
     sll_real64, dimension(1:6) :: elec_params
     type(sll_t_particle_group_4d),  pointer :: part_group
     type(sll_t_cartesian_mesh_2d),    pointer :: m2d
     type(sll_t_particle_sorter_2d), pointer :: sorter
     type(sll_t_charge_accumulator_2d_ptr), dimension(:), pointer  :: q_accumulator
     type(sll_t_electric_field_accumulator), pointer :: E_accumulator
     logical :: use_cubic_splines
     type(sll_t_charge_accumulator_2d_cs_ptr), dimension(:), pointer  :: q_accumulator_CS
     type(sll_t_electric_field_accumulator_cs), pointer :: E_accumulator_CS
     sll_real64, dimension(:,:), pointer :: rho
     type(sll_t_poisson_2d_fft_solver), pointer :: poisson
     sll_real64, dimension(:,:), pointer :: E1, E2
     sll_int32 :: my_rank
     sll_int32 :: world_size
     sll_int32 :: n_threads
   contains
     procedure, pass(sim) :: init_from_file => init_4d_pic_cartesian
     procedure, pass(sim) :: run => run_4d_pic_cartesian
  end type sll_t_pic_simulation_4d_cartesian

  interface sll_o_delete
     module procedure delete_4d_pic_cart
  end interface sll_o_delete

!!$  interface initialize
!!$     module procedure initialize_4d_qns_general
!!$  end interface initialize

contains

  subroutine init_4d_pic_cartesian( sim, filename )
    intrinsic :: trim
    class(sll_t_pic_simulation_4d_cartesian), intent(inout) :: sim
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
    sll_real64  :: XMIN, XMAX, YMIN, YMAX
    sll_real64  :: KX, KY
    sll_real64  :: er, psi, omega_i, omega_r
    sll_int32, parameter  :: input_file = 99
    sll_int32, dimension(:), allocatable  :: rand_seed
    sll_int32   :: rand_seed_size
    sll_int32  :: thread_id, it
    type(sll_t_particle_group_4d), pointer  :: pa_gr

    namelist /sim_params/ NUM_PARTICLES, GUARD_SIZE, &
                          PARTICLE_ARRAY_SIZE, &
                          THERM_SPEED, dt, number_iterations, &
                          QoverM, ALPHA, UseCubicSplines
    namelist /grid_dims/  NC_X, NC_Y, XMIN, KX, YMIN, YMAX
!!$    namelist /grid_dims/  NC_X, NC_Y, XMIN, KX, KY, YMIN, YMAX! le cas LANDAU 2D en X
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

    sim%world_size = sll_f_get_collective_size(sll_v_world_collective)
    sim%my_rank    = sll_f_get_collective_rank(sll_v_world_collective)


    XMAX = (2._f64*sll_p_pi/KX)
!!$    YMAX = (2._f64*sll_p_pi/KY)! le cas LANDAU 2D en X
    sim%use_cubic_splines = UseCubicSplines
    sim%thermal_speed_parts = THERM_SPEED
    sim%parts_number = NUM_PARTICLES
    sim%guard_size = GUARD_SIZE  
    sim%array_size = PARTICLE_ARRAY_SIZE
    sim%dt = dt
    sim%num_iterations = number_iterations

    sim%elec_params = (/KX, ALPHA, er, psi, omega_r, omega_i /)
    
    sim%m2d =>  sll_f_new_cartesian_mesh_2d( NC_X, NC_Y, &
                XMIN, XMAX, YMIN, YMAX )

    sim%part_group => sll_f_new_particle_4d_group( &
         NUM_PARTICLES, &
         PARTICLE_ARRAY_SIZE, &
         GUARD_SIZE, &
         QoverM,     &
         sim%m2d )

    sim%sorter => sll_f_new_particle_sorter_2d( sim%m2d )

    sim%poisson => sll_f_new_poisson_2d_fft_solver( sim%m2d%eta1_min,    &
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
!!$    call sll_initial_particles_4d_L2d( sim%thermal_speed_parts, & 
!!$                                   ALPHA, KX, KY, sim%m2d,     &
    call sll_s_initial_particles_4d( sim%thermal_speed_parts, & 
                                   ALPHA, KX, sim%m2d,     &
                                   sim%parts_number,        &
                                   pa_gr, &
                                   rand_seed, sim%my_rank, &
                                   sim%world_size ) 
    SLL_DEALLOCATE_ARRAY( rand_seed, ierr )

    !$omp parallel
    !$omp do
    do j=1,sim%parts_number
       sim%part_group%p_list(j) = pa_gr%p_list(j)
    enddo
    !$omp end do

    !$omp single
    call sll_s_sort_particles_2d( sim%sorter, sim%part_group )
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
       SLL_ALLOCATE(sim%q_accumulator_CS(1:sim%n_threads), ierr)
       thread_id = 0
       !$omp parallel PRIVATE(thread_id)
#ifdef _OPENMP
       thread_id = OMP_GET_THREAD_NUM()
#endif 
       sim%q_accumulator_CS(thread_id+1)%q => sll_f_new_charge_accumulator_2d_cs( sim%m2d )       
       !$omp end parallel
       sim%E_accumulator_CS => sll_f_new_field_accumulator_cs_2d( sim%m2d )           
       call sll_s_first_charge_accumulation_2d_cs( sim%part_group, sim%q_accumulator_CS )

    else
       
       SLL_ALLOCATE(sim%q_accumulator(1:sim%n_threads), ierr)
       thread_id = 0
       !$omp parallel PRIVATE(thread_id)
#ifdef _OPENMP
       thread_id = OMP_GET_THREAD_NUM()
#endif 
       sim%q_accumulator(thread_id+1)%q => sll_f_new_charge_accumulator_2d( sim%m2d )
       !$omp end parallel
       sim%E_accumulator => sll_f_new_field_accumulator_2d( sim%m2d )
       call sll_s_first_charge_accumulation_2d( sim%part_group, sim%q_accumulator)!(1)%q )
    endif
    
  end subroutine init_4d_pic_cartesian

  subroutine run_4d_pic_cartesian( sim )
    intrinsic :: trim
!!!   calls of routines with '_CS' mean use of Cubic Splines  !!!
!     for deposition or interpolation step                      !
    class(sll_t_pic_simulation_4d_cartesian), intent(inout)  :: sim
    sll_int32  :: ierr, it, jj, counter
    sll_int32  :: i, j
    sll_real64 :: tmp5, tmp6, temp
    sll_real64 :: ttmp(1:4,1:2)
    sll_real64, dimension(:,:), pointer :: phi
    sll_int32  :: ncx, ncy, ic_x,ic_y
    sll_int32  :: ic_x1,ic_y1
    sll_real64 :: xmin, ymin
    sll_real64 :: dtrdx, dtrdy
    sll_int32  :: gi ! counter index for guard list
    sll_real64 :: Ex, Ey, Ex1, Ey1
    sll_real64 :: dtqom! dt * qoverm
    sll_real64 :: x, x1 ! for global position
    sll_real64 :: y, y1 ! for global position
    sll_real64 :: dt, ttime
    sll_real64 :: pp_vx, pp_vy
    type(sll_t_particle_4d), dimension(:), pointer :: p
    type(sll_t_field_accumulator_cell), dimension(:), pointer :: accumE
    type(sll_t_field_accumulator_cs), dimension(:), pointer :: accumE_CS
    type(sll_t_particle_4d_guard), dimension(:), pointer :: p_guard
    sll_real64, dimension(:,:), allocatable  ::  diag_energy!, diag_energyL2D
    sll_real64, dimension(:),   allocatable  ::  diag_TOTmoment
    sll_real64, dimension(:),   allocatable  ::  diag_TOTenergy
    sll_real64, dimension(:,:), allocatable :: diag_AccMem! a memory buffer
    sll_real64  ::  t2, t3!  type(sll_time_mark)  :: t2, t3
    sll_real64, dimension(:), allocatable :: rho1d_send
    sll_real64, dimension(:), allocatable :: rho1d_receive
    sll_real64   :: t_init, t_fin, time
    sll_int32 :: save_nb
    sll_int32 :: thread_id
    sll_int32 :: n_threads
    type(sll_t_charge_accumulator_2d),    pointer :: q_accum
    type(sll_t_charge_accumulator_2d_cs), pointer :: q_accum_CS
    sll_int32  :: sort_nb
    sll_real64 :: some_val, une_cst
    sll_real64 :: val_lee, exval_ee
    sll_real64 :: tot_ee, val_ee
    sll_real64 :: omega_i, omega_r, psi
    character(len=2) :: nom_thnb
    sll_real64 :: bors
    
    ncx = sim%m2d%num_cells1
    ncy = sim%m2d%num_cells2
    n_threads = sim%n_threads
    thread_id = 0
    save_nb = sim%num_iterations/2
    if (sim%world_size > 1 ) then
       SLL_ALLOCATE( rho1d_send(1:(ncx+1)*(ncy+1)),    ierr)
       SLL_ALLOCATE( rho1d_receive(1:(ncx+1)*(ncy+1)), ierr)
    endif
    SLL_ALLOCATE(sim%rho(1:ncx+1,1:ncy+1),ierr)
    SLL_ALLOCATE( sim%E1(1:ncx+1,1:ncy+1), ierr )
    SLL_ALLOCATE( sim%E2(1:ncx+1,1:ncy+1), ierr )
    SLL_ALLOCATE(phi(1:ncx+1, 1:ncy+1), ierr)
    SLL_ALLOCATE(diag_energy(1:save_nb, 1:5), ierr)
!    SLL_ALLOCATE(diag_energyL2D(1:save_nb, 1:3), ierr)
!!$    SLL_ALLOCATE(diag_TOTmoment(1:save_nb), ierr)
!!$    SLL_ALLOCATE(diag_TOTenergy(0:sim%num_iterations-1), ierr)
!    SLL_ALLOCATE(diag_AccMem(0:sim%num_iterations-1, 1:2), ierr)

    sort_nb = 10
    p => sim%part_group%p_list
!!!    p_guard => sim%part_group%p_guard
    dt = sim%dt
    dtqom = dt*sim%part_group%qoverm
    xmin = sim%m2d%eta1_min
    ymin = sim%m2d%eta2_min
    dtrdx = dt/sim%m2d%delta_eta1
    dtrdy = dt/sim%m2d%delta_eta2

    if (sim%use_cubic_splines) then
       accumE_CS => sim%E_accumulator_CS%e_acc
       call sll_s_sum_accumulators_cs( sim%q_accumulator_CS, n_threads, ncx*ncy )
       call sll_s_convert_charge_to_rho_2d_per_per_cs( sim%q_accumulator_CS(1)%q, sim%rho ) 
    else
       accumE => sim%E_accumulator%e_acc
       call sll_s_sum_accumulators( sim%q_accumulator, n_threads, ncx*ncy )
       call sll_s_convert_charge_to_rho_2d_per_per( sim%q_accumulator(1)%q, sim%rho ) 
    endif

    if (sim%world_size > 1 ) then
       do j = 1, ncy+1
          do i = 1, ncx+1
             rho1d_send(i+(j-1)*(ncx+1)) = sim%rho(i, j)
             rho1d_receive(i+(j-1)*(ncx+1)) = 0._f64
          enddo
       enddo
       call sll_o_collective_allreduce( sll_v_world_collective, rho1d_send, (ncx+1)*(ncy+1), &
            MPI_SUM, rho1d_receive   )
       do j = 1, ncy+1
          do i = 1, ncx+1
             sim%rho(i, j) = rho1d_receive(i+(j-1)*(ncx+1))
          enddo
       enddo
    endif
!!$    it = 0
!!$    if (sim%my_rank == 0) then
!!$    call sll_gnuplot_2d(sim%m2d%eta1_min, sim%m2d%eta1_max, &
!!$            sim%m2d%num_cells1+1, sim%m2d%eta2_min, &
!!$            sim%m2d%eta2_max, sim%m2d%num_cells2+1, &
!!$            sim%rho, 'rho_init', it, ierr )
!!$    endif

            !!!     REMEMBER: POISSON changes RHO    !!!!
    call sim%poisson%compute_E_from_rho(sim%E1, sim%E2, -sim%rho)

!!$    call electric_energy( tot_ee, sim%E1, sim%E2, sim%m2d%num_cells1, &
!!$         sim%m2d%num_cells2, sim%m2d%eta1_max, sim%m2d%eta2_max )
!!$    bors = 0.0_f64
!!$    !$omp parallel
!!$    !$omp do reduction(+:bors)
!!$    do i =1, sim%parts_number,2
!!$       bors = bors + p(i)%vx**2 + p(i)%vy**2 &
!!$            + p(i+1)%vx**2 + p(i+1)%vy**2
!!$    enddo
!!$    !$omp end do
!!$    !$omp end parallel
!!$    diag_TOTenergy(0) = bors * 0.5_f64*(sim%m2d%eta1_max - sim%m2d%eta1_min)  &
!!$         * (sim%m2d%eta2_max - sim%m2d%eta2_min)/( sim%world_size*sim%parts_number)  &
!!$         + tot_ee * 0.5_f64! ATTENTION : resultats faux avec MPI 
    
    if (sim%use_cubic_splines) then 
       call sll_s_reset_field_accumulator_cs_to_zero( sim%E_accumulator_CS )
       call sll_s_accumulate_field_cs( sim%E1, sim%E2, sim%E_accumulator_CS )
       !$omp parallel do PRIVATE (pp_vx, pp_vy, Ex, Ey, ttmp)
       !$&omp FIRSTPRIVATE(dtqom)
       do i = 1, sim%parts_number
          pp_vx = p(i)%vx
          pp_vy = p(i)%vy
          SLL_INTERPOLATE_FIELD_CS(Ex,Ey,accumE_CS,p(i),ttmp)
          p(i)%vx = pp_vx - 0.5_f64 * dtqom * Ex
          p(i)%vy = pp_vy - 0.5_f64 * dtqom * Ey
       end do! half-step advection of the velocities by -dt/2 here
       !$omp end parallel do

    else
       
       call sll_s_reset_field_accumulator_to_zero( sim%E_accumulator )
       call sll_s_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
       !$omp parallel do PRIVATE (pp_vx, pp_vy, Ex, Ey, tmp5, tmp6)
       !$&omp FIRSTPRIVATE(dtqom)
       do i = 1, sim%parts_number
          pp_vx = p(i)%vx
          pp_vy = p(i)%vy
          SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,p(i),tmp5,tmp6)
          p(i)%vx = pp_vx - 0.5_f64 * dtqom * Ex
          p(i)%vy = pp_vy - 0.5_f64 * dtqom * Ey
       enddo
       !$omp end parallel do   
    endif

    une_cst = 16._f64*(sim%elec_params(2))**2 *(sim%elec_params(3))**2 &
         * sll_p_pi/sim%elec_params(1) !pi/KX*(4._f64*ALPHA*er)**2
    omega_i = sim%elec_params(6)
    omega_r = sim%elec_params(5)
    psi = sim%elec_params(4)

    if (sim%my_rank ==0) then
       write(nom_thnb,'(i2.2)') n_threads
       nom_thnb = trim(nom_thnb)
       open(65,file='logE_OptPush_th'//nom_thnb//'.dat')
    endif
#ifdef _OPENMP
    t2 = omp_get_wtime()!!  call sll_set_time_mark(t2)    
#endif

!  -------------------------
!  ------  TIME LOOP  ------
!  -------------------------
    do it = 0, sim%num_iterations-1
       if (sim%my_rank == 0) then
          exval_ee = une_cst * exp(2._f64*omega_i*real(it,f64)*sim%dt) * &
               ( 0.5_f64 + 0.5_f64*cos(2._f64*(omega_r*real(it,f64)*sim%dt-psi)) )
!!          exval_ee = log( abs(exval_ee) )
          call normL2_field_Ex ( val_lee, val_ee, ncx, ncy, &
                                 sim%E1,  &
                                 sim%m2d%delta_eta1, sim%m2d%delta_eta2 )
!!          call electric_energy( val_ee, sim%E1, sim%E2, ncx, ncy, &
!!                                sim%m2d%delta_eta1, sim%m2d%delta_eta2 )

          counter = 1 + mod(it,save_nb)
!!          diag_energyL2D(counter,:) = (/it*sim%dt, val_ee, log(sqrt(val_ee)) /)
          diag_energy(counter,:) = (/ it*sim%dt, val_lee, 0.5_f64*log(exval_ee), &
                                      val_ee, exval_ee /)

          if (mod(it+1,save_nb)==0) then
             do jj=1,save_nb
!                write(65,*) diag_energyL2D(jj,:)
                write(65,*) diag_energy(jj,:)
             enddo
          endif
       endif

       if (mod(it+1,sort_nb)==0) then 
          call sll_s_sort_particles_2d( sim%sorter, sim%part_group )
       endif

       ! *******************************************************************
       !
       !                   ---- PUSH PARTICLES ----
       !
       ! *******************************************************************
!       some_val = 0.0_f64
#ifdef _OPENMP
       t3 = omp_get_wtime()!!  call sll_set_time_mark(t3)
#endif
       if (sim%use_cubic_splines) then 

          !$omp parallel PRIVATE(x,y,x1,y1,Ex,Ey,gi,temp,ttmp,ic_x,ic_y,thread_id,p_guard,q_accum_CS)
          !$&omp FIRSTPRIVATE(dtqom,ncx,ncy,xmin,ymin,dtrdx,dtrdy)
#ifdef _OPENMP
          thread_id = OMP_GET_THREAD_NUM()
#endif
          call sll_s_reset_charge_accumulator_to_zero_cs( sim%q_accumulator_CS(thread_id+1)%q )
          q_accum_CS => sim%q_accumulator_CS(thread_id+1)%q
          p_guard => sim%part_group%p_guard(thread_id+1)%g_list
          gi = 0
          !$omp do!!  reduction(+:some_val)
          do i = 1, sim%parts_number,2
             SLL_INTERPOLATE_FIELD_CS(Ex,Ey,accumE_CS,p(i),ttmp)
             p(i)%vx = p(i)%vx + dtqom * Ex
             p(i)%vy = p(i)%vy + dtqom * Ey
             SLL_INTERPOLATE_FIELD_CS(Ex,Ey,accumE_CS,p(i+1),ttmp)
             p(i+1)%vx = p(i+1)%vx + dtqom * Ex
             p(i+1)%vy = p(i+1)%vy + dtqom * Ey
             !some_val = some_val + p(i)%vx*p(i)%vx + p(i)%vy*p(i)%vy &
             !              + p(i+1)%vx*p(i+1)%vx + p(i+1)%vy*p(i+1)%vy
             ic_x  = mod( p(i)%ic-1, ncx  )
             ic_y  = int( (p(i)%ic-1)/ncx )
             x = real(ic_x, f64) + p(i)%dx + dtrdx * p(i)%vx! the PUSH
             y = real(ic_y, f64) + p(i)%dy + dtrdy * p(i)%vy
             !
             ic_x  = mod( p(i+1)%ic-1, ncx)
             ic_y  = int( (p(i+1)%ic-1)/ncx)
             x1 = real(ic_x, f64) + p(i+1)%dx + dtrdx * p(i+1)%vx! the PUSH
             y1 = real(ic_y, f64) + p(i+1)%dy + dtrdy * p(i+1)%vy
!
             if (in_bounds(x, y, ncx, ncy)) then ! finish push
                ic_x = int(x)
                ic_y = int(y)
                p(i)%dx = x - real(ic_x, f64)
                p(i)%dy = y - real(ic_y, f64)
                p(i)%ic = ic_x + 1 + ic_y*ncx
                SLL_ACCUMULATE_PARTICLE_CHARGE_CS(q_accum_CS,p(i),ttmp,temp)
             else ! store reference for later processing
                gi = gi + 1
                p_guard(gi)%p => p(i)
             end if
             if (in_bounds( x1, y1, ncx, ncy )) then ! finish push
                ic_x = int(x1)
                ic_y = int(y1)
                p(i+1)%dx = x1 - real(ic_x, f64)
                p(i+1)%dy = y1 - real(ic_y, f64)
                p(i+1)%ic = ic_x + 1 + ic_y * ncx
                SLL_ACCUMULATE_PARTICLE_CHARGE_CS(q_accum_CS,p(i+1),ttmp,temp)
             else ! store reference for later processing
                gi = gi + 1
                p_guard(gi)%p => p(i+1)
             end if
          enddo
          !$omp end do
          sim%part_group%num_postprocess_particles(thread_id+1) = gi
          !$omp end parallel
!       diag_TOTenergy(mod(counter,save_nb)) = some_val

!         ttime = omp_get_wtime()!! ttime = sll_time_elapsed_since(t3)
!          diag_AccMem(it,:) = (/ (it+1)*dt, (32*sim%parts_number*2 + gi*2*8 + &
!               2*128*ncx*ncy + 2*128*ncx*ncy)/ttime/1e9 /)! access to memory in GB/sec
!!$            2*sizeof(sim%q_accumulator_CS%q_acc) + sizeof(sim%E_accumulator_CS%e_acc))

 
          ! Process the particles in the guard list. In the periodic case, no
          ! destruction of particles is needed, so this is simple.
          
          !$omp parallel PRIVATE(x,y,ic_x,ic_y,ttmp,temp,p_guard,q_accum_CS,thread_id,i)
          !$&omp FIRSTPRIVATE(ncx,ncy,xmin,ymin,dtrdx,dtrdy)
#ifdef _OPENMP
          thread_id = OMP_GET_THREAD_NUM()
#endif
          q_accum_CS => sim%q_accumulator_CS(thread_id+1)%q
          p_guard => sim%part_group%p_guard(thread_id+1)%g_list
!          p => sim%part_group%p_list
          do i=1, sim%part_group%num_postprocess_particles(thread_id+1)
             ic_x  = mod( p_guard(i)%p%ic-1, ncx )
             ic_y  = int( (p_guard(i)%p%ic-1)/ncx)
             x = real(ic_x, f64) + p_guard(i)%p%dx + dtrdx * p_guard(i)%p%vx! the push
             y = real(ic_y, f64) + p_guard(i)%p%dy + dtrdy * p_guard(i)%p%vy    
             call apply_periodic_bc( x, y, ncx, ncy)
             ic_x = int(x)
             ic_y = int(y)
             p_guard(i)%p%dx = x - real(ic_x, f64)
             p_guard(i)%p%dy = y - real(ic_y, f64)
             p_guard(i)%p%ic = ic_x + 1 + ic_y*ncx
             SLL_ACCUMULATE_PARTICLE_CHARGE_CS(q_accum_CS,p_guard(i)%p,ttmp,temp)
          end do
          !$omp end parallel
          !       ! reset any counters
          gi = 0
          
       else 

          !$omp parallel PRIVATE(i,x,y,Ex,Ey,gi,tmp5,tmp6,ic_x,ic_y,thread_id,p_guard,q_accum)
          !$&omp FIRSTPRIVATE(dtqom,ncx,ncy,xmin,ymin,dtrdx,dtrdy,accumE)
#ifdef _OPENMP
          thread_id = OMP_GET_THREAD_NUM()
#endif
          call sll_s_reset_charge_accumulator_to_zero ( sim%q_accumulator(thread_id+1)%q )
          q_accum => sim%q_accumulator(thread_id+1)%q
          p_guard => sim%part_group%p_guard(thread_id+1)%g_list
          gi = 0
          !$omp do schedule(runtime)
          do i = 1, sim%parts_number!!$,2
             SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,p(i),tmp5,tmp6)
             p(i)%vx = p(i)%vx + dtqom * Ex
             p(i)%vy = p(i)%vy + dtqom * Ey
!!$             SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,p(i+1),tmp5,tmp6)
!!$             p(i+1)%vx = p(i+1)%vx + dtqom * Ex
!!$             p(i+1)%vy = p(i+1)%vy + dtqom * Ey
             ic_x  = mod( p(i)%ic-1, ncx  )
             ic_y  = int( (p(i)%ic-1)/ncx )
             x = real(ic_x, f64) + p(i)%dx + dtrdx * p(i)%vx! the PUSH
             y = real(ic_y, f64) + p(i)%dy + dtrdy * p(i)%vy
             !
!!$             ic_x  = mod( p(i+1)%ic-1, ncx)
!!$             ic_y  = int( (p(i+1)%ic-1)/ncx)
!!$             x1 = real(ic_x, f64) + p(i+1)%dx + dtrdx * p(i+1)%vx! the PUSH
!!$             y1 = real(ic_y, f64) + p(i+1)%dy + dtrdy * p(i+1)%vy
             !
             if (in_bounds(x, y, ncx, ncy)) then ! finish push
                ic_x = int(x)
                ic_y = int(y)
                p(i)%dx = x - real(ic_x, f64)
                p(i)%dy = y - real(ic_y, f64)
                p(i)%ic = ic_x + 1 + ic_y*ncx
                SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p(i),tmp5,tmp6)
             else ! store reference for later processing
                gi = gi + 1
                p_guard(gi)%p => p(i)
             end if
!!$             if (in_bounds( x1, y1, ncx, ncy )) then ! finish push
!!$                ic_x = int(x1)
!!$                ic_y = int(y1)
!!$                p(i+1)%dx = x1 - real(ic_x, f64)
!!$                p(i+1)%dy = y1 - real(ic_y, f64)
!!$                p(i+1)%ic = ic_x + 1 + ic_y * ncx
!!$                SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p(i+1),tmp5,tmp6)
!!$             else ! store reference for later processing
!!$                gi = gi + 1
!!$                p_guard(gi)%p => p(i+1)
!!$             end if
          enddo
          !$omp end do
          sim%part_group%num_postprocess_particles(thread_id+1) = gi

          ! Process the particles in the guard list. In the periodic case, no
          ! destruction of particles is needed, so this is simple.

!!          !$omp parallel PRIVATE(i,x,y,ic_x,ic_y,tmp5,tmp6,p_guard,q_accum,thread_id)
!!          !$&omp FIRSTPRIVATE(ncx,ncy,xmin,ymin,dtrdx,dtrdy)
!!#ifdef _OPENMP
!!          thread_id = OMP_GET_THREAD_NUM()
!!#endif
          q_accum => sim%q_accumulator(thread_id+1)%q
          p_guard => sim%part_group%p_guard(thread_id+1)%g_list
          do i=1, sim%part_group%num_postprocess_particles(thread_id+1)
             ic_x  = mod( p_guard(i)%p%ic-1, ncx )
             ic_y  = int( (p_guard(i)%p%ic-1)/ncx)
             x = real(ic_x, f64) + p_guard(i)%p%dx + dtrdx * p_guard(i)%p%vx! the push
             y = real(ic_y, f64) + p_guard(i)%p%dy + dtrdy * p_guard(i)%p%vy    
             call apply_periodic_bc( x, y, ncx, ncy)
             ic_x = int(x)
             ic_y = int(y)
             p_guard(i)%p%dx = x - real(ic_x, f64)
             p_guard(i)%p%dy = y - real(ic_y, f64)
             p_guard(i)%p%ic = ic_x + 1 + ic_y*ncx
             SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p_guard(i)%p,tmp5,tmp6)
          end do
          !$omp end parallel
!!$#ifdef _OPENMP
!!$          ttime = omp_get_wtime()!! ttime = sll_time_elapsed_since(t3)
!!$          diag_AccMem(it,:) = (/ (it+1)*dt, &
!!$               (real(32*sim%parts_number*2,f64)+real(gi*2*8,f64) + real(4*32*ncx*ncy,f64) ) &
!!$               /(ttime-t3)/1e9 /)! access to memory in GB/sec
!!$!            2*sizeof(sim%q_accumulator%q_acc) + sizeof(sim%E_accumulator%e_acc))
!!$#endif
          !       ! reset any counters
          gi = 0
       endif
       ! ----- END PUSH PARTICLES -----
       ! -----------------------------

       if (sim%use_cubic_splines) then
          call sll_s_sum_accumulators_cs( sim%q_accumulator_CS, n_threads, ncx*ncy )
          call sll_s_convert_charge_to_rho_2d_per_per_cs( sim%q_accumulator_CS(1)%q, sim%rho )
       else
          call sll_s_sum_accumulators( sim%q_accumulator, n_threads, ncx*ncy )
          call sll_s_convert_charge_to_rho_2d_per_per( sim%q_accumulator(1)%q, sim%rho )
       endif
       if (sim%world_size > 1 ) then
          do j = 1, ncy+1
             do i = 1, ncx+1
                rho1d_send(i+(j-1)*(ncx+1)) = sim%rho(i, j)
                rho1d_receive(i+(j-1)*(ncx+1)) = 0._f64
             enddo
          enddo
          call sll_o_collective_allreduce( sll_v_world_collective, rho1d_send, (ncx+1)*(ncy+1), &
               MPI_SUM, rho1d_receive   )
          do j = 1, ncy+1
             do i = 1, ncx+1
                sim%rho(i, j) = rho1d_receive(i+(j-1)*(ncx+1))
             enddo
          enddo
       endif

       call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, -sim%rho )
       
       if (sim%use_cubic_splines) then
          call sll_s_reset_field_accumulator_cs_to_zero( sim%E_accumulator_CS )
          call sll_s_accumulate_field_cs( sim%E1, sim%E2, sim%E_accumulator_CS )
       else
          call sll_s_reset_field_accumulator_to_zero( sim%E_accumulator )
          call sll_s_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
!!$          some_val = 0.0_f64
!!$          !$omp parallel PRIVATE(Ex,Ey,Ex1,Ey1,tmp5,tmp6)
!!$          !$omp do reduction(+:some_val)
!!$          do i =1,sim%parts_number,2
!!$             SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,p(i),tmp5,tmp6)
!!$             SLL_INTERPOLATE_FIELD(Ex1,Ey1,accumE,p(i+1),tmp5,tmp6)
!!$             some_val = some_val + (p(i)%vx + 0.5_f64 * dtqom * Ex)**2 &
!!$                  + (p(i)%vy + 0.5_f64 * dtqom * Ey)**2 &
!!$                  + (p(i+1)%vx + 0.5_f64 * dtqom * Ex1)**2  &
!!$                  + (p(i+1)%vy + 0.5_f64 * dtqom * Ey1)**2
!!$          enddo
!!$          !$omp end do
!!$          !$omp end parallel
!!$          call electric_energy( tot_ee, sim%E1, sim%E2, sim%m2d%num_cells1, &
!!$               sim%m2d%num_cells2, sim%m2d%eta1_max, sim%m2d%eta2_max )
!!$          diag_TOTenergy(it) = some_val* 0.5_f64 *(sim%m2d%eta1_max - sim%m2d%eta1_min) * &
!!$               (sim%m2d%eta2_max - sim%m2d%eta2_min)/( sim%world_size*sim%parts_number)  &
!!$               + tot_ee * 0.5_f64
       endif

    enddo
!  ---  ---  - - -   END TIME LOOP  - - -  --- -----
#ifdef _OPENMP
    time = omp_get_wtime()!! time = sll_time_elapsed_since(t2)
    if (sim%my_rank ==0) then 
       open(93,file='time_OptPush_omp_WITH_SCHEDULE.dat',position='append')
       write(93,*) '# Nb of threads   ||   time (sec)   ||   average pushes/sec'
       write(93,*) n_threads, time-t2, int(sim%num_iterations,i64)*int(sim%parts_number,i64)/(time-t2)
       close(93)
    endif
#endif

    if (sim%my_rank ==0) then 
       close(65)
!!$       open(65,file='Energie_totale_OptPush.dat')
!!$       do it=0,sim%num_iterations-1
!!$          write(65,*) it*dt,  diag_TOTenergy(it), (diag_TOTenergy(it)-diag_TOTenergy(0))/diag_TOTenergy(0)
!!$       enddo
!!$       close(65)
    endif

!!$#ifdef _OPENMP
!!$    if (sim%my_rank==0) then! .and. thread_id==0) then
!!$       write(nom_thnb,'(i2.2)') n_threads
!!$       nom_thnb = trim(nom_thnb)
!!$       open(65,file='AccesstoMemory_rk0_th'//nom_thnb//'.dat')
!!$!    !!$    if (sim%my_rank==1) open(66,file='AccesstoMemory_rk1.dat')
!!$       write(65,*) diag_AccMem(0,:), diag_AccMem(0,2)
!!$       do jj = 1, sim%num_iterations-1
!!$          if ( mod(jj+1,sort_nb) == 0 ) then
!!$             write(65,*) diag_AccMem(jj,:), diag_AccMem(jj,2)
!!$!    !!$          if (sim%my_rank==1) write(66,*) diag_AccMem(jj,:), diag_AccMem(jj,2)
!!$          else
!!$             write(65,*) diag_AccMem(jj,:)
!!$!    !!$          if (sim%my_rank==1) write(66,*) diag_AccMem(jj,:)
!!$          endif
!!$       enddo
!!$       close(65) 
!!$!    !!$       close(66)
!!$    endif
!!$#endif

    if (sim%world_size > 1 ) then
       SLL_DEALLOCATE_ARRAY(rho1d_send, ierr)
       SLL_DEALLOCATE_ARRAY(rho1d_receive, ierr)
    endif

    SLL_DEALLOCATE(sim%rho,   ierr)
    SLL_DEALLOCATE(sim%E1,    ierr)
    SLL_DEALLOCATE(sim%E2,    ierr)
    SLL_DEALLOCATE(phi, ierr)
    SLL_DEALLOCATE_ARRAY(diag_energy, ierr)
!    SLL_DEALLOCATE_ARRAY(diag_energyL2D, ierr)
!!$    SLL_DEALLOCATE_ARRAY(diag_TOTenergy, ierr)
!!$    SLL_DEALLOCATE_ARRAY(diag_TOTmoment, ierr)
!    SLL_DEALLOCATE_ARRAY(diag_AccMem, ierr)
  end subroutine run_4d_pic_cartesian





!!$  ! Note that the following function has no local variables, which is silly...
!!$  ! This just happened since the guts of the unit test were transplanted here
!!$  ! directly, but this should be cleaned up.
!!$  subroutine run_4d_pic_cart(sim)
!!$    class(sll_t_pic_simulation_4d_cartesian), intent(inout) :: sim
!!$
!!$
!!$
!!$  end subroutine run_4d_pic_cart


  subroutine delete_4d_pic_cart( sim )
    type(sll_t_pic_simulation_4d_cartesian) :: sim
  end subroutine delete_4d_pic_cart


  function in_bounds( x, y, nc_x, nc_y ) result(res)
    logical :: res
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_int32,  intent(in) :: nc_x, nc_y

    res = (x >= 0._f64) .and. (x <= real(nc_x, f64)) .and. &
          (y >= 0._f64) .and. (y <= real(nc_y, f64))
  end function in_bounds

  subroutine apply_periodic_bc( x, y, nc_x, nc_y)
    sll_real64, intent(inout) :: x
    sll_real64, intent(inout) :: y
    sll_int32,  intent(in) :: nc_x, nc_y

    x = modulo( x, real(nc_x, f64) )! Ca marche quand xmin=ymin=0
    y = modulo( y, real(nc_y, f64) )! Sinon, il faut modulo(...) + xmin/mesh%delta_x
  end subroutine apply_periodic_bc

  subroutine normL2_field_Ex (lee,ee,nx,ny,e,dx,dy)
    ! lee = log of the L2-norm of the eletric field
    ! ee  = the electric energy of Ex
    sll_real64, intent(out) ::  lee, ee
    sll_real64, intent(in) :: dx,dy
    sll_int32, intent(in) :: nx,ny
    sll_real64, dimension(1:nx+1,1:ny+1),intent(in) :: e
    sll_int32 :: i,j
    
    lee = 0._f64
    do j=1,ny
       do i=1,nx
          lee = lee + e(i,j)**2
       enddo
    enddo
    lee = lee*dx*dy
    ee  = lee
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
  
end module sll_m_sim_pic_vp_2d2v_cart_optim_push
