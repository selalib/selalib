module sll_m_sim_pic_gc_2d0v_cart_optim

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_accumulators.h"
#include "particle_representation.h"

  use sll_m_accumulators, only: &
    electric_field_accumulator, &
    electric_field_accumulator_cs, &
    field_accumulator_cell, &
    field_accumulator_cs, &
    new_charge_accumulator_2d, &
    new_charge_accumulator_2d_cs, &
    new_field_accumulator_2d, &
    new_field_accumulator_cs_2d, &
    reset_charge_accumulator_to_zero, &
    reset_charge_accumulator_to_zero_cs, &
    sll_charge_accumulator_2d, &
    sll_charge_accumulator_2d_cs, &
    sll_charge_accumulator_2d_cs_ptr, &
    sll_charge_accumulator_2d_ptr, &
    sum_accumulators, &
    sum_accumulators_cs

  use sll_m_cartesian_meshes, only: &
    new_cartesian_mesh_2d, &
    sll_cartesian_mesh_2d

  use sll_m_charge_to_density, only: &
    sll_accumulate_field, &
    sll_accumulate_field_cs, &
    sll_convert_charge_to_rho_2d_per_per, &
    sll_convert_charge_to_rho_2d_per_per_cs

  use sll_m_collective, only: &
    sll_collective_allreduce, &
    sll_get_collective_rank, &
    sll_get_collective_size, &
    sll_world_collective

  use sll_m_constants, only: &
    sll_pi

  use sll_m_gnuplot, only: &
    sll_gnuplot_2d

  use sll_m_particle_group_2d, only: &
    new_particle_2d_group, &
    sll_particle_group_2d

  use sll_m_particle_initializers_2d, only: &
    sll_initial_particles_2d_kh

  use sll_m_particle_representations, only: &
    sll_particle_2d, &
    sll_particle_2d_guard

  use sll_m_particle_sort, only: &
    sll_new_particle_sorter_2d, &
    sll_particle_sorter_2d, &
    sll_sort_gc_particles_2d

  use sll_m_pic_utilities, only: &
    sll_first_gc_charge_accumulation_2d, &
    sll_first_gc_charge_accumulation_2d_cs

  use sll_m_poisson_2d_fft, only: &
    new_poisson_2d_fft_solver, &
    poisson_2d_fft_solver

  use sll_m_sim_base, only: &
    sll_simulation_base_class

  use sll_m_timer, only: &
    sll_set_time_mark, &
    sll_time_mark

  use sll_mpi, only: &
    mpi_sum

#ifdef _OPENMP
  use omp_lib, only: &
    omp_get_num_threads, &
    omp_get_thread_num

#endif
  implicit none

  public :: &
    sll_delete, &
    sll_pic_simulation_2d_gc_cartesian

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  type, extends(sll_simulation_base_class) :: sll_pic_simulation_2d_gc_cartesian
     ! Physics/numerical parameters
     sll_real64 :: dt
     sll_int32  :: num_iterations
     sll_real64 :: thermal_speed_parts
     sll_int32  :: parts_number
     sll_int32  :: guard_size
     sll_int32  :: array_size
     type(sll_particle_group_2d),  pointer :: part_group
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
     procedure, pass(sim) :: init_from_file => init_2d_pic_cartesian
     procedure, pass(sim) :: run => run_2d_pic_cartesian
  end type sll_pic_simulation_2d_gc_cartesian

  interface sll_delete
     module procedure delete_4d_pic_cart
  end interface sll_delete

!!$  interface initialize
!!$     module procedure initialize_4d_qns_general
!!$  end interface initialize

contains

  subroutine init_2d_pic_cartesian( sim, filename )
    intrinsic :: trim
    class(sll_pic_simulation_2d_gc_cartesian), intent(inout) :: sim
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
    sll_int32, parameter  :: input_file = 99
    sll_int32, dimension(:), allocatable  :: rand_seed
    sll_int32   :: rand_seed_size
    sll_int32   :: thread_id
    type(sll_particle_group_2d), pointer  :: pa_gr

    namelist /sim_params/ NUM_PARTICLES, GUARD_SIZE, &
                          PARTICLE_ARRAY_SIZE, &
                          THERM_SPEED, number_iterations, &
                          dt, QoverM, ALPHA, UseCubicSplines
    namelist /grid_dims/  NC_X, NC_Y, XMIN, KX, YMIN
    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_vp4d_par_cart() failed to open file ', filename
       STOP
    end if
    read(input_file, sim_params)
    read(input_file, grid_dims)
    close(input_file)

    sim%world_size = sll_get_collective_size(sll_world_collective)
    sim%my_rank    = sll_get_collective_rank(sll_world_collective)

    XMAX = 2._f64*sll_pi/KX
    YMAX = 2._f64*sll_pi
    sim%use_cubic_splines = UseCubicSplines
    sim%thermal_speed_parts = THERM_SPEED
    sim%parts_number = NUM_PARTICLES
    sim%guard_size = GUARD_SIZE  
    sim%array_size = PARTICLE_ARRAY_SIZE
    sim%num_iterations = number_iterations
    sim%dt = dt
    sim%m2d =>  new_cartesian_mesh_2d( NC_X, NC_Y, &
                XMIN, XMAX, YMIN, YMAX )

    sim%part_group => new_particle_2d_group( &
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
    call sll_initial_particles_2d_KH( ALPHA, KX, sim%m2d,     &
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
    call sll_sort_gc_particles_2d( sim%sorter, sim%part_group )
    sim%n_threads =  1
    !$omp end single

#ifdef _OPENMP
    if (OMP_GET_THREAD_NUM() == 0) then
       sim%n_threads =  OMP_GET_NUM_THREADS()
    endif
#endif
    !$omp end parallel

    if (sim%use_cubic_splines) then
       SLL_ALLOCATE(sim%q_accumulator_CS(1:sim%n_threads), ierr)
       thread_id = 0
       !$omp parallel PRIVATE(thread_id)
#ifdef _OPENMP
       thread_id = OMP_GET_THREAD_NUM()
#endif 
       sim%q_accumulator_CS(thread_id+1)%q => new_charge_accumulator_2d_CS( sim%m2d )       
       !$omp end parallel
       sim%E_accumulator_CS => new_field_accumulator_CS_2d( sim%m2d )           
       call sll_first_gc_charge_accumulation_2d_CS( sim%part_group, sim%q_accumulator_CS )

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
       call sll_first_gc_charge_accumulation_2d( sim%part_group, sim%q_accumulator)!(1)%q )
    endif

  end subroutine init_2d_pic_cartesian


  subroutine run_2d_pic_cartesian( sim )!_RK2
!!!   calls of routines with '_CS' mean use of Cubic Splines  !!!
!     for deposition or interpolation step                      !
    class(sll_pic_simulation_2d_gc_cartesian), intent(inout)  :: sim
    sll_int32  :: ierr, it, jj, counter
    sll_int32  :: i, j
    sll_real64 :: tmp3, tmp4, tmp5, tmp6
    sll_real64 :: tmp7, tmp8
    sll_real64 :: ee_val, enst
    sll_real64, dimension(:,:), pointer :: phi
    sll_int32  :: ncx, ncy, ic_x,ic_y
    sll_int32  :: ic_x1,ic_y1
    sll_real32 :: off_x, off_y,off_x1,off_y1
    sll_real64 :: xmin, ymin, rdx, rdy
    sll_real64 :: xmax, ymax
    sll_int32  :: gi ! counter index for guard list
    sll_real64 :: Ex, Ey, Ex1, Ey1
    sll_real64 :: qoverm
    sll_real64 :: x, x1  ! for global position
    sll_real64 :: y, y1  ! for global position
    sll_real64 :: pp_vx,  pp_vy
    sll_real64 :: dt, tfin
    type(sll_time_mark) :: tinit
    sll_real64 :: temp
    sll_real64 :: temp1(1:4,1:2), temp2(1:4,1:2)
    type(sll_particle_2d), dimension(:), pointer :: p
    type(sll_particle_2d), dimension(:), allocatable :: ploc
    type(field_accumulator_cell), dimension(:), pointer :: accumE
    type(field_accumulator_CS), dimension(:), pointer :: accumE_CS
    type(sll_particle_2d_guard), dimension(:), pointer :: p_guard
    sll_real64, dimension(:,:), allocatable  ::  diag_energy! a memory buffer
    sll_real64, dimension(:),   allocatable  ::  diag_enstrophy
    sll_real64, dimension(:),   allocatable  ::  diag_TOTenergy
!!$    sll_real64, dimension(:,:), allocatable :: diag_AccMem! a memory buffer
    type(sll_time_mark)  :: t2, t3, t8
    sll_real64, dimension(:), allocatable :: rho1d_send
    sll_real64, dimension(:), allocatable :: rho1d_receive
    sll_real64   :: t_init, t_fin, time
    sll_int32 :: save_nb
    sll_int32 :: thread_id
    sll_int32 :: n_threads
    type(sll_charge_accumulator_2d),    pointer :: q_accum
    type(sll_charge_accumulator_2d_CS), pointer :: q_accum_CS
    sll_int32  :: sort_nb
    sll_real64  :: some_val
    sll_real64  :: x2, y2
    character(6)  :: it_name
    character(22) :: nnnom

    ncx = sim%m2d%num_cells1
    ncy = sim%m2d%num_cells2
    n_threads = sim%n_threads
    thread_id = 0
    save_nb = sim%num_iterations/2

    if (sim%world_size > 1 ) then
       SLL_ALLOCATE( rho1d_send(1:(ncx+1)*(ncy+1)),    ierr)
       SLL_ALLOCATE( rho1d_receive(1:(ncx+1)*(ncy+1)), ierr)
    endif
    SLL_ALLOCATE(ploc(1:sim%parts_number),ierr )
    SLL_ALLOCATE(sim%rho(1:ncx+1,1:ncy+1),ierr)
    SLL_ALLOCATE( sim%E1(1:ncx+1,1:ncy+1), ierr )
    SLL_ALLOCATE( sim%E2(1:ncx+1,1:ncy+1), ierr )
    SLL_ALLOCATE(phi(1:ncx+1, 1:ncy+1), ierr)
!!$    SLL_ALLOCATE(diag_energy(1:save_nb, 1:3), ierr)
    SLL_ALLOCATE(diag_TOTenergy(0:sim%num_iterations-1), ierr)
    SLL_ALLOCATE(diag_enstrophy(0:sim%num_iterations), ierr)
!!$    SLL_ALLOCATE(diag_AccMem(0:sim%num_iterations-1, 1:2), ierr)

    sort_nb = 10
    p => sim%part_group%p_list
    qoverm = sim%part_group%qoverm
!!!    p_guard => sim%part_group%p_guard
    dt   = sim%dt
    xmin = sim%m2d%eta1_min
    ymin = sim%m2d%eta2_min
    xmax = sim%m2d%eta1_max
    ymax = sim%m2d%eta2_max

    rdx  = 1._f64/sim%m2d%delta_eta1
    rdy  = 1._f64/sim%m2d%delta_eta2

    if (sim%use_cubic_splines) then
       accumE_CS => sim%E_accumulator_CS%e_acc
       call sum_accumulators_CS( sim%q_accumulator_CS, n_threads, ncx*ncy )
       call sll_convert_charge_to_rho_2d_per_per_CS( sim%q_accumulator_CS(1)%q, sim%rho )
       if (sim%world_size > 1 ) then
       print*, sim%world_size, 'mpi nodes'
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
       endif
       if (sim%my_rank == 0) then
          call normL2_field(enst,ncx,ncy,sim%rho,sim%m2d%delta_eta1,sim%m2d%delta_eta2)
          diag_enstrophy(0) = enst
          it = 1
          call sll_gnuplot_2d(xmin, sim%m2d%eta1_max, ncx+1, ymin, &
               sim%m2d%eta2_max, ncy+1, &
               sim%rho, 'RHO_init_CS', it, ierr )
       endif
          ! POISSON changes rho
       call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, -sim%rho )
       call sll_accumulate_field_CS( sim%E1, sim%E2, sim%E_accumulator_CS )

    else

       accumE => sim%E_accumulator%e_acc
       call sum_accumulators( sim%q_accumulator, n_threads, ncx*ncy )
       call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator(1)%q, sim%rho ) 
       if (sim%world_size > 1 ) then
       print*, sim%world_size, 'mpi nodes'
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
       endif
       if (sim%my_rank == 0) then
          call normL2_field(enst,ncx,ncy,sim%rho,sim%m2d%delta_eta1,sim%m2d%delta_eta2)
          diag_enstrophy(0) = enst
!          open(50,file='rhoGC_1d7parts_INIT.dat')
!          do i = 1, ncx+1
!             do j = 1, ncy+1
!                write(50,*)  sim%m2d%delta_eta1*real(ncx-1,f64), &
!                                sim%m2d%delta_eta2*real(ncy-1,f64), &
!                                sim%rho(i,j)! XMIN=0, YMIN=0 !!!
!             enddo
!             write(50,*)
!          enddo
!          close(50)
          it = 1
          call sll_gnuplot_2d(xmin, sim%m2d%eta1_max, ncx+1, ymin, &
               sim%m2d%eta2_max, ncy+1, &
               sim%rho, 'RHO_init', it, ierr )
       endif
       ! POISSON changes rho
       call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, -sim%rho )
       call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
    endif
    
! -----------------------------
! --------  TIME LOOP  --------
! -----------------------------
!    call sll_set_time_mark( tinit )
    do it = 0, sim%num_iterations-1

       if (sim%world_size == 1 ) then
          call normL2_field_E(ee_val,ncx,ncy, sim%E1, sim%E2, &
               sim%m2d%delta_eta1, sim%m2d%delta_eta2)
          diag_TOTenergy(it) = ee_val
       endif
       
       if (mod(it+1,sort_nb)==0) then 
          call sll_sort_gc_particles_2d( sim%sorter, sim%part_group )
       endif
       ! *******************************************************************
       !
       !                   ---- PUSH PARTICLES ----
       !
       ! *******************************************************************  
       if (sim%use_cubic_splines) then 
          call sll_accumulate_field_CS( sim%E1, sim%E2, sim%E_accumulator_CS )
          !$omp parallel PRIVATE(x,y,Ex,Ey,ic_x,ic_y,off_x,off_y,temp1,temp2,temp,tmp5,tmp6,thread_id,q_accum_CS)
          !$&omp FIRSTPRIVATE(dt,xmin,xmax,ymin,ymax,ncx,rdx,rdy)
#ifdef _OPENMP
          thread_id = OMP_GET_THREAD_NUM()
#endif    
          call reset_charge_accumulator_to_zero_CS ( sim%q_accumulator_CS(thread_id+1)%q )
          q_accum_CS => sim%q_accumulator_CS(thread_id+1)%q
          !$omp do
          do i = 1, sim%parts_number
             SLL_INTERPOLATE_FIELD_CS(Ex,Ey,accumE_CS,p(i),temp1)
             GET_PARTICLE_POSITION(p(i),sim%m2d,x,y)
             ploc(i) = p(i)
             x = modulo(x + 0.5_f64*dt*Ey, xmax - xmin)
             y = modulo(y - 0.5_f64*dt*Ex, ymax - ymin)
             SET_PARTICLE_POSITION(p(i),xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp5,tmp6)
             SLL_ACCUMULATE_PARTICLE_CHARGE_CS(q_accum_CS,p(i),temp2,temp)
          enddo
          !$omp end do   
          !$omp end parallel
          call sum_accumulators_CS( sim%q_accumulator_CS, n_threads, ncx*ncy )
          call sll_convert_charge_to_rho_2d_per_per_CS( sim%q_accumulator_CS(1)%q, sim%rho ) 

       else

          call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
          !$omp parallel PRIVATE(x,y,Ex,Ey,ic_x,ic_y,off_x,off_y,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,thread_id,q_accum)
          !$&omp FIRSTPRIVATE(dt,xmin,xmax,ymin,ymax,ncx,rdx,rdy)
#ifdef _OPENMP
          thread_id = OMP_GET_THREAD_NUM()
#endif    
          call reset_charge_accumulator_to_zero ( sim%q_accumulator(thread_id+1)%q )
          q_accum => sim%q_accumulator(thread_id+1)%q
          !$omp do
          do i = 1, sim%parts_number
             SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,p(i),tmp3,tmp4)
             GET_PARTICLE_POSITION(p(i),sim%m2d,x,y)
             ploc(i) = p(i)
             x = modulo(x + 0.5_f64*dt*Ey, xmax - xmin)
             y = modulo(y - 0.5_f64*dt*Ex, ymax - ymin)
             SET_PARTICLE_POSITION(p(i),xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp5,tmp6)
             SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p(i),tmp7,tmp8)
          enddo
          !$omp end do   
          !$omp end parallel
          call sum_accumulators( sim%q_accumulator, n_threads, ncx*ncy )
          call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator(1)%q, sim%rho ) 
       endif

       if (sim%world_size > 1 ) then
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
       endif
       call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, -sim%rho )

       if (sim%use_cubic_splines) then 
          call sll_accumulate_field_CS( sim%E1, sim%E2, sim%E_accumulator_CS )
          !$omp parallel PRIVATE(x,y,Ex,Ey,ic_x,ic_y,off_x,off_y,temp1,temp2,temp,tmp5,tmp6,thread_id,q_accum_CS)
          !$&omp FIRSTPRIVATE(dt,xmin,xmax,ymin,ymax,ncx,rdx,rdy)
#ifdef _OPENMP
          thread_id = OMP_GET_THREAD_NUM()
#endif    
          call reset_charge_accumulator_to_zero_CS ( sim%q_accumulator_CS(thread_id+1)%q )
          q_accum_CS => sim%q_accumulator_CS(thread_id+1)%q
          !$omp do
          do i = 1, sim%parts_number
             SLL_INTERPOLATE_FIELD_CS(Ex,Ey,accumE_CS,p(i),temp1)
             GET_PARTICLE_POSITION(ploc(i),sim%m2d,x,y)
             x = modulo(x + dt*Ey, xmax - xmin)
             y = modulo(y - dt*Ex, ymax - ymin)
             SET_PARTICLE_POSITION(p(i),xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp5,tmp6)
             SLL_ACCUMULATE_PARTICLE_CHARGE_CS(q_accum_CS,p(i),temp2,temp)
          enddo
          !$omp end do   
          !$omp end parallel
          call sum_accumulators_CS( sim%q_accumulator_CS, n_threads, ncx*ncy )
          call sll_convert_charge_to_rho_2d_per_per_CS( sim%q_accumulator_CS(1)%q, sim%rho ) 

       else

          call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
          !$omp parallel PRIVATE(x,y,Ex,Ey,ic_x,ic_y,off_x,off_y,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,thread_id,q_accum)
          !$&omp FIRSTPRIVATE(dt,xmin,xmax,ymin,ymax,ncx,rdx,rdy)
#ifdef _OPENMP
          thread_id = OMP_GET_THREAD_NUM()
#endif    
          call reset_charge_accumulator_to_zero ( sim%q_accumulator(thread_id+1)%q )
          q_accum => sim%q_accumulator(thread_id+1)%q
          !$omp do
          do i = 1, sim%parts_number
             SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,p(i),tmp3,tmp4)
             GET_PARTICLE_POSITION(ploc(i),sim%m2d,x,y)
             x = modulo(x + dt*Ey, xmax - xmin)
             y = modulo(y - dt*Ex, ymax - ymin)
             SET_PARTICLE_POSITION(p(i),xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp5,tmp6)
             SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p(i),tmp7,tmp8)
          enddo
          !$omp end do   
          !$omp end parallel
          call sum_accumulators( sim%q_accumulator, n_threads, ncx*ncy )
          call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator(1)%q, sim%rho ) 
       endif
       ! ----- END PUSH PARTICLES -----
       ! -----------------------------
       if (sim%world_size > 1 ) then
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
       endif

       if (mod(it+1, 200)==0 .and. sim%my_rank == 0) then
!          tfin = sll_time_elapsed_since(tinit)
!!$          write(it_name,'(i5.5)') it+1
!!$          nnnom = 'parts_at'//trim(adjustl(it_name))//'.dat'
!!$          open(50,file=nnnom)
!!$          do i = 1, sim%parts_number
!!$             GET_PARTICLE_POSITION( p(i),sim%m2d,x,y )
!!$             write(50,*) x, y
!!$          enddo
!!$          close(50)
!          filename = 'rhoGC_1d7parts_'//trim(adjustl(it_name))//'.dat'
!          open(50,file=filename)
!          do i = 1, ncx+1
!             do j = 1, ncy+1
!                write(50,*)  sim%m2d%delta_eta1*real(ncx-1,f64), &
!                                sim%m2d%delta_eta2*real(ncy-1,f64), &
!                                sim%rho(i,j)! XMIN=0, YMIN=0 !!!
!             enddo
!             write(50,*)
!          enddo
!          close(50)
          call sll_gnuplot_2d(xmin, sim%m2d%eta1_max, ncx+1, ymin, &
               sim%m2d%eta2_max, ncy+1, &
               sim%rho, 'rhoGC_1d7p_it', it+1, ierr )
       endif

       if (sim%my_rank == 0) then
          call normL2_field(enst,ncx,ncy,sim%rho,sim%m2d%delta_eta1,sim%m2d%delta_eta2)
          diag_enstrophy(it+1) = enst
       endif

       call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, -sim%rho )

    enddo
    ! ----------   END TIME LOOP   -----------
    ! ----------------------------------------
    if (sim%world_size == 1 ) then
       open(65,file='elecenergy.dat')
       open(66,file='enstrophy.dat')
       do it=0,sim%num_iterations-1
          write(65,*) it*dt,  diag_TOTenergy(it), &
               (diag_TOTenergy(it)-diag_TOTenergy(0))/diag_TOTenergy(0)
          write(66,*) it*dt, diag_enstrophy(it), &
               (diag_enstrophy(it)-diag_enstrophy(0))/diag_enstrophy(0)
       enddo
       close(65) ; close(66)
    endif

    SLL_DEALLOCATE_ARRAY(ploc,ierr)    
    SLL_DEALLOCATE(sim%rho,   ierr)
    SLL_DEALLOCATE(sim%E1,    ierr)
    SLL_DEALLOCATE(sim%E2,    ierr)
    SLL_DEALLOCATE(phi, ierr)
    if (sim%world_size > 1 ) then
       SLL_DEALLOCATE_ARRAY(rho1d_send, ierr)
       SLL_DEALLOCATE_ARRAY(rho1d_receive, ierr)
    endif
    SLL_DEALLOCATE_ARRAY(diag_TOTenergy, ierr)
    SLL_DEALLOCATE_ARRAY(diag_enstrophy, ierr)
!!$    SLL_DEALLOCATE_ARRAY(diag_energy, ierr)

  end subroutine run_2d_pic_cartesian!_RK2

  subroutine run_2d_pic_cartesian_ExpEuler( sim )!( sim )!
!!!   calls of routines with '_CS' mean use of Cubic Splines  !!!
!     for deposition or interpolation step                      !
    class(sll_pic_simulation_2d_gc_cartesian), intent(inout)  :: sim
    sll_int32  :: ierr, it, jj, counter
    sll_int32  :: i, j
    sll_real64 :: tmp3, tmp4, tmp5, tmp6
    sll_real64 :: tmp7, tmp8
    sll_real64, dimension(:,:), pointer :: phi
    sll_int32  :: ncx, ncy, ic_x,ic_y
    sll_int32  :: ic_x1,ic_y1
    sll_real32 :: off_x, off_y,off_x1,off_y1
    sll_real64 :: xmin, ymin, rdx, rdy
    sll_real64 :: xmax, ymax
    sll_int32  :: gi ! counter index for guard list
    sll_real64 :: Ex, Ey, Ex1, Ey1
    sll_real64 :: qoverm
    sll_real64 :: x, x1  ! for global position
    sll_real64 :: y, y1  ! for global position
    sll_real64 :: pp_vx,  pp_vy
    sll_real64 :: dt, tfin
    type(sll_time_mark) :: tinit
    sll_real64 :: temp
    sll_real64 :: temp1(1:4,1:2), temp2(1:4,1:2)
    type(sll_particle_2d), dimension(:), pointer :: p
    type(field_accumulator_cell), dimension(:), pointer :: accumE
    type(field_accumulator_CS), dimension(:), pointer :: accumE_CS
    type(sll_particle_2d_guard), dimension(:), pointer :: p_guard
    sll_real64, dimension(:,:), allocatable  ::  diag_energy! a memory buffer
    sll_real64, dimension(:),   allocatable  ::  diag_TOTenergy
!!$    sll_real64, dimension(:,:), allocatable :: diag_AccMem! a memory buffer
    type(sll_time_mark)  :: t2, t3, t8
    sll_real64, dimension(:), allocatable :: rho1d_send
    sll_real64, dimension(:), allocatable :: rho1d_receive
    sll_real64   :: t_init, t_fin, time
    sll_int32 :: save_nb
    sll_int32 :: thread_id
    sll_int32 :: n_threads
    type(sll_charge_accumulator_2d),    pointer :: q_accum
    type(sll_charge_accumulator_2d_CS), pointer :: q_accum_CS
    sll_int32  :: sort_nb
    sll_real64  :: some_val
    sll_real64  :: x2, y2
    character(6)  :: it_name
    character(22) :: nnnom

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
!!$    SLL_ALLOCATE(diag_energy(1:save_nb, 1:3), ierr)
!!$    SLL_ALLOCATE(diag_TOTenergy(0:save_nb), ierr)
!!$    SLL_ALLOCATE(diag_AccMem(0:sim%num_iterations-1, 1:2), ierr)

    sort_nb = 10
    p => sim%part_group%p_list
    qoverm = sim%part_group%qoverm
!!!    p_guard => sim%part_group%p_guard
    dt   = sim%dt
    xmin = sim%m2d%eta1_min
    ymin = sim%m2d%eta2_min
    xmax = sim%m2d%eta1_max
    ymax = sim%m2d%eta2_max

    rdx  = 1._f64/sim%m2d%delta_eta1
    rdy  = 1._f64/sim%m2d%delta_eta2

    if (sim%use_cubic_splines) then
       accumE_CS => sim%E_accumulator_CS%e_acc
       call sum_accumulators_CS( sim%q_accumulator_CS, n_threads, ncx*ncy )
       call sll_convert_charge_to_rho_2d_per_per_CS( sim%q_accumulator_CS(1)%q, sim%rho )
       if (sim%world_size > 1 ) then
          print*, sim%world_size, 'mpi nodes'
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
       endif
       if (sim%my_rank == 0) then
          it = 1
          call sll_gnuplot_2d(xmin, sim%m2d%eta1_max, ncx+1, ymin, &
               sim%m2d%eta2_max, ncy+1, &
               sim%rho, 'RHO_init_CS', it, ierr )
       endif
       ! POISSON changes rho
       call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, -sim%rho )
       call sll_accumulate_field_CS( sim%E1, sim%E2, sim%E_accumulator_CS )
    else
       accumE => sim%E_accumulator%e_acc
       call sum_accumulators( sim%q_accumulator, n_threads, ncx*ncy )
       call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator(1)%q, sim%rho ) 
       if (sim%world_size > 1 ) then
          print*, sim%world_size, 'mpi nodes'
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
       endif
       if (sim%my_rank == 0) then
          it = 1
          call sll_gnuplot_2d(xmin, sim%m2d%eta1_max, ncx+1, ymin, &
               sim%m2d%eta2_max, ncy+1, &
               sim%rho, 'RHO_init_CS', it, ierr )
       endif
       ! POISSON changes rho
       call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, -sim%rho )
       call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
    endif

! -----------------------------
! --------  TIME LOOP  --------
! -----------------------------
    call sll_set_time_mark( tinit )
    do it = 0, sim%num_iterations-1

       if (mod(it+1,sort_nb)==0) then 
          call sll_sort_gc_particles_2d( sim%sorter, sim%part_group )
       endif
       ! *******************************************************************
       !
       !                   ---- PUSH PARTICLES ----
       !
       ! *******************************************************************  
       if (sim%use_cubic_splines) then 
          call sll_accumulate_field_CS( sim%E1, sim%E2, sim%E_accumulator_CS )
          !$omp parallel default(SHARED) PRIVATE(x,y,Ex,Ey,ic_x,ic_y,off_x,off_y,temp1,temp2,temp,tmp5,tmp6,thread_id,q_accum_CS)
          !$&omp FIRSTPRIVATE(dt,xmin,xmax,ymin,ymax,ncx,rdx,rdy)
#ifdef _OPENMP
          thread_id = OMP_GET_THREAD_NUM()
#endif    
          call reset_charge_accumulator_to_zero_CS ( sim%q_accumulator_CS(thread_id+1)%q )
          q_accum_CS => sim%q_accumulator_CS(thread_id+1)%q
          !$omp do
          do i = 1, sim%parts_number
             SLL_INTERPOLATE_FIELD_CS(Ex,Ey,accumE_CS,p(i),temp1)
             GET_PARTICLE_POSITION(p(i),sim%m2d,x,y)
             x = modulo(x + dt*Ey, xmax - xmin)
             y = modulo(y - dt*Ex, ymax - ymin)
             SET_PARTICLE_POSITION(p(i),xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp5,tmp6)
             SLL_ACCUMULATE_PARTICLE_CHARGE_CS(q_accum_CS,p(i),temp2,temp)
          enddo
          !$omp end do   
          !$omp end parallel
          call sum_accumulators_CS( sim%q_accumulator_CS, n_threads, ncx*ncy )
          call sll_convert_charge_to_rho_2d_per_per_CS( sim%q_accumulator_CS(1)%q, sim%rho ) 

       else

          call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
          !$omp parallel default(SHARED) PRIVATE(x,y,Ex,Ey,ic_x,ic_y,off_x,off_y,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,thread_id,q_accum)
          !$&omp FIRSTPRIVATE(dt,xmin,xmax,ymin,ymax,ncx,rdx,rdy)
#ifdef _OPENMP
          thread_id = OMP_GET_THREAD_NUM()
#endif    
          call reset_charge_accumulator_to_zero ( sim%q_accumulator(thread_id+1)%q )
          q_accum => sim%q_accumulator(thread_id+1)%q
          !$omp do
          do i = 1, sim%parts_number
             SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,p(i),tmp3,tmp4)
             GET_PARTICLE_POSITION(p(i),sim%m2d,x,y)
             x = modulo(x + dt*Ey, xmax - xmin)
             y = modulo(y - dt*Ex, ymax - ymin)
             SET_PARTICLE_POSITION(p(i),xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp5,tmp6)
             SLL_ACCUMULATE_PARTICLE_CHARGE(q_accum,p(i),tmp7,tmp8)
          enddo
          !$omp end do   
          !$omp end parallel
          call sum_accumulators( sim%q_accumulator, n_threads, ncx*ncy )
          call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator(1)%q, sim%rho ) 
       endif
       ! ----- END PUSH PARTICLES -----
       ! -----------------------------

       if (sim%world_size > 1 ) then
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
       endif
       if  (mod(it+1, 100)==0 .and. sim%my_rank == 0) then
!          tfin = sll_time_elapsed_since(tinit)
!!$          write(it_name,'(i4.4)') it+1
!!$          print*, 'tfin at iter',it_name, 'is', tfin
!!$          nnnom = 'parts_at'//trim(adjustl(it_name))//'.dat'
!!$          open(50,file=nnnom)
!!$          do i = 1, sim%parts_number
!!$             GET_PARTICLE_POSITION( p(i),sim%m2d,x,y )
!!$             write(50,*) x, y
!!$          enddo
!!$          close(50)
          call sll_gnuplot_2d(xmin, sim%m2d%eta1_max, ncx+1, ymin, &
               sim%m2d%eta2_max, ncy+1, &
               sim%rho, 'rho_GC_eeuler_it', it+1, ierr )
       endif
       call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, -sim%rho )

    enddo
    ! ----------   END TIME LOOP   -----------
    ! ----------------------------------------
    
    SLL_DEALLOCATE(sim%rho,   ierr)
    SLL_DEALLOCATE(sim%E1,    ierr)
    SLL_DEALLOCATE(sim%E2,    ierr)
    SLL_DEALLOCATE(phi, ierr)
    if (sim%world_size > 1 ) then
       SLL_DEALLOCATE_ARRAY(rho1d_send, ierr)
       SLL_DEALLOCATE_ARRAY(rho1d_receive, ierr)
    endif
!!$    SLL_DEALLOCATE_ARRAY(diag_TOTenergy, ierr)
!!$    SLL_DEALLOCATE_ARRAY(diag_energy, ierr)

  end subroutine run_2d_pic_cartesian_ExpEuler


  subroutine delete_4d_pic_cart( sim )
    type(sll_pic_simulation_2d_gc_cartesian) :: sim
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

    x = modulo(x,mesh%eta1_max - mesh%eta1_min)
    y = modulo(y,mesh%eta2_max - mesh%eta2_min)
  end subroutine apply_periodic_bc

  subroutine normL2_field(res,nx,ny,e,dx,dy)
    sll_real64, intent(out) :: res
    sll_real64, intent(in) :: dx,dy
    sll_int32, intent(in) :: nx,ny
    sll_real64, dimension(1:nx+1,1:ny+1),intent(in) :: e
    sll_int32 :: i,j
    
    res = 0._f64
    do j=1,ny
       do i=1,nx
          res = res + e(i,j)**2
       enddo
    enddo
    res = sqrt(res*dx*dy)
  end subroutine normL2_field

  subroutine normL2_field_E(res,nx,ny,e1,e2,dx,dy)
    sll_real64, intent(out) :: res
    sll_real64, intent(in) :: dx,dy
    sll_int32, intent(in) :: nx,ny
    sll_real64, dimension(1:nx+1,1:ny+1),intent(in) :: e1
    sll_real64, dimension(1:nx+1,1:ny+1),intent(in) :: e2
    sll_int32 :: i,j
    
    res = 0._f64
    do j=1,ny
       do i=1,nx
          res = res + e1(i,j)**2 + e2(i,j)**2
       enddo
    enddo
    res = sqrt(res*dx*dy)
  end subroutine normL2_field_E

  subroutine lognormL2_field_Ex (res,nx,ny,e,dx,dy)
    sll_real64, intent(out) :: res
    sll_real64, intent(in) :: dx,dy
    sll_int32, intent(in) :: nx,ny
    sll_real64, dimension(1:nx+1,1:ny+1),intent(in) :: e
    sll_int32 :: i,j
    
    res = 0._f64
    do j=1,ny
       do i=1,nx
          res = res + e(i,j)*e(i,j)
       enddo
    enddo
    res=res*dx*dy
    res = log(res)*0.5_f64
  end subroutine lognormL2_field_Ex
  
end module sll_m_sim_pic_gc_2d0v_cart_optim
