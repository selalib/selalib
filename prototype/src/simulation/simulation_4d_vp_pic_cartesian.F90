module sll_pic_simulation_4d_cartesian_module

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_utilities.h"
#include "sll_accumulators.h" 
#include "particle_representation.h"


  use sll_constants
  use sll_simulation_base
  use sll_logical_meshes
  use sll_timer
  use sll_particle_group_2d_module
  use sll_particle_initializers
  use sll_particle_sort_module
!  use sll_accumulators
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
     type(sll_particle_group_2d),  pointer :: part_group
     type(sll_logical_mesh_2d),    pointer :: m2d
     type(sll_particle_sorter_2d), pointer :: sorter
     type(sll_charge_accumulator_2d), pointer  :: q_accumulator
     type(electric_field_accumulator), pointer :: E_accumulator
     logical :: cubicsplines
     type(sll_charge_accumulator_2d_CS), pointer  :: q_accumulator_CS
     type(electric_field_accumulator_CS), pointer :: E_accumulator_CS
     sll_real64, dimension(:,:), pointer :: rho
     type(poisson_2d_fft_solver), pointer :: poisson
     sll_real64, dimension(:,:), pointer :: E1, E2
     sll_int32 :: my_rank
     sll_int32 :: world_size
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
    logical     :: CubicSplines
    sll_int32   :: NC_X,  NC_Y
    sll_real64  :: XMIN, KX, XMAX, YMIN, YMAX
    sll_int32, parameter  :: input_file = 99
    sll_int32, dimension(:), allocatable  :: rand_seed
    sll_int32   :: rand_seed_size

    namelist /sim_params/ NUM_PARTICLES, GUARD_SIZE, &
                          PARTICLE_ARRAY_SIZE, &
                          THERM_SPEED, dt, number_iterations, &
                          QoverM, ALPHA, CubicSplines
    namelist /grid_dims/  NC_X, NC_Y, XMIN, KX, YMIN, YMAX
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

    XMAX = (2._f64*sll_pi/KX)
    sim%cubicsplines = CubicSplines
    sim%thermal_speed_ions = THERM_SPEED
    sim%ions_number = NUM_PARTICLES
    sim%guard_size = GUARD_SIZE  
    sim%array_size = PARTICLE_ARRAY_SIZE
    sim%dt = dt
    sim%num_iterations = number_iterations

    sim%m2d =>  new_logical_mesh_2d( NC_X, NC_Y, &
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

    call sll_initial_particles_4d( sim%thermal_speed_ions, & 
                                   ALPHA, KX, sim%m2d,     &
                                   sim%ions_number,        &
                                   sim%part_group,         &
                                   rand_seed, sim%my_rank  ) 
    SLL_DEALLOCATE_ARRAY( rand_seed, ierr )

    call sll_sort_particles_2d( sim%sorter, sim%part_group )

    if (sim%cubicsplines) then
       sim%q_accumulator_CS => new_charge_accumulator_2d_CS( sim%m2d )       
       sim%E_accumulator_CS => new_field_accumulator_CS_2d( sim%m2d )           
       call sll_first_charge_accumulation_2d_CS( sim%part_group, sim%q_accumulator_CS )
    else
       print*, 'First order splines are used'  
       sim%q_accumulator => new_charge_accumulator_2d( sim%m2d )    
       sim%E_accumulator => new_field_accumulator_2d( sim%m2d )
       call sll_first_charge_accumulation_2d( sim%part_group, sim%q_accumulator )
    endif

    !$omp parallel
#ifdef _OPENMP
    print*, 'in the MAIN', OMP_GET_NUM_THREADS()
#endif
    !$omp end parallel
    
  end subroutine init_4d_pic_cartesian

  subroutine run_4d_pic_cartesian( sim )
!!!   calls of routines with '_CS' mean use of Cubic Splines  !!!
!     for deposition or interpolation step                      !
    class(sll_pic_simulation_4d_cartesian), intent(inout)  :: sim
    sll_int32  :: ierr, it, jj, counter
    sll_int32  :: i, j
    sll_real64 :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
    sll_real64 :: valeur
    sll_real64 :: ttmp(1:4,1:2), ttmp1(1:4,1:2), ttmp2(1:4,1:2)
    sll_real64, dimension(:,:), pointer :: phi
    sll_int32  :: ncx, ncy, ic_x,ic_y
    sll_int32  :: ic_x1,ic_y1
    sll_real32 :: off_x, off_y,off_x1,off_y1
    sll_real64 :: xmin, ymin, rdx, rdy
    sll_int32  :: gi ! counter index for guard list
    sll_real64 :: Ex, Ey, Ex1, Ey1,Ex_CS, Ey_CS
    sll_real64 :: qoverm
    sll_real64 :: x, x1  ! for global position
    sll_real64 :: y, y1  ! for global position
    sll_real64 :: dt, time, ttime, pp_vx, pp_vy, temp
    type(sll_particle_2d), dimension(:), pointer :: p
    type(field_accumulator_cell), dimension(:), pointer :: accumE
    type(field_accumulator_CS), dimension(:), pointer :: accumE_CS
    type(sll_particle_2d_guard), dimension(:), pointer :: p_guard
    sll_real64, dimension(:,:), allocatable :: diag_energy! a memory buffer
!!$    sll_real64, dimension(:,:), allocatable :: diag_AccMem! a memory buffer
    type(sll_time_mark)  :: t2, t3
    sll_real64, dimension(:), allocatable :: rho1d_send
    sll_real64, dimension(:), allocatable :: rho1d_receive
    sll_real64   :: t_init, t_fin

    ncx = sim%m2d%num_cells1
    ncy = sim%m2d%num_cells2

    SLL_ALLOCATE( rho1d_send(1:(ncx+1)*(ncy+1)),    ierr)
    SLL_ALLOCATE( rho1d_receive(1:(ncx+1)*(ncy+1)), ierr)

    SLL_ALLOCATE(sim%rho(ncx+1,ncy+1),ierr)
    SLL_ALLOCATE( sim%E1(1:ncx+1,1:ncy+1), ierr )
    SLL_ALLOCATE( sim%E2(1:ncx+1,1:ncy+1), ierr )
    SLL_ALLOCATE(phi(1:ncx+1, 1:ncy+1), ierr)
    SLL_ALLOCATE(diag_energy(1:500, 1:2), ierr)
!!$    SLL_ALLOCATE(diag_AccMem(0:sim%num_iterations-1, 1:2), ierr)

    p => sim%part_group%p_list
    qoverm = sim%part_group%qoverm
    p_guard => sim%part_group%p_guard

    dt = sim%dt
    gi = 0
    xmin = sim%m2d%eta1_min
    ymin = sim%m2d%eta2_min
    rdx = 1._f64/sim%m2d%delta_eta1
    rdy = 1._f64/sim%m2d%delta_eta2

    if (sim%cubicsplines) then
       accumE_CS => sim%E_accumulator_CS%e_acc
       call sll_convert_charge_to_rho_2d_per_per_CS( sim%q_accumulator_CS, sim%rho ) 
    else
       print*, 'First order splines are used'  
       accumE => sim%E_accumulator%e_acc
       call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator, sim%rho ) 
    endif

    it = 0
    call sll_gnuplot_corect_2d(xmin, sim%m2d%eta1_max, ncx+1, ymin, &
            sim%m2d%eta2_max, ncy+1, &
            sim%rho, 'rho_init', it, ierr )

    call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, sim%rho )

!!$    if (sim%my_rank == 0) then
!!$       call sll_gnuplot_corect_2d(xmin, sim%m2d%eta1_max, ncx+1, ymin, &
!!$            sim%m2d%eta2_max, ncy+1, &
!!$            sim%rho, 'rho_init', it, ierr )
!!$       call sll_gnuplot_corect_2d(xmin, sim%m2d%eta1_max, ncx+1, ymin, &
!!$            sim%m2d%eta2_max, ncy+1, &
!!$            sim%E1, 'E1_init', it, ierr )
!!$       call sll_gnuplot_corect_2d(xmin, sim%m2d%eta1_max, ncx+1, ymin, &
!!$            sim%m2d%eta2_max, ncy+1, &
!!$            sim%E2, 'E2_init', it, ierr )
!!$    endif

    if (sim%cubicsplines) then 
       call reset_field_accumulator_CS_to_zero( sim%E_accumulator_CS )
       call sll_accumulate_field_CS( sim%E1, sim%E2, sim%E_accumulator_CS )
    else
       call reset_field_accumulator_to_zero( sim%E_accumulator )
       call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
    endif

    !$omp parallel do PRIVATE (pp_vx, pp_vy, tmp3, tmp4, ttmp)
    !$&omp FIRSTPRIVATE(qoverm, dt, sim%cubicsplines)
    do i = 1, sim%ions_number
       pp_vx = p(i)%vx
       pp_vy = p(i)%vy
       if (sim%cubicsplines) then 
          SLL_INTERPOLATE_FIELD_CS(Ex_CS,Ey_CS,accumE_CS,p(i),ttmp)
          p(i)%vx = pp_vx - 0.5_f64 *dt* Ex_CS* qoverm
          p(i)%vy = pp_vy - 0.5_f64 *dt* Ey_CS* qoverm
       else
          SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,p(i),tmp3,tmp4)
          p(i)%vx = pp_vx - 0.5_f64 *dt* Ex* qoverm
          p(i)%vy = pp_vy - 0.5_f64 *dt* Ey* qoverm
       endif
    end do! half-step advection of the velocities by -dt/2 here
    !$omp end parallel do
    
    open(65,file='logE_vals.dat')
!!$    call sll_set_time_mark(t2)    
#ifdef _OPENMP
    t_init = omp_get_wtime()
#else
    call cpu_time(t_init)
#endif
!  ----  TIME LOOP  ----
    do it = 0, sim%num_iterations-1
       call normL2_field_Ex ( valeur, ncx, ncy, &
                              sim%E1,  &
                              sim%m2d%delta_eta1, sim%m2d%delta_eta2 )
       counter = 1 + modulo(it,500)
       diag_energy(counter,:) = (/ it*sim%dt, valeur /)

       if ( mod(it+1,500)==0 .and. (sim%my_rank == 0)) then
          do jj=1,500
             write(65,*) diag_energy(jj,:)
          enddo
       endif

       if (mod(it+1,10)==0) then 
           if (sim%my_rank == 0) print*, 'iter=', it+1
          call sll_sort_particles_2d( sim%sorter, sim%part_group )
       endif
       if (sim%cubicsplines) then 
          call reset_charge_accumulator_to_zero_CS( sim%q_accumulator_CS )
       else
          call reset_charge_accumulator_to_zero ( sim%q_accumulator)
       endif

!!$       call sll_set_time_mark(t3)
       ! ---- PUSH PARTICLES ----
       !$omp parallel do PRIVATE(x,y,x1,y1,Ex,Ey,Ex1,Ey1,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,ttmp1,ttmp2,off_x,off_y,ic_x,ic_y)
       !$&omp FIRSTPRIVATE(qoverm,dt,ncx,xmin,ymin,rdx,rdy,sim%cubicsplines)
! ! SOMETHING like this should be done for gi       !$omp do REDUCTION (+:gi)
       do i = 1, sim%ions_number,2
          if (sim%cubicsplines) then 
             SLL_INTERPOLATE_FIELD_CS(Ex,Ey,accumE_CS,p(i),ttmp1)
             SLL_INTERPOLATE_FIELD_CS(Ex1,Ey1,accumE_CS,p(i+1),ttmp2)
          else
             SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,p(i),tmp3,tmp4)
             SLL_INTERPOLATE_FIELD(Ex1,Ey1,accumE,p(i+1),tmp5,tmp6)
          endif
          p(i)%vx = p(i)%vx + dt * Ex* qoverm
          p(i)%vy = p(i)%vy + dt * Ey* qoverm
          p(i+1)%vx = p(i+1)%vx + dt * Ex1* qoverm
          p(i+1)%vy = p(i+1)%vy + dt * Ey1* qoverm
          GET_PARTICLE_POSITION(p(i),sim%m2d,x,y)
          GET_PARTICLE_POSITION(p(i+1),sim%m2d,x1,y1)
          x = x + dt * p(i)%vx
          y = y + dt * p(i)%vy
          x1 = x1 + dt * p(i+1)%vx
          y1 = y1 + dt * p(i+1)%vy
          if(in_bounds( x, y, sim%m2d )) then ! finish push
             SET_PARTICLE_POSITION(p(i),xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
             if (sim%cubicsplines) then 
                SLL_ACCUMULATE_PARTICLE_CHARGE_CS(sim%q_accumulator_CS,p(i),ttmp1,temp)
             else  
                SLL_ACCUMULATE_PARTICLE_CHARGE(sim%q_accumulator,p(i),tmp1,tmp2)
             endif
          else ! store reference for later processing
             gi = gi + 1
             p_guard(gi)%p => p(i)
          end if

          if(in_bounds( x1, y1, sim%m2d )) then ! finish push
             SET_PARTICLE_POSITION(p(i+1),xmin,ymin,ncx,x1,y1,ic_x1,ic_y1,off_x1,off_y1,rdx,rdy,tmp3,tmp4)
             if (sim%cubicsplines) then 
                SLL_ACCUMULATE_PARTICLE_CHARGE_CS(sim%q_accumulator_CS,p(i+1),ttmp2,temp)
             else
                SLL_ACCUMULATE_PARTICLE_CHARGE(sim%q_accumulator,p(i+1),tmp3,tmp4)                
             endif
          else ! store reference for later processing
             gi = gi + 1
             p_guard(gi)%p => p(i+1)
          end if
       enddo
       !$omp end parallel do
       ! ---- END PUSH PARTICLES ----

!!$       ttime = sll_time_elapsed_since(t3)
!!$       diag_AccMem(it,:) = (/ (it+1)*dt, (32*sim%ions_number*2 + gi*2*8 + &
!!$            2*sizeof(sim%q_accumulator_CS%q_acc) + sizeof(sim%E_accumulator_CS%e_acc))/ttime/1e9 /)! access to memory in GB/sec
!!$            2*sizeof(sim%q_accumulator%q_acc) + sizeof(sim%E_accumulator%e_acc))/ttime/1e9 /)! access to memory in GB/sec
!  64*ncx*ncy + 2*32*ncx*ncy)/ttime/1e9 /)! access to memory in GB/sec

!!$       if (mod(it+1,10)==0) then 
!!$          print*, (32*sim%ions_number*2 + gi*2*8 + 64*ncx*ncy + 2*32*ncx*ncy)/ttime/1e6,'MeB/sec: access to memory'
!!$       endif


       ! Process the particles in the guard list. In the periodic case, no
       ! destruction of particles is needed, so this is simple.

       !$omp parallel do PRIVATE(x,y,ic_x,ic_y,off_x,off_y,tmp1,tmp2,ttmp,temp)
       !$&omp FIRSTPRIVATE(dt,ncx,xmin,ymin,rdx,rdy,sim%cubicsplines)
       do i=1, gi
          GET_PARTICLE_POSITION(p_guard(i)%p,sim%m2d,x,y)
          x = x + dt * p_guard(i)%p%vx
          y = y + dt * p_guard(i)%p%vy
          call apply_periodic_bc( sim%m2d, x, y)
          SET_PARTICLE_POSITION(p_guard(i)%p,xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)
          if (sim%cubicsplines) then
             SLL_ACCUMULATE_PARTICLE_CHARGE_CS(sim%q_accumulator_CS,p_guard(i)%p,ttmp,temp)
          else
             SLL_ACCUMULATE_PARTICLE_CHARGE(sim%q_accumulator,p_guard(i)%p,tmp1,tmp2)
          endif
       end do
       !$omp end parallel do
!       ! reset any counters
       gi = 0!  --- ?????  gi = 1

       if (sim%cubicsplines) then
          call sll_convert_charge_to_rho_2d_per_per_CS( sim%q_accumulator_CS, sim%rho )
       else
          call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator, sim%rho )
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
     
!       if ( (sim%my_rank == 0).and.mod(it,3)==0) &
!            call sll_gnuplot_corect_2d(xmin, sim%m2d%eta1_max, ncx+1, ymin, sim%m2d%eta2_max, ncy+1, &
!            sim%rho, 'rhototal', it, ierr )
!
       call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, sim%rho )

       if (sim%cubicsplines) then
          call reset_field_accumulator_CS_to_zero( sim%E_accumulator_CS )
          call sll_accumulate_field_CS( sim%E1, sim%E2, sim%E_accumulator_CS )
       else
          call reset_field_accumulator_to_zero( sim%E_accumulator )
          call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
       endif

    enddo ! END TIME LOOP
#ifdef _OPENMP
    t_fin = omp_get_wtime()
#else
    call cpu_time(t_fin)
#endif
    close(65)
    print*, 'time is', t_fin-t_init, 'sec; and', int(sim%num_iterations,i64)*&
         int(sim%ions_number,i64)/(t_fin-t_init),'average pushes/sec for Proc',sim%my_rank

!!$    time = sll_time_elapsed_since(t2)
!!$    print*, int(sim%num_iterations,i64)*int(sim%ions_number,i64)/time, 'average pushes/sec for Proc', sim%my_rank

    if (sim%my_rank == 0) print*, 'END --- write diags'
!!$    if (sim%my_rank==0) open(65,file='AccesstoMemory_rk0.dat')
!!$! URGENT d'utiliser the rank_name !
!!$    if (sim%my_rank==1) open(66,file='AccesstoMemory_rk1.dat')
!!$    do jj = 0, sim%num_iterations-1
!!$       if ( mod(jj+1,10) == 0 ) then
!!$          if (sim%my_rank==0) write(65,*) diag_AccMem(jj,:), diag_AccMem(jj,2)
!!$          if (sim%my_rank==1) write(66,*) diag_AccMem(jj,:), diag_AccMem(jj,2)
!!$       else
!!$          if (sim%my_rank==0) write(65,*) diag_AccMem(jj,:)
!!$          if (sim%my_rank==1) write(66,*) diag_AccMem(jj,:)
!!$       endif
!!$    enddo
!!$    close(65) ; close(66)

    SLL_DEALLOCATE(sim%rho,   ierr)
    SLL_DEALLOCATE(sim%E1,    ierr)
    SLL_DEALLOCATE(sim%E2,    ierr)
    SLL_DEALLOCATE(phi, ierr)
    SLL_DEALLOCATE_ARRAY(rho1d_send, ierr)
    SLL_DEALLOCATE_ARRAY(rho1d_receive, ierr)
!!$    SLL_DEALLOCATE_ARRAY(diag_AccMem, ierr)
  end subroutine run_4d_pic_cartesian





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
    type(sll_logical_mesh_2d), pointer :: mesh

    res = (x >= mesh%eta1_min) .and. (x <= mesh%eta1_max) .and. &
          (y >= mesh%eta2_min) .and. (y <= mesh%eta2_max)
!!$    if((x >= mesh%eta1_min) .and. (x <= mesh%eta1_max) .and. &
!!$       (y >= mesh%eta2_min) .and. (y <= mesh%eta2_max)) then
!!$       res = .true.
!!$    else
!!$       res = .false.
!!$    end if
  end function in_bounds

  subroutine apply_periodic_bc( mesh, x, y )
    type(sll_logical_mesh_2d), pointer :: mesh
    sll_real64, intent(inout) :: x
    sll_real64, intent(inout) :: y
!!$    sll_real64 :: xmin
!!$    sll_real64 :: xmax
!!$    sll_real64 :: ymin
!!$    sll_real64 :: ymax

!!$    xmin = mesh%eta1_min
!!$    xmax = mesh%eta1_max
!!$    ymin = mesh%eta2_min
!!$    ymax = mesh%eta2_max
!!$    if( x < xmin ) x = x + xmax-xmin! create Branch with MOD instead of that
!!$    if( x > xmax ) x = x - xmax-xmin
!!$    if( y < ymin ) y = y + ymax-ymin
!!$    if( y > ymax ) y = y - ymax-ymin

    x = modulo(x,mesh%eta1_max - mesh%eta1_min)
    y = modulo(y,mesh%eta2_max - mesh%eta2_min)
    ! and the condition that the particle is in-bounds should trigger some
    ! alarm as this would not be supposed to happen here!
  end subroutine apply_periodic_bc

  subroutine normL2_field_Ex (res,nx,ny,e,dx,dy)
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
  end subroutine normL2_field_Ex
  
!!$  subroutine norme_champs_x_etsin( e_val, sin_val, nx, ny, e , dx, dy, temps )
!!$    sll_real64, intent(out) :: e_val, sin_val
!!$    sll_real64, intent(in) :: dx, dy, temps
!!$    sll_int32, intent(in)  :: nx, ny
!!$    sll_real64, dimension(1:nx+1,1:ny+1),intent(in) :: e
!!$    sll_int32  :: i,j
!!$    sll_real64 :: xxx
!!$    
!!$    e_val = 0._f64
!!$    sin_val = 0._f64
!!$!    print*, 'temps=', temps
!!$    do i=1,nx
!!$       xxx = XMIN + real(i-1,f64)*dx
!!$       sin_val = sin_val + ( sin(KX*xxx - OMEGA*temps) )**2
!!$       do j=1,ny
!!$          e_val   = e_val + e(i,j)*e(i,j)
!!$       enddo
!!$    enddo
!!$    e_val   = e_val*dx*dy
!!$    e_val = log(e_val)*0.5_f64
!!$
!!$    sin_val = sin_val*dx
!    print*, sin_val
!    print*, log(sin_val)*0.5_f64
!!$
!!$    sin_val = GAMMA*temps + log(sin_val)*0.5_f64
!!$    
!!$  end subroutine norme_champs_x_etsin
  
end module sll_pic_simulation_4d_cartesian_module
