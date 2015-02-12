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
     type(sll_charge_accumulator_2d_CS), pointer  :: q_accumulator_CS
!!$     type(electric_field_accumulator), pointer :: E_accumulator
     type(electric_field_accumulator_CS), pointer :: E_accumulator_CS
     sll_real64, dimension(:,:), pointer :: rho
     type(poisson_2d_fft_solver), pointer :: poisson
     sll_real64, dimension(:,:), pointer :: E1, E2
   contains
     procedure, pass(sim) :: run => run_4d_pic_cartesian
     procedure, pass(sim) :: init_from_file => init_4d_pic_cartesian
!     procedure, pass(sim) :: run_manually => run_4d_pic_cartesian
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
    sll_real64  :: dt
    sll_int32   :: number_iterations
    sll_int32   :: NUM_PARTICLES, GUARD_SIZE, PARTICLE_ARRAY_SIZE
    sll_real64  :: THERM_SPEED
    sll_real64  :: QoverM, ALPHA
    sll_int32   :: NC_X,  NC_Y
    sll_real64  :: XMIN, KX, XMAX, YMIN, YMAX
    sll_int32, parameter  :: input_file = 99

    namelist /sim_params/ NUM_PARTICLES, GUARD_SIZE, &
                          PARTICLE_ARRAY_SIZE, &
                          THERM_SPEED, dt, number_iterations, &
                          QoverM, ALPHA
    namelist /grid_dims/  NC_X, NC_Y, XMIN, KX, YMIN, YMAX
    open(unit = input_file, file=trim(filename),IOStat=IO_stat)
    if( IO_stat /= 0 ) then
       print *, 'init_vp4d_par_cart() failed to open file ', filename
       STOP
    end if
    read(input_file, sim_params)
    read(input_file, grid_dims)
    close(input_file)

    XMAX = (2._f64*sll_pi/KX)
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
                                              sim%m2d%num_cells2  )
    call sll_initial_particles_4d( sim%thermal_speed_ions, & 
                                   ALPHA, KX, sim%m2d, &
                                   sim%ions_number,  &
                                   sim%part_group ) 

    call sll_sort_particles_2d( sim%sorter, sim%part_group )

!!$    sim%q_accumulator => new_charge_accumulator_2d( sim%m2d )    
    sim%q_accumulator_CS => new_charge_accumulator_2d_CS( sim%m2d )


!!$    sim%E_accumulator => new_field_accumulator_2d( sim%m2d )
    sim%E_accumulator_CS => new_field_accumulator_CS_2d( sim%m2d )
    
!!$    call sll_first_charge_accumulation_2d( sim%part_group, sim%q_accumulator )
    call sll_first_charge_accumulation_2d_CS( sim%part_group, sim%q_accumulator_CS )
    
  end subroutine init_4d_pic_cartesian

  ! Tentative function to RUN the simulation object 'manually'.

  subroutine run_4d_pic_cartesian( sim )

    class(sll_pic_simulation_4d_cartesian), intent(inout)  :: sim
    sll_int32  :: ierr, it, jj, counter
    sll_int32  :: i
    sll_real64 :: tmp1, tmp2, tmp3, tmp4, valeur
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
!!$    type(field_accumulator_cell), dimension(:), pointer :: accumE
    type(field_accumulator_CS), dimension(:), pointer :: accumE_CS
    type(sll_particle_2d_guard), dimension(:), pointer :: p_guard
    sll_real64, dimension(:,:), allocatable :: diag_energy! a memory buffer
    sll_real64, dimension(:,:), allocatable :: diag_AccMem! a memory buffer
    type(sll_time_mark)  :: t2, t3

    ncx = sim%m2d%num_cells1
    ncy = sim%m2d%num_cells2

    SLL_ALLOCATE(sim%rho(ncx+1,ncy+1),ierr)
    SLL_ALLOCATE( sim%E1(1:ncx+1,1:ncy+1), ierr )
    SLL_ALLOCATE( sim%E2(1:ncx+1,1:ncy+1), ierr )
    SLL_ALLOCATE(phi(1:ncx+1, 1:ncy+1), ierr)
    SLL_ALLOCATE(diag_energy(1:500, 1:2), ierr)
    SLL_ALLOCATE(diag_AccMem(0:sim%num_iterations-1, 1:2), ierr)

    p => sim%part_group%p_list
    qoverm = sim%part_group%qoverm
    p_guard => sim%part_group%p_guard
!!$    accumE => sim%E_accumulator%e_acc
    accumE_CS => sim%E_accumulator_CS%e_acc
    dt = sim%dt
    gi = 0
    xmin = sim%m2d%eta1_min
    ymin = sim%m2d%eta2_min
    rdx = 1._f64/sim%m2d%delta_eta1
    rdy = 1._f64/sim%m2d%delta_eta2

!!$    call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator, sim%rho ) 
    call sll_convert_charge_to_rho_2d_per_per_CS( sim%q_accumulator_CS, sim%rho ) 
    call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, sim%rho )

!!$    call reset_field_accumulator_to_zero( sim%E_accumulator )
    call reset_field_accumulator_CS_to_zero( sim%E_accumulator_CS )
!!$    call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
    call sll_accumulate_field_CS( sim%E1, sim%E2, sim%E_accumulator_CS )

!    open(58,file='verif_S1.dat')
!    open(68,file='verif_S3.dat')
    do i = 1, sim%ions_number
       pp_vx = p(i)%vx
       pp_vy = p(i)%vy
!       particle => p(i)
!!$       SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,p(i),tmp3,tmp4)
!!$       p(i)%vx = p(i)%vx - 0.5_f64 *dt* Ex* qoverm
!!$       p(i)%vy = p(i)%vy - 0.5_f64 *dt* Ey* qoverm
!!$       write(58,*) p(i)%vx, p(i)%vy
!
       SLL_INTERPOLATE_FIELD_CS(Ex_CS,Ey_CS,accumE_CS,p(i),ttmp)
       p(i)%vx = pp_vx - 0.5_f64 *dt* Ex_CS* qoverm
       p(i)%vy = pp_vy - 0.5_f64 *dt* Ey_CS* qoverm
!       write(68,*) p(i)%vx, p(i)%vy
    end do! half-step advection of the velocities by -dt/2 here
    
    open(65,file='logE_vals.dat')
    call sll_set_time_mark(t2)    

!  ----  TIME LOOP  ----
    do it = 0, sim%num_iterations-1

       call normL2_field_Ex ( valeur, sim%m2d%num_cells1, &
                             sim%m2d%num_cells2, sim%E1,  &
                             sim%m2d%delta_eta1, sim%m2d%delta_eta2 )
!!$       call norme_champs_x_etsin ( valeur1,  valeur2,  &
!!$            sim%m2d%num_cells1, sim%m2d%num_cells2, sim%E1,   &
!!$            sim%m2d%delta_eta1, sim%m2d%delta_eta2, it*sim%dt )
!!$       time = sll_time_elapsed_since(t0)
!!$       print*, 'time for norm calculation=', time

       counter = 1 + modulo(it,500)
       diag_energy(counter,:) = (/ it*sim%dt, valeur /)
!       diag_energy(counter,:) = (/ it*sim%dt, valeur1, valeur2, GAMMA*it*sim%dt /)

       if ( mod(it+1,500)==0 ) then
!!$          time = sll_time_elapsed_since(t2)
!!$          print*, 'iter=',it+1, 'TIME=', time, 100*sim%ions_number/time, 'average pushes/sec'
!!$          call sll_set_time_mark(t2)
          do jj=1,500
             write(65,*) diag_energy(jj,:)
          enddo
       endif

       if (mod(it+1,10)==0) then 
          print*, 'iter=', it+1
          call sll_sort_particles_2d( sim%sorter, sim%part_group )
       endif

       call reset_charge_accumulator_to_zero_CS( sim%q_accumulator_CS )

       call sll_set_time_mark(t3)
       ! ---- PUSH PARTICLES ----
       do i = 1, sim%ions_number,2
!!$          SLL_INTERPOLATE_FIELD(Ex,Ey,accumE,p(i),tmp3,tmp4)
!!$          SLL_INTERPOLATE_FIELD(Ex1,Ey1,accumE,p(i+1),tmp5,tmp6)
          SLL_INTERPOLATE_FIELD_CS(Ex,Ey,accumE_CS,p(i),ttmp1)
          SLL_INTERPOLATE_FIELD_CS(Ex1,Ey1,accumE_CS,p(i+1),ttmp2)
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
!!$             SLL_ACCUMULATE_PARTICLE_CHARGE(sim%q_accumulator,p(i),tmp1,tmp2)
             SLL_ACCUMULATE_PARTICLE_CHARGE_CS(sim%q_accumulator_CS,p(i),ttmp1,temp)
          else ! store reference for later processing
             gi = gi + 1
             p_guard(gi)%p => p(i)
          end if

          if(in_bounds( x1, y1, sim%m2d )) then ! finish push
             SET_PARTICLE_POSITION(p(i+1),xmin,ymin,ncx,x1,y1,ic_x1,ic_y1,off_x1,off_y1,rdx,rdy,tmp3,tmp4)
!!$             SLL_ACCUMULATE_PARTICLE_CHARGE(sim%q_accumulator,p(i+1),tmp3,tmp4)
             SLL_ACCUMULATE_PARTICLE_CHARGE_CS(sim%q_accumulator_CS,p(i+1),ttmp2,temp)
          else ! store reference for later processing
             gi = gi + 1
             p_guard(gi)%p => p(i+1)
          end if
       enddo
       ! ---- END PUSH PARTICLES ----
       ttime = sll_time_elapsed_since(t3)
       diag_AccMem(it,:) = (/ (it+1)*dt, (32*sim%ions_number*2 + gi*2*8 + &
            2*sizeof(sim%q_accumulator_CS%q_acc) + sizeof(sim%E_accumulator_CS%e_acc))/ttime/1e9 /)! access to memory in GB/sec
!  64*ncx*ncy + 2*32*ncx*ncy)/ttime/1e9 /)! access to memory in GB/sec

!!$       if (mod(it+1,10)==0) then 
!!$          print*, (32*sim%ions_number*2 + gi*2*8 + 64*ncx*ncy + 2*32*ncx*ncy)/ttime/1e6,'MeB/sec: access to memory'
!!$       endif

       ! Process the particles in the guard list. In the periodic case, no
       ! destruction of particles is needed, so this is simple.
       do i=1, gi
!          particle => p_guard(i)%p
          GET_PARTICLE_POSITION(p_guard(i)%p,sim%m2d,x,y)!(particle,sim%m2d,x,y)!
          x = x + dt * p_guard(i)%p%vx!particle%vx! 
          y = y + dt * p_guard(i)%p%vy!particle%vy! 
          call apply_periodic_bc( sim%m2d, x, y)
          SET_PARTICLE_POSITION(p_guard(i)%p,xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2)!particle
!!$          SLL_ACCUMULATE_PARTICLE_CHARGE(sim%q_accumulator,p_guard(i)%p,tmp1,tmp2)!particle
          SLL_ACCUMULATE_PARTICLE_CHARGE_CS(sim%q_accumulator_CS,p_guard(i)%p,ttmp,temp)!particle
       end do
       ! reset any counters
       gi = 1

!!$       call sll_convert_charge_to_rho_2d_per_per( sim%q_accumulator, sim%rho )
       call sll_convert_charge_to_rho_2d_per_per_CS( sim%q_accumulator_CS, sim%rho )
       call sim%poisson%compute_E_from_rho( sim%E1, sim%E2, sim%rho )

!!$       call reset_field_accumulator_to_zero( sim%E_accumulator )
       call reset_field_accumulator_CS_to_zero( sim%E_accumulator_CS )
!!$       call sll_accumulate_field( sim%E1, sim%E2, sim%E_accumulator )
       call sll_accumulate_field_CS( sim%E1, sim%E2, sim%E_accumulator_CS )

    enddo ! END TIME LOOP
    time = sll_time_elapsed_since(t2)
    print*, sim%num_iterations*sim%ions_number/time, 'average pushes/sec'
    close(65)

    print*, 'END --- write diags'

!    open(65,file='AccesstoMemory_withoutModuloforPerBC.dat')
    open(65,file='AccesstoMemory.dat')
    do jj = 0, sim%num_iterations-1
       if ( mod(jj+1,10) == 0 ) then
          write(65,*) diag_AccMem(jj,:), diag_AccMem(jj,2)
       else
          write(65,*) diag_AccMem(jj,:)
       endif
    enddo
    close(65)

    SLL_DEALLOCATE(sim%rho,   ierr)
    SLL_DEALLOCATE(sim%E1,    ierr)
    SLL_DEALLOCATE(sim%E2,    ierr)
    SLL_DEALLOCATE(phi, ierr)
    
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
