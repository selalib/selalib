! Sample computation with the following characteristics:
! - vlasov-ampere
! - 1Dx1D cartesian: x1=x, x2=vx
! - parallel

program vlasov_ampere_2d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_utilities.h"
#include "sll_constants.h"

use sll_simulation_2d_vlasov_ampere_cartesian
use sll_common_array_initializers_module
use sll_parallel_array_initializer_module
use sll_collective
use sll_buffer_loader_utilities_module
use sll_remapper
use sll_fft
use sll_timer
use sll_module_advection_1d_ampere

!$ use omp_lib

implicit none

character(len=*), parameter :: this_prog_name = 'vlasov_ampere_2d'
character(len=128)          :: err_msg

class(sll_simulation_2d_vlasov_ampere_cart), pointer :: sim

character(len=256)  :: filename
character(len=256)  :: filename_local
type(sll_time_mark) :: t0
sll_real64          :: time
sll_int32           :: ierr
sll_int32           :: i

procedure(sll_scalar_initializer_2d), pointer :: init_func

logical                           :: init_from_unit_test  

sll_real64, dimension(:,:), pointer :: f_x1
sll_real64, dimension(:,:), pointer :: f_x2
sll_real64, dimension(:,:), pointer :: f_x1_equil

sll_real64, dimension(:),   pointer :: rho
sll_real64, dimension(:),   pointer :: efield
sll_real64, dimension(:),   pointer :: e_app
sll_real64, dimension(:),   pointer :: rho_loc
sll_real64, dimension(:),   pointer :: current

sll_int32 :: rhotot_id
sll_int32 :: efield_id     
sll_int32 :: adr_id
sll_int32 :: Edr_id
sll_int32 :: deltaf_id
sll_int32 :: t_id
sll_int32 :: th_diag_id
sll_int32 :: restart_id

type(layout_2D),            pointer :: layout_x1
type(layout_2D),            pointer :: layout_x2
type(remap_plan_2D_real64), pointer :: remap_plan_x1_x2
type(remap_plan_2D_real64), pointer :: remap_plan_x2_x1
sll_real64, dimension(:),   pointer :: f1d
sll_real64, dimension(:,:), pointer :: f1d_omp
sll_real64, dimension(:,:), pointer :: f1d_omp_in
sll_real64, dimension(:,:), pointer :: f1d_omp_out
sll_int32                           :: np_x1
sll_int32                           :: np_x2
sll_int32                           :: nproc_x1
sll_int32                           :: nproc_x2
sll_int32                           :: global_indices(2)
sll_int32                           :: local_size_x1
sll_int32                           :: local_size_x2
sll_real64                          :: adr
sll_real64                          :: tmp_loc(5)
sll_real64                          :: tmp(5)
sll_int32                           :: istep
sll_int32                           :: ig
sll_int32                           :: k

sll_real64  ::   mass, momentum, l1norm, l2norm
sll_real64  ::   kinetic_energy,potential_energy

sll_real64, dimension(:), allocatable :: x2_array_unit
sll_real64, dimension(:), allocatable :: x2_array_middle
sll_real64, dimension(:), allocatable :: node_positions_x2

sll_int32                             :: file_id

type(sll_fft_plan), pointer           :: pfwd
sll_real64, dimension(:), allocatable :: buf_fft
sll_comp64, dimension(:), allocatable :: rho_mode

sll_int32  :: nb_mode 
sll_int32  :: split_istep
sll_int32  :: num_dof_x2 
sll_real64 :: time_init
sll_real64 :: t_step
 
logical :: split_T

! for parallelization (output of distribution function in one single file)
sll_int32                              :: collective_size
sll_int32, dimension(:),   allocatable :: collective_displs
sll_int32, dimension(:),   allocatable :: collective_recvcnts
sll_real64,dimension(:,:), pointer     :: f_visu 
sll_real64,dimension(:),   pointer     :: f_visu_buf1d
sll_real64,dimension(:),   pointer     :: f_x1_buf1d
sll_real64,dimension(:),   pointer     :: f_hat_x2_loc
sll_real64,dimension(:),   pointer     :: f_hat_x2
sll_int32                              :: iplot
character(len=4)                       :: cproc
character(len=4)                       :: cplot
sll_int32                              :: iproc
logical                                :: file_exists
sll_int32                              :: tid
sll_int32                              :: i_omp
sll_int32                              :: ig_omp
sll_real64                             :: alpha_omp
sll_real64                             :: mean_omp

logical                                :: MPI_MASTER

init_from_unit_test = .false.

call sll_boot_collective()

if(sll_get_collective_rank(sll_world_collective)==0)then

  print *, '#Start time mark t0'
  call sll_set_time_mark(t0)
  print *, '#Booting parallel environment...'

endif

call get_command_argument(1, filename)
filename_local = trim(filename)
sim => new_va2d_par_cart( filename_local )
      
if (sll_get_collective_rank(sll_world_collective)==0) then
  MPI_MASTER = .true.
else
  MPI_MASTER = .false.
end if

iplot = 1

nb_mode          = sim%nb_mode
time_init        = sim%time_init
np_x1            = sim%mesh2d%num_cells1+1
np_x2            = sim%mesh2d%num_cells2+1
num_dof_x2       = sim%num_dof_x2

if (MPI_MASTER) then
  print *,'#collective_size=',sll_get_collective_size(sll_world_collective)
  SLL_ALLOCATE(f_visu(np_x1,num_dof_x2),ierr)
  SLL_ALLOCATE(f_visu_buf1d(np_x1*num_dof_x2),ierr)
else
  SLL_ALLOCATE(f_visu(1:1,1:1),ierr)          
  SLL_ALLOCATE(f_visu_buf1d(1:1),ierr)
endif

collective_size = sll_get_collective_size(sll_world_collective)
SLL_ALLOCATE(collective_displs(collective_size),ierr)
SLL_ALLOCATE(collective_recvcnts(collective_size),ierr)

SLL_ALLOCATE(buf_fft(np_x1-1),ierr)
pfwd => fft_new_plan(np_x1-1,buf_fft,buf_fft,FFT_FORWARD,FFT_NORMALIZE)
SLL_ALLOCATE(rho_mode(0:nb_mode),ierr)      

layout_x1       => new_layout_2D( sll_world_collective )
layout_x2       => new_layout_2D( sll_world_collective )    
nproc_x1 = sll_get_collective_size( sll_world_collective )
nproc_x2 = 1
call initialize_layout_with_distributed_array( &
  np_x1, num_dof_x2, nproc_x1, nproc_x2, layout_x2 )
call initialize_layout_with_distributed_array( &
  np_x1, num_dof_x2, nproc_x2, nproc_x1, layout_x1 )

call sll_view_lims( layout_x1 )
call sll_view_lims( layout_x2 )

call compute_local_sizes( layout_x2, local_size_x1, local_size_x2 )
SLL_ALLOCATE(f_x2(local_size_x1,local_size_x2),ierr)

call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices(1:2) = local_to_global( layout_x1, (/1, 1/) )
SLL_ALLOCATE(f_x1(local_size_x1,local_size_x2),ierr)    
SLL_ALLOCATE(f_x1_equil(local_size_x1,local_size_x2),ierr)    
SLL_ALLOCATE(f_x1_buf1d(local_size_x1*local_size_x2),ierr)    

remap_plan_x1_x2 => NEW_REMAP_PLAN(layout_x1, layout_x2, f_x1)
remap_plan_x2_x1 => NEW_REMAP_PLAN(layout_x2, layout_x1, f_x2)

SLL_ALLOCATE(rho(np_x1),ierr)
SLL_ALLOCATE(rho_loc(np_x1),ierr)
SLL_ALLOCATE(current(np_x1),ierr)
SLL_ALLOCATE(efield(np_x1),ierr)
SLL_ALLOCATE(e_app(np_x1),ierr)
SLL_ALLOCATE(f1d(max(np_x1,np_x2)),ierr)
SLL_ALLOCATE(f1d_omp(max(np_x1,np_x2),sim%num_threads),ierr)
SLL_ALLOCATE(f1d_omp_in(max(np_x1,np_x2),sim%num_threads),ierr)
SLL_ALLOCATE(f1d_omp_out(max(np_x1,np_x2),sim%num_threads),ierr)
SLL_ALLOCATE(x2_array_unit(np_x2),ierr)
SLL_ALLOCATE(x2_array_middle(np_x2),ierr)
SLL_ALLOCATE(node_positions_x2(num_dof_x2),ierr)
SLL_ALLOCATE(f_hat_x2_loc(nb_mode+1),ierr)
SLL_ALLOCATE(f_hat_x2(nb_mode+1),ierr)

x2_array_unit(1:np_x2) = &
  (sim%x2_array(1:np_x2)-sim%x2_array(1))/(sim%x2_array(np_x2)-sim%x2_array(1))
do i = 1, np_x2-1
   x2_array_middle(i) = 0.5_f64*(sim%x2_array(i)+sim%x2_array(i+1))
end do
x2_array_middle(np_x2) = x2_array_middle(1)+sim%x2_array(np_x2)-sim%x2_array(1)

select case (sim%advection_form_x2)
  case (SLL_ADVECTIVE)
    node_positions_x2(1:num_dof_x2) = sim%x2_array(1:num_dof_x2)
  case (SLL_CONSERVATIVE)
    node_positions_x2(1:num_dof_x2) = x2_array_middle(1:num_dof_x2)
  case default
    write(err_msg,*) '#sim%advection_form_x2=', sim%advection_form_x2, &
                     ' not implemented'
    SLL_ERROR( this_prog_name, err_msg )
end select  
    
call sll_2d_parallel_array_initializer_cartesian( &
   layout_x1,                                     &
   sim%x1_array,                                  &
   node_positions_x2,                             &
   f_x1,                                          &
   sim%init_func,                                 &
   sim%params)

iproc = sll_get_collective_rank(sll_world_collective)
call int2string(iproc, cproc)
call int2string(iplot, cplot)    

call check_restart()
time_init = sim%time_init

!ES initialise f_x1_equil that is used to define deltaf: 
!ES deltaf =  f_x1 - f_x1_equil
call sll_2d_parallel_array_initializer_cartesian( &
   layout_x1,                                     &
   sim%x1_array,                                  &
   node_positions_x2,                             &
   f_x1_equil,                                    &
   sim%equil_func,                                &
   sim%equil_func_params)

call compute_displacements_array_2d( layout_x1,       &
                                     collective_size, &
                                     collective_displs )

collective_recvcnts = receive_counts_array_2d( layout_x1, &
                                               collective_size )

!ES replace f_x1 by f_x1_equil and write f_x1_equil into 
!ES f0.bdat instead of f_x1  

call load_buffer_2d( layout_x1, f_x1_equil, f_x1_buf1d )

call sll_collective_gatherv_real64( sll_world_collective,        &
                                    f_x1_buf1d,                  &
                                    local_size_x1*local_size_x2, &
                                    collective_recvcnts,         &
                                    collective_displs,           &
                                    0,                           &
                                    f_visu_buf1d )
    
f_visu = reshape(f_visu_buf1d, shape(f_visu))

call compute_rho()

call sim%poisson%compute_E_from_rho( efield, rho )
    
! Ponderomotive force at initial time. We use a sine wave
! with parameters k_dr and omega_dr.
istep = 0

if (sim%driven) then
  call set_e_app(time_init)
else
  e_app = 0._f64
end if

if (MPI_MASTER) call write_init_files()

call load_buffer_2d( layout_x1, f_x1-f_x1_equil, f_x1_buf1d )

call sll_collective_gatherv_real64( sll_world_collective,        &
                                    f_x1_buf1d,                  &
                                    local_size_x1*local_size_x2, &
                                    collective_recvcnts,         &
                                    collective_displs,           &
                                    0,                           &
                                    f_visu_buf1d )

f_visu = reshape(f_visu_buf1d, shape(f_visu))

if (MPI_MASTER) then
  call sll_binary_write_array_2d(deltaf_id,f_visu(1:np_x1-1,1:np_x2-1),ierr)
  print *,'#step=',0,time_init,'iplot=',iplot
endif

iplot = iplot+1  

do istep = 1, sim%num_iterations

   split_T = sim%split%split_begin_T
   t_step = real(istep-1,f64)

   do split_istep=1,sim%split%nb_split_step

     if (split_T) then

       if (sim%ampere) then
         call advection_ampere_x(sim%split%split_step(split_istep)*sim%dt)
       else
         call advection_poisson_x(sim%split%split_step(split_istep)*sim%dt)
       end if
       t_step = t_step+sim%split%split_step(split_istep)

     else

       if (sim%driven) call set_e_app(time_init+(istep-1)*sim%dt)
       call transpose_xv()
       call advection_v(sim%split%split_step(split_istep)*sim%dt)
       call transpose_vx()

     endif
     split_T = .not.(split_T)
   enddo

  if (mod(istep,sim%freq_diag_time)==0) then

    call diagnostics()
    
    if (mod(istep,sim%freq_diag_restart)==0) then          
      call save_for_restart()
    endif 

    if (mod(istep,sim%freq_diag)==0) then          

      if (MPI_MASTER) then        
        print *,'#step=',istep,time_init+real(istep,f64)*sim%dt,'iplot=',iplot
      endif

      call gnuplot_write(f_x1-f_x1_equil, 'deltaf', 'intdeltafdx')
      call gnuplot_write(f_x1,                 'f',      'intfdx')

      iplot = iplot+1  
                
    endif

 end if

enddo

if (MPI_MASTER) then
  call sll_ascii_file_close(th_diag_id,ierr) 
  call sll_binary_file_close(deltaf_id,ierr) 
  call sll_binary_file_close(efield_id,ierr)
  call sll_binary_file_close(rhotot_id,ierr)
  call sll_binary_file_close(t_id,ierr)
  if (sim%driven) then
    call sll_binary_file_close(Edr_id,ierr)
    call sll_binary_file_close(adr_id,ierr)
  endif   
endif

if(sll_get_collective_rank(sll_world_collective)==0)then

  print *, '#reached end of va2d test'
  time = sll_time_elapsed_since(t0)
  print *, '#time elapsed since t0 : ',time
  print *, '#PASSED'

endif

call sll_halt_collective()

contains


subroutine advection_poisson_x(delta_t)
sll_real64 :: delta_t

call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x1, (/1, 1/) )

tid=1          
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid) 
!advection in x
!$ tid = omp_get_thread_num()+1
!$OMP DO
do i_omp = 1, local_size_x2

  ig_omp    = i_omp+global_indices(2)-1
  alpha_omp = sim%factor_x1*node_positions_x2(ig_omp)
  f1d_omp_in(1:np_x1,tid) = f_x1(1:np_x1,i_omp)
  
  call sim%advect_x1(tid)%ptr%advect_1d_constant(  &
    alpha_omp,                                     &
    delta_t,                                       &
    f1d_omp_in(1:np_x1,tid),                       &
    f1d_omp_out(1:np_x1,tid))

  f_x1(1:np_x1,i_omp)=f1d_omp_out(1:np_x1,tid)

end do
!$OMP END DO          
!$OMP END PARALLEL

call compute_rho()
call sim%poisson%compute_E_from_rho( efield, rho )

end subroutine advection_poisson_x

subroutine advection_ampere_x(delta_t)

sll_real64 :: delta_t
sll_int32  :: nc_x1
sll_comp64 :: s

call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x1, (/1, 1/) )

nc_x1 = np_x1-1

tid=1          


!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid) 
!advection in x
!$ tid = omp_get_thread_num()+1

sim%advect_ampere_x1(tid)%ptr%rk = cmplx(0.0,0.0,kind=f64)
!$OMP DO 
do i_omp = 1, local_size_x2

  ig_omp    = i_omp+global_indices(2)-1
  alpha_omp = sim%factor_x1*node_positions_x2(ig_omp)*delta_t
  f1d_omp_in(1:np_x1,tid) = f_x1(1:np_x1,i_omp)
  
  sim%advect_ampere_x1(tid)%ptr%d_dx = f1d_omp_in(1:nc_x1,tid)

  call fft_apply_plan(sim%advect_ampere_x1(tid)%ptr%fwx,  &
                      sim%advect_ampere_x1(tid)%ptr%d_dx, &
                      sim%advect_ampere_x1(tid)%ptr%fk)
  do i = 2, nc_x1/2+1
    sim%advect_ampere_x1(tid)%ptr%fk(i) = &
       sim%advect_ampere_x1(tid)%ptr%fk(i) & 
       * cmplx(cos(sim%advect_ampere_x1(tid)%ptr%kx(i)*alpha_omp), &
              -sin(sim%advect_ampere_x1(tid)%ptr%kx(i)*alpha_omp),kind=f64)
  end do

  sim%advect_ampere_x1(tid)%ptr%rk(2:nc_x1/2+1) = &
       sim%advect_ampere_x1(tid)%ptr%rk(2:nc_x1/2+1) &
     + sim%advect_ampere_x1(tid)%ptr%fk(2:nc_x1/2+1) * sim%integration_weight(ig_omp)

  call fft_apply_plan(sim%advect_ampere_x1(tid)%ptr%bwx, &
                      sim%advect_ampere_x1(tid)%ptr%fk,  &
                      sim%advect_ampere_x1(tid)%ptr%d_dx)

  f1d_omp_out(1:nc_x1, tid) = sim%advect_ampere_x1(tid)%ptr%d_dx/nc_x1
  f1d_omp_out(np_x1, tid)   = f1d_omp_out(1, tid) 

  f_x1(1:np_x1,i_omp)=f1d_omp_out(1:np_x1,tid)

end do
!$OMP END DO          

!$OMP END PARALLEL


sim%advect_ampere_x1(tid)%ptr%d_dx = efield(1:nc_x1)
call fft_apply_plan(sim%advect_ampere_x1(1)%ptr%fwx,  &
                    sim%advect_ampere_x1(1)%ptr%d_dx, &
                    sim%advect_ampere_x1(1)%ptr%ek)

do i = 2, nc_x1/2+1
  s = cmplx(0.0,0.0,kind=f64)
  do tid = 1, sim%num_threads
    s = s + sim%advect_ampere_x1(tid)%ptr%rk(i)
  end do
  sim%advect_ampere_x1(1)%ptr%rk(i) = s
end do

do i = 2, nc_x1/2+1
  sim%advect_ampere_x1(1)%ptr%ek(i) =  &
     - sim%advect_ampere_x1(1)%ptr%rk(i) * sim%L / (2*sll_pi*cmplx(0.,i-1,kind=f64))
end do

call fft_apply_plan(sim%advect_ampere_x1(1)%ptr%bwx, &
                    sim%advect_ampere_x1(1)%ptr%ek,  &
                    efield)

efield(1:nc_x1) = efield(1:nc_x1) / nc_x1
efield(np_x1) = efield(1)


!call compute_rho()
!call sim%poisson%compute_E_from_rho( efield, rho )

end subroutine advection_ampere_x

subroutine compute_rho()

call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x1, (/1, 1/) )
!computation of electric field
rho_loc = 0._f64
ig = global_indices(2)-1
do i=1,np_x1
  rho_loc(i)=rho_loc(i)                              &
    +sum(f_x1(i,1:local_size_x2)                     &
    *sim%integration_weight(1+ig:local_size_x2+ig))
end do
    
call sll_collective_allreduce( sll_world_collective, &
                               rho_loc,              &
                               np_x1,                &
                               MPI_SUM,              &
                               rho )

rho = sim%factor_x2_1*1._f64-sim%factor_x2_rho*rho

end subroutine compute_rho

subroutine compute_current()

sll_int32  :: global_indices(2)
sll_int32  :: local_size_x1
sll_int32  :: local_size_x2
sll_real64 :: v
sll_int32  :: j
sll_int32  :: gj

call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x1, (/1, 1/) )

do i = 1,local_size_x1
  rho_loc(i) = 0._f64
  do j = 1,local_size_x2
    global_indices = local_to_global( layout_x1, (/i, j/) )
    gj = global_indices(2)
    v  = node_positions_x2(gj)
    rho_loc(i) = rho_loc(i)+f_x1(i,j)*sim%integration_weight(gj)*v
  end do
end do

call sll_collective_allreduce( sll_world_collective, &
                               rho_loc,              &
                               np_x1,                &
                               MPI_SUM,              &
                               current )

end subroutine compute_current


subroutine advection_v(delta_t)
sll_real64 :: delta_t

call compute_local_sizes( layout_x2, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x2, (/1, 1/) )
tid = 1
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid,mean_omp,f1d) 
!advection in v
!$ tid = omp_get_thread_num()+1
!advection in v
!$OMP DO
do i_omp = 1,local_size_x1

  ig_omp=i_omp+global_indices(1)-1

  alpha_omp = -(efield(ig_omp)+e_app(ig_omp))

  f1d_omp_in(1:num_dof_x2,tid) = f_x2(i_omp,1:num_dof_x2)

  if (sim%advection_form_x2==SLL_CONSERVATIVE) then

    call function_to_primitive(f1d_omp_in(:,tid),    &
                               x2_array_unit,        &
                               np_x2-1,mean_omp)
  endif

  call sim%advect_x2(tid)%ptr%advect_1d_constant(    &
    alpha_omp,                                       &
    delta_t,                                         &
    f1d_omp_in(1:num_dof_x2,tid),                    &
    f1d_omp_out(1:num_dof_x2,tid))

  if (sim%advection_form_x2==SLL_CONSERVATIVE) then

    call primitive_to_function(f1d_omp_out(:,tid),   &
                               x2_array_unit,        &
                               np_x2-1,              &
                               mean_omp)
  endif
  f_x2(i_omp,1:num_dof_x2) = f1d_omp_out(1:num_dof_x2,tid)
end do
!$OMP END DO          
!$OMP END PARALLEL

end subroutine advection_v

subroutine transpose_xv()

call apply_remap_2D( remap_plan_x1_x2, f_x1, f_x2 )

end subroutine transpose_xv

subroutine transpose_vx()

call apply_remap_2D( remap_plan_x2_x1, f_x2, f_x1 )

end subroutine transpose_vx

subroutine gnuplot_write( f, fname, intfname)

sll_real64, dimension(:,:) :: f
character(len=*)           :: fname
character(len=*)           :: intfname

call load_buffer_2d( layout_x1, f, f_x1_buf1d )

call sll_collective_gatherv_real64( &
  sll_world_collective,             &
  f_x1_buf1d,                       &
  local_size_x1*local_size_x2,      &
  collective_recvcnts,              &
  collective_displs,                &
  0,                                &
  f_visu_buf1d )

f_visu = reshape(f_visu_buf1d, shape(f_visu))

if (MPI_MASTER) then

  do i=1,num_dof_x2
    f_visu_buf1d(i) = sum(f_visu(1:np_x1-1,i))*sim%mesh2d%delta_eta1
  enddo

  call sll_gnuplot_1d(         &
    f_visu_buf1d(1:num_dof_x2),      &
    node_positions_x2(1:num_dof_x2), &
    intfname,                        &
    iplot )

  call sll_gnuplot_1d(         &
    f_visu_buf1d(1:num_dof_x2),      &
    node_positions_x2(1:num_dof_x2), &
    intfname )                        

  call sll_binary_write_array_2d(deltaf_id,                      &
                                 f_visu(1:np_x1-1,1:np_x2-1),    &
                                 ierr)  

  call sll_plot_f_cartesian( &
        iplot,               &
        f_visu,              &
        sim%x1_array,        &
        np_x1,               &
        node_positions_x2,   &
        sim%num_dof_x2,      &
        fname,time)                    


endif

end subroutine gnuplot_write

subroutine diagnostics()

call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
global_indices = local_to_global( layout_x1, (/1, 1/) )

time             = time_init+real(istep,f64)*sim%dt
mass             = 0._f64
momentum         = 0._f64
l1norm           = 0._f64
l2norm           = 0._f64
kinetic_energy   = 0._f64
potential_energy = 0._f64
tmp_loc          = 0._f64
ig               = global_indices(2)-1               


do i = 1, np_x1-1        
  tmp_loc(1)= tmp_loc(1)+sum(f_x1(i,1:local_size_x2) &
    *sim%integration_weight(1+ig:local_size_x2+ig))
  tmp_loc(2)= tmp_loc(2)+sum(abs(f_x1(i,1:local_size_x2)) &
    *sim%integration_weight(1+ig:local_size_x2+ig))
  tmp_loc(3)= tmp_loc(3)+sum((f_x1(i,1:local_size_x2))**2 &
    *sim%integration_weight(1+ig:local_size_x2+ig))
  tmp_loc(4)= tmp_loc(4) +sum(f_x1(i,1:local_size_x2) &
    *sim%x2_array(global_indices(2)-1+1:global_indices(2)-1+local_size_x2) &
    *sim%integration_weight(1+ig:local_size_x2+ig))          
  tmp_loc(5)= tmp_loc(5)+sum(f_x1(i,1:local_size_x2) &
    *sim%x2_array(global_indices(2)-1+1:global_indices(2)-1+local_size_x2)**2 &
    *sim%integration_weight(1+ig:local_size_x2+ig) )          
end do

call sll_collective_allreduce( sll_world_collective, &
                               tmp_loc,              &
                               5,                    &
                               MPI_SUM,              &
                               tmp )

mass           = tmp(1)  * sim%mesh2d%delta_eta1 !* sim%mesh2d%delta_eta2
l1norm         = tmp(2)  * sim%mesh2d%delta_eta1 !* sim%mesh2d%delta_eta2
l2norm         = tmp(3)  * sim%mesh2d%delta_eta1 !* sim%mesh2d%delta_eta2
momentum       = tmp(4)  * sim%mesh2d%delta_eta1 !* sim%mesh2d%delta_eta2
kinetic_energy = 0.5_f64 * tmp(5) * sim%mesh2d%delta_eta1 !*sim%mesh2d%delta_eta2
potential_energy = 0._f64
do i=1, np_x1-1
  potential_energy = potential_energy+(efield(i)+e_app(i))**2
enddo
potential_energy = 0.5_f64*potential_energy* sim%mesh2d%delta_eta1

f_hat_x2_loc(1:nb_mode+1) = 0._f64
do i=1,local_size_x2
  buf_fft = f_x1(1:np_x1-1,i)
  call fft_apply_plan(pfwd,buf_fft,buf_fft)
  do k=0,nb_mode
    f_hat_x2_loc(k+1) = f_hat_x2_loc(k+1) &
      +abs(fft_get_mode(pfwd,buf_fft,k))**2 &
      *sim%integration_weight(ig+i)
  enddo
enddo

call sll_collective_allreduce( sll_world_collective, &
                               f_hat_x2_loc,         &
                               nb_mode+1,            &
                               MPI_SUM,              &
                               f_hat_x2 )

if (sll_get_collective_rank(sll_world_collective)==0) then                  

  buf_fft = rho(1:np_x1-1)
  call fft_apply_plan(pfwd,buf_fft,buf_fft)

  do k=0,nb_mode
    rho_mode(k)=fft_get_mode(pfwd,buf_fft,k)
  enddo  

  write(th_diag_id,'(f12.5,7g20.12)',advance='no') &
    time,                                          &
    mass,                                          &
    l1norm,                                        &
    momentum,                                      &
    l2norm,                                        &
    kinetic_energy,                                &
    potential_energy,                              &
    kinetic_energy + potential_energy

  do k=0,nb_mode
    write(th_diag_id,'(g20.12)',advance='no') abs(rho_mode(k))
  enddo

  do k=0,nb_mode-1
    write(th_diag_id,'(g20.12)',advance='no') f_hat_x2(k+1)
  enddo

  write(th_diag_id,'(g20.12)') f_hat_x2(nb_mode+1)

  call sll_binary_write_array_1d(efield_id,efield(1:np_x1-1),ierr)
  call sll_binary_write_array_1d(rhotot_id,rho(1:np_x1-1),ierr)
  call sll_binary_write_array_0d(t_id,time,ierr)
  if (sim%driven) then
    call sll_binary_write_array_1d(Edr_id,e_app(1:np_x1-1),ierr)
    call sll_binary_write_array_0d(adr_id,adr,ierr)
  endif   
endif

end subroutine diagnostics

subroutine write_init_files()

call sll_binary_file_create('f0.bdat', file_id, ierr)
call sll_binary_write_array_2d(file_id,f_visu(1:np_x1-1,1:np_x2-1),ierr)
call sll_binary_file_close(file_id,ierr)

call sll_plot_f_cartesian( iplot,             &
                           f_visu,            &
                           sim%x1_array,      &
                           np_x1,             &
                           node_positions_x2, &
                           sim%num_dof_x2,    &
                           'f',               &
                           time_init )        

print *,'#maxf',maxval(f_visu), minval(f_visu) 

call sll_binary_file_create("x.bdat", file_id, ierr)
call sll_binary_write_array_1d(file_id,sim%x1_array(1:np_x1-1),ierr)
call sll_binary_file_close(file_id,ierr)                    
call sll_binary_file_create("v.bdat", file_id, ierr)
call sll_binary_write_array_1d(file_id,node_positions_x2(1:np_x2-1),ierr)
call sll_binary_file_close(file_id,ierr)                                             

call sll_ascii_file_create(sim%thdiag_filename, th_diag_id, ierr)

call sll_binary_file_create('deltaf.bdat', deltaf_id, ierr)
call sll_binary_file_create('rhotot.bdat', rhotot_id, ierr)
call sll_binary_file_create('efield.bdat', efield_id, ierr)
call sll_binary_file_create('t.bdat', t_id, ierr)
call sll_binary_write_array_1d(efield_id,efield(1:np_x1-1),ierr)
call sll_binary_write_array_1d(rhotot_id,rho(1:np_x1-1),ierr)
call sll_binary_write_array_0d(t_id,real(istep,f64)*sim%dt,ierr)

if (sim%driven) then
  call sll_binary_file_create('adr.bdat', adr_id, ierr)
  call sll_binary_file_create('Edr.bdat', Edr_id, ierr)
  call sll_binary_write_array_1d(Edr_id,e_app(1:np_x1-1),ierr)
  call sll_binary_write_array_0d(adr_id,adr,ierr)
endif                    

end subroutine write_init_files

subroutine set_e_app(t)
sll_real64 :: t

  call PFenvelope(adr,                    &
                  t,                      &
                  sim%tflat,              &
                  sim%tL,                 &
                  sim%tR,                 &
                  sim%twL,                &
                  sim%twR,                &
                  sim%t0,                 &
                  sim%turn_drive_off)

  do i = 1, np_x1
    e_app(i) = sim%Edrmax*adr*sim%kx                          &
             * sin(sim%kx*real(i-1,f64)*sim%mesh2d%delta_eta1 &
             - sim%omegadr*t)
  enddo

end subroutine set_e_app

subroutine save_for_restart()

  call int2string(iplot,cplot) 
  call sll_binary_file_create('f_plot_'//cplot//'_proc_'//cproc//'.rst', &
                              restart_id, ierr )
  call sll_binary_write_array_0d(restart_id,time,ierr)
  call sll_binary_write_array_2d(restart_id, &
                                 f_x1(1:local_size_x1,1:local_size_x2),ierr)
  call sll_binary_file_close(restart_id,ierr)    

end subroutine save_for_restart

subroutine check_restart()
  character(len=*), parameter :: this_sub_name = 'check_restart'

if (trim(sim%restart_file) /= "no_restart_file" ) then


  INQUIRE(FILE=trim(sim%restart_file)//'_proc_'//cproc//'.rst', EXIST=file_exists)

  if (.not. file_exists) then
    err_msg = '#file '//trim(sim%restart_file)//'_proc_'//cproc//'.rst &
              & does not exist'
    SLL_ERROR( this_sub_name, err_msg )
  endif

  open(unit=restart_id, &
       file=trim(sim%restart_file)//'_proc_'//cproc//'.rst', ACCESS="STREAM", &
       form='unformatted', IOStat=ierr)      

  if ( ierr .ne. 0 ) then
    err_msg = 'ERROR while opening file &
               &'//trim(sim%restart_file)//'_proc_'//cproc//'.rst &
               & Called from run_va2d_cartesian().'
    SLL_ERROR( this_sub_name, err_msg )
  end if

  print *,'#read restart file '//trim(sim%restart_file)//'_proc_'//cproc//'.rst'      
  call sll_binary_read_array_0d(restart_id,time_init,ierr)
  call sll_binary_read_array_2d(restart_id,f_x1,ierr)
  call sll_binary_file_close(restart_id,ierr)

endif      

if (sim%time_init_from_restart_file .eqv. .true.) then
  sim%time_init = time_init  
endif

end subroutine check_restart

end program vlasov_ampere_2d
