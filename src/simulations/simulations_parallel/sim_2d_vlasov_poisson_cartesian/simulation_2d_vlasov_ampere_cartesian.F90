!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
!
! Vlasov-Poisson simulation in 1Dx1D
! Landau and KEEN waves test cases
!
! contact: Pierre Navaro (navaro AT unistra.fr)
!
! current investigations:
!   High order splitting in time
!   KEEN waves with uniform and non uniform grid in velocity

module sll_simulation_2d_vlasov_ampere_cartesian

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"
#include "sll_utilities.h"
#include "sll_poisson_solvers.h"

use sll_collective
use sll_remapper
use sll_buffer_loader_utilities_module
use sll_constants
use sll_cartesian_meshes  
use sll_gnuplot_parallel
use sll_coordinate_transformation_2d_base_module
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations
use sll_common_array_initializers_module
use sll_parallel_array_initializer_module
use sll_module_advection_1d_periodic
use sll_module_advection_1d_non_uniform_cubic_splines
use sll_module_advection_1d_spectral
use sll_poisson_1d_periodic  
use sll_fft
use sll_simulation_base
use sll_time_splitting_coeff_module
use sll_module_poisson_1d_periodic_solver
use sll_module_poisson_1d_polar_solver
use sll_module_ampere_1d_pstd
use sll_simulation_2d_vlasov_poisson_cartesian

#ifdef _OPENMP
use omp_lib
#endif

implicit none

type, extends(sll_simulation_2d_vlasov_poisson_cart) :: &
     sll_simulation_2d_vlasov_ampere

 class(sll_ampere_1d_pstd),  pointer :: ampere

contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 procedure :: run 
         
end type sll_simulation_2d_vlasov_ampere

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function new_va2d_par_cart( filename ) result(sim)    

  type(sll_simulation_2d_vlasov_ampere), pointer :: sim    
  character(len=*), intent(in), optional               :: filename

  sll_int32  :: ierr
  sll_int32  :: num_cells_x1  
  sll_real64 :: x1_min
  sll_real64 :: x1_max
  sll_int32  :: tid

  SLL_ALLOCATE(sim, ierr)
  call init_vp2d_par_cart( sim, filename)

  num_cells_x1 = sim%mesh2d%num_cells1
  x1_min       = sim%mesh2d%eta1_min
  x1_max       = sim%mesh2d%eta1_max

  tid = 1

  sim%advect_x1(tid)%ptr => new_spectral_1d_advector( &
    num_cells_x1,                                     &
    x1_min,                                           &
    x1_max,                                           &
    SLL_PERIODIC)

  sim%ampere => new_ampere_1d_pstd( x1_min,       &
                                    x1_max,       &
                                    num_cells_x1)

     
end function new_va2d_par_cart

subroutine run(sim)

  class(sll_simulation_2d_vlasov_ampere), intent(inout) :: sim

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
  sll_int32                           :: np_x1_light
  sll_int32                           :: nproc_x1
  sll_int32                           :: nproc_x2
  sll_int32                           :: global_indices(2)
  sll_int32                           :: ierr
  sll_int32                           :: local_size_x1
  sll_int32                           :: local_size_x2
  type(poisson_1d_periodic)           :: poisson_1d
  sll_real64                          :: adr
  sll_real64                          :: tmp_loc(5)
  sll_real64                          :: tmp(5)
  sll_int32                           :: i
  sll_int32                           :: j
  sll_int32                           :: istep
  sll_int32                           :: ig
  sll_int32                           :: k
  
  sll_real64  ::   time, mass, momentum, l1norm, l2norm
  sll_real64  ::   kinetic_energy,potential_energy

  sll_real64, dimension(:), allocatable :: x2_array_unit
  sll_real64, dimension(:), allocatable :: x2_array_middle
  sll_real64, dimension(:), allocatable :: node_positions_x2
  sll_real64, dimension(:), allocatable :: node_positions_x2_light

  sll_int32                             :: file_id
  
  type(sll_fft_plan), pointer           :: pfwd
  sll_real64, dimension(:), allocatable :: buf_fft
  sll_comp64,dimension(:),allocatable   :: rho_mode
  
  sll_int32  :: nb_mode 
  sll_real64 :: t_step
  sll_real64 :: time_init
  sll_int32  :: split_istep
  sll_int32  :: num_dof_x2 
  sll_int32  :: num_dof_x2_light
   
  logical :: split_T
  
  ! for parallelization (output of distribution function in one single file)
  sll_int32                              :: collective_size
  sll_int32, dimension(:),   allocatable :: collective_displs
  sll_int32, dimension(:),   allocatable :: collective_recvcnts
  sll_real64,dimension(:,:), pointer     :: f_visu 
  sll_real64,dimension(:,:), pointer     :: f_visu_light 
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
  num_dof_x2_light = sim%num_dof_x2_light
  np_x1_light      = sim%light_size_x1
  
  if (MPI_MASTER) then

    print *,'#collective_size=',sll_get_collective_size(sll_world_collective)
    SLL_ALLOCATE(f_visu(np_x1,num_dof_x2),ierr)
    SLL_ALLOCATE(f_visu_buf1d(np_x1*num_dof_x2),ierr)
    SLL_ALLOCATE(f_visu_light(np_x1_light,num_dof_x2_light),ierr)

  else

    SLL_ALLOCATE(f_visu(1:1,1:1),ierr)          
    SLL_ALLOCATE(f_visu_buf1d(1:1),ierr)
    SLL_ALLOCATE(f_visu_light(1:1,1:1),ierr)          

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
  SLL_ALLOCATE(node_positions_x2_light(num_dof_x2),ierr)
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
      print *,'#sim%advection_form_x2=',sim%advection_form_x2
      SLL_ERROR('#not implemented')
  end select  
      

  if (MPI_MASTER) then

    call sll_binary_file_create("x.bdat", file_id, ierr)
    call sll_binary_write_array_1d(file_id,sim%x1_array(1:np_x1-1),ierr)
    call sll_binary_file_close(file_id,ierr)                    
    call sll_binary_file_create("v.bdat", file_id, ierr)
    call sll_binary_write_array_1d(file_id,node_positions_x2(1:np_x2-1),ierr)
    call sll_binary_file_close(file_id,ierr)                                             
  endif

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

  if (trim(sim%restart_file) /= "no_restart_file" ) then

    INQUIRE(FILE=trim(sim%restart_file)//'_proc_'//cproc//'.rst', EXIST=file_exists)

    if (.not. file_exists) then
      SLL_ERROR('#file '//trim(sim%restart_file)//'_proc_'//cproc//'.rst &
                & does not exist')
    endif

    open(unit=restart_id, &
         file=trim(sim%restart_file)//'_proc_'//cproc//'.rst', ACCESS="STREAM", &
         form='unformatted', IOStat=ierr)      

    if ( ierr .ne. 0 ) then
      SLL_ERROR('ERROR while opening file &
                 &'//trim(sim%restart_file)//'_proc_'//cproc//'.rst &
                 & Called from run_vp2d_cartesian().')
    end if

    print *,'#read restart file '//trim(sim%restart_file)//'_proc_'//cproc//'.rst'      
    call sll_binary_read_array_0d(restart_id,time_init,ierr)
    call sll_binary_read_array_2d(restart_id,f_x1,ierr)
    call sll_binary_file_close(restart_id,ierr)

  endif      
  
  if (sim%time_init_from_restart_file .eqv. .true.) then
    sim%time_init = time_init  
  endif
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

  if (MPI_MASTER) then

    call sll_binary_file_create('f0.bdat', file_id, ierr)
    call sll_binary_write_array_2d(file_id,f_visu(1:np_x1-1,1:np_x2-1),ierr)
    call sll_binary_file_close(file_id,ierr)

    call plot_f_cartesian( iplot,             &
                           f_visu,            &
                           sim%x1_array,      &
                           np_x1,             &
                           node_positions_x2, &
                           sim%num_dof_x2,    &
                           'f',               &
                           time_init )        

    print *,'#maxf',maxval(f_visu), minval(f_visu) 

  endif

  call initialize(poisson_1d,sim%mesh2d%eta1_min,sim%mesh2d%eta1_max,np_x1-1,ierr)

  !computation of electric field
  rho_loc = 0._f64
  ig = global_indices(2)-1
  do i=1,np_x1
    rho_loc(i) = rho_loc(i)&
      +sum(f_x1(i,:)*sim%integration_weight(1+ig:local_size_x2+ig))
  end do

  call sll_collective_allreduce( sll_world_collective, &
                                 rho_loc,              &
                                 np_x1,                &
                                 MPI_SUM,              &
                                 rho )

  rho = sim%factor_x2_1*1._f64-sim%factor_x2_rho*rho        
  
  call sim%poisson%compute_E_from_rho( efield, rho )
      
  ! Ponderomotive force at initial time. We use a sine wave
  ! with parameters k_dr and omega_dr.
  istep = 0
  e_app = 0._f64

  if (sim%driven) then

    call PFenvelope(adr,                    &
                    time_init+istep*sim%dt, &
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
               - sim%omegadr*(time_init+real(istep,f64)*sim%dt))
    enddo

  endif

  if (MPI_MASTER) then

    call sll_ascii_file_create('thdiag.dat', th_diag_id, ierr)

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

  endif
  
  !write also initial deltaf function
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
    print *,'#step=',0,time_init+real(0,f64)*sim%dt,'iplot=',iplot
  endif

  iplot = iplot+1  

  do istep = 1, sim%num_iterations

    if (mod(istep,sim%freq_diag)==0) then
      if (MPI_MASTER) then        
        print *,'#step=',istep,time_init+real(istep,f64)*sim%dt,'iplot=',iplot
      endif
    endif  

    split_T = sim%split%split_begin_T
    t_step = real(istep-1,f64)

    do split_istep=1,sim%split%nb_split_step

      if (split_T) then
        !! T ADVECTION 
        tid=1          
        !$OMP PARALLEL DEFAULT(SHARED) &
        !$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid) 
        !advection in x
#ifdef _OPENMP
          tid = omp_get_thread_num()+1
#endif
        !$OMP DO
        do i_omp = 1, local_size_x2

          ig_omp    = i_omp+global_indices(2)-1
          alpha_omp = sim%factor_x1*node_positions_x2(ig_omp) &
                    * sim%split%split_step(split_istep) 
          f1d_omp_in(1:np_x1,tid) = f_x1(1:np_x1,i_omp)
          
          call sim%advect_x1(tid)%ptr%advect_1d_constant(  &
            alpha_omp,                                     &
            sim%dt,                                        &
            f1d_omp_in(1:np_x1,tid),                       &
            f1d_omp_out(1:np_x1,tid))

          f_x1(1:np_x1,i_omp)=f1d_omp_out(1:np_x1,tid)

        end do
        !$OMP END DO          
        !$OMP END PARALLEL

        t_step = t_step+sim%split%split_step(split_istep)

          stop
        if (associated(sim%ampere)) then
          !Compute current
          rho_loc = 0._f64
          do i=1,np_x1
            do j= 1, local_size_x2
              global_indices = local_to_global( layout_x1, (/i, j/) )
              rho_loc(i)=rho_loc(i)                                &
              +f_x1(i,j)*sim%integration_weight(global_indices(2)) &
              *node_positions_x2(global_indices(2))
            end do
          end do
              
          call sll_collective_allreduce( sll_world_collective, &
                                         rho_loc,              &
                                         np_x1,                &
                                         MPI_SUM,              &
                                         current )

          call sim%ampere%compute_e_from_j( sim%dt, current, efield )

        else

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

          call sim%poisson%compute_E_from_rho( efield, rho )

        end if
        
        if (sim%driven) then

          call PFenvelope(adr,                     &
                          time_init+t_step*sim%dt, &
                          sim%tflat,               &
                          sim%tL,                  &
                          sim%tR,                  &
                          sim%twL,                 &
                          sim%twR,                 &
                          sim%t0,                  &
                          sim%turn_drive_off)

          do i = 1, np_x1
            e_app(i) = sim%Edrmax*adr*sim%kx      &
                     * sin(sim%kx*real(i-1,f64)   &
                     * sim%mesh2d%delta_eta1      &
                     - sim%omegadr*(time_init+t_step*sim%dt))
          enddo

        endif

      else

        !! V ADVECTION 
        !transposition
        call apply_remap_2D( remap_plan_x1_x2, f_x1, f_x2 )
        call compute_local_sizes( layout_x2, local_size_x1, local_size_x2 )
        global_indices(1:2) = local_to_global( layout_x2, (/1, 1/) )
        tid = 1

        !$OMP PARALLEL DEFAULT(SHARED) &
        !$OMP PRIVATE(i_omp,ig_omp,alpha_omp,tid,mean_omp,f1d) 
        !advection in v
#ifdef _OPENMP
          tid = omp_get_thread_num()+1
#endif
        !advection in v
        !$OMP DO
        do i_omp = 1,local_size_x1

          ig_omp=i_omp+global_indices(1)-1

          alpha_omp = -(efield(ig_omp)+e_app(ig_omp)) &
                    * sim%split%split_step(split_istep)

          f1d_omp_in(1:num_dof_x2,tid) = f_x2(i_omp,1:num_dof_x2)

          if (sim%advection_form_x2==SLL_CONSERVATIVE) then

            call function_to_primitive(f1d_omp_in(:,tid),    &
                                       x2_array_unit,        &
                                       np_x2-1,mean_omp)
          endif

          call sim%advect_x2(tid)%ptr%advect_1d_constant(    &
            alpha_omp,                                       &
            sim%dt,                                          &
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

        !transposition
        call apply_remap_2D( remap_plan_x2_x1, f_x2, f_x1 )
        call compute_local_sizes( layout_x1, local_size_x1, local_size_x2 )
        global_indices(1:2) = local_to_global( layout_x1, (/1, 1/) )

      endif
      !split_x= 1-split_x
      split_T = .not.(split_T)
    enddo
    
   !!DIAGNOSTICS
   if (mod(istep,sim%freq_diag_time)==0) then

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
      
      if (mod(istep,sim%freq_diag_restart)==0) then          
        call int2string(iplot,cplot) 
        call sll_binary_file_create('f_plot_'//cplot//'_proc_'//cproc//'.rst', &
                                    restart_id, ierr )
        call sll_binary_write_array_0d(restart_id,time,ierr)
        call sll_binary_write_array_2d(restart_id, &
                                       f_x1(1:local_size_x1,1:local_size_x2),ierr)
        call sll_binary_file_close(restart_id,ierr)    
      endif 

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
        
      if (mod(istep,sim%freq_diag)==0) then          
        !we substract f0
        !we gather in one file
        call load_buffer_2d( layout_x1, f_x1-f_x1_equil, f_x1_buf1d )

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

          call sll_gnuplot_write_1d(         &
            f_visu_buf1d(1:num_dof_x2),      &
            node_positions_x2(1:num_dof_x2), &
            'intdeltafdx',                   &
            iplot )

          call sll_gnuplot_write_1d(         &
            f_visu_buf1d(1:num_dof_x2),      &
            node_positions_x2(1:num_dof_x2), &
            'intdeltafdx')                        

          call sll_binary_write_array_2d(deltaf_id,           &
                                         f_visu(1:np_x1-1,1:np_x2-1),ierr)  

          call plot_f_cartesian( &
            iplot,               &
            f_visu,              &
            sim%x1_array,        &
            np_x1,               &
            node_positions_x2,   &
            sim%num_dof_x2,      &
            'deltaf',time)                    
        
        
        endif
        !we store f for visu
        call load_buffer_2d( layout_x1, f_x1, f_x1_buf1d )
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

          call sll_gnuplot_write_1d(         &
            f_visu_buf1d(1:num_dof_x2),      &
            node_positions_x2(1:num_dof_x2), &
            'intfdx',                        &
            iplot )

          call sll_gnuplot_write_1d(         &
            f_visu_buf1d(1:num_dof_x2),      &
            node_positions_x2(1:num_dof_x2), &
            'intfdx')

          call plot_f_cartesian( iplot,             &
                                 f_visu,            &
                                 sim%x1_array,      &
                                 np_x1,             &
                                 node_positions_x2, &
                                 sim%num_dof_x2,    &
                                 'f', time)                    
        endif

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

end subroutine run

end module sll_simulation_2d_vlasov_ampere_cartesian
