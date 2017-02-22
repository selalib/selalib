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

program test_qn_solver_3d_polar_parallel_x1_wrapper

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_memory.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann_mode_0

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_poisson_3d_base, only: &
    sll_c_poisson_3d_base

  use sll_m_qn_solver_3d_polar_parallel_x1_wrapper, only: &
    sll_f_new_qn_solver_3d_polar_parallel_x1_wrapper

  use sll_m_remapper, only: &
    sll_o_compute_local_sizes, &
    sll_o_initialize_layout_with_distributed_array, &
    sll_f_new_layout_2d, &
    sll_t_layout_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  class(sll_c_poisson_3d_base), pointer :: poisson 
  sll_real64 :: err
  sll_real64 :: x1_min
  sll_real64 :: x1_max
  sll_int32 :: Nc_x1
  sll_int32 :: Nc_x2
  sll_int32 :: Nc_x3
  sll_real64,dimension(:,:,:),allocatable :: phi
  sll_real64,dimension(:,:,:),allocatable :: rho
  sll_int32 :: ierr
  sll_int32 :: num_proc
  sll_int32 :: nproc_x1
  sll_int32 :: nproc_x2
  type(sll_t_layout_2d), pointer :: layout2d_parx1
  type(sll_t_layout_2d), pointer :: layout2d_parx2
  sll_int32 :: loc2d_sz_x1
  sll_int32 :: loc2d_sz_x2
  
  x1_min = 0.2_f64
  x1_max = 0.8_f64

  
  Nc_x1 = 32
  Nc_x2 = 64
  Nc_x3 = 8
  
  
  call sll_s_boot_collective()
   
  num_proc = sll_f_get_collective_size(sll_v_world_collective)  
  nproc_x1 = num_proc
  nproc_x2 = 1
  layout2d_parx1  => sll_f_new_layout_2d( sll_v_world_collective )
  call sll_o_initialize_layout_with_distributed_array( &
    Nc_x1+1, & 
    Nc_x2, & 
    nproc_x1, &
    nproc_x2, &
    layout2d_parx1 )

  call sll_o_compute_local_sizes( &
    layout2d_parx1, &
    loc2d_sz_x1, &
    loc2d_sz_x2)



  SLL_ALLOCATE(phi(loc2d_sz_x1,loc2d_sz_x2,Nc_x3+1),ierr)
  SLL_ALLOCATE(rho(loc2d_sz_x1,loc2d_sz_x2,Nc_x3+1),ierr)

  nproc_x1 = 1
  nproc_x2 = num_proc
  layout2d_parx2  => sll_f_new_layout_2d( sll_v_world_collective )
  call sll_o_initialize_layout_with_distributed_array( &
    Nc_x1+1, & 
    Nc_x2, & 
    nproc_x1, &
    nproc_x2, &
    layout2d_parx2 )



  
  rho = 1._f64
  
  err = 0._f64




  poisson => sll_f_new_qn_solver_3d_polar_parallel_x1_wrapper( &
    layout2d_parx2, &
    layout2d_parx1, &
    x1_min, &
    x1_max, &
    Nc_x1, &
    Nc_x2, &
    Nc_x3, &
    sll_p_neumann_mode_0, &
    sll_p_dirichlet)
  
  call poisson%compute_phi_from_rho( phi, rho )
  print *,'#bounds on proc', sll_f_get_collective_rank(sll_v_world_collective),'=',maxval(phi),minval(phi)

  if(sll_f_get_collective_rank(sll_v_world_collective)==0)then
    print *,'#bounds on proc 0 =',maxval(phi),minval(phi)
    print *,'#bounds on proc 0.1 =',maxval(phi(:,:,1)),minval(phi(:,:,1))
    print *,'#bounds on proc 0.2 =',maxval(phi(:,:,Nc_x3)),minval(phi(:,:,Nc_x3))

    if(err==0)then    
      print *, '#PASSED'
    endif
  endif
  
  call sll_s_halt_collective()
  
end program test_qn_solver_3d_polar_parallel_x1_wrapper
