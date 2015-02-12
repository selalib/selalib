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

program unit_test_qn_solver_3d_polar_parallel_x1_wrapper
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_qn_solver_3d_polar_parallel_x1_wrapper_module
use sll_collective

!use sll_boundary_condition_descriptors

implicit none
  
  class(sll_poisson_3d_base), pointer :: poisson 
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
  type(layout_2D), pointer :: layout2d_parx1
  type(layout_2D), pointer :: layout2d_parx2
  sll_int32 :: loc2d_sz_x1
  sll_int32 :: loc2d_sz_x2
  
  x1_min = 0.2_f64
  x1_max = 0.8_f64

  
  Nc_x1 = 32
  Nc_x2 = 64
  Nc_x3 = 8
  
  
  call sll_boot_collective()
   
  num_proc = sll_get_collective_size(sll_world_collective)  
  nproc_x1 = num_proc
  nproc_x2 = 1
  layout2d_parx1  => new_layout_2D( sll_world_collective )
  call initialize_layout_with_distributed_2D_array( &
    Nc_x1+1, & 
    Nc_x2, & 
    nproc_x1, &
    nproc_x2, &
    layout2d_parx1 )

  call compute_local_sizes_2d( &
    layout2d_parx1, &
    loc2d_sz_x1, &
    loc2d_sz_x2)



  SLL_ALLOCATE(phi(loc2d_sz_x1,loc2d_sz_x2,Nc_x3+1),ierr)
  SLL_ALLOCATE(rho(loc2d_sz_x1,loc2d_sz_x2,Nc_x3+1),ierr)

  nproc_x1 = 1
  nproc_x2 = num_proc
  layout2d_parx2  => new_layout_2D( sll_world_collective )
  call initialize_layout_with_distributed_2D_array( &
    Nc_x1+1, & 
    Nc_x2, & 
    nproc_x1, &
    nproc_x2, &
    layout2d_parx2 )



  
  rho = 1._f64
  
  err = 0._f64




  poisson => new_qn_solver_3d_polar_parallel_x1_wrapper( &
    layout2d_parx2, &
    layout2d_parx1, &
    x1_min, &
    x1_max, &
    Nc_x1, &
    Nc_x2, &
    Nc_x3, &
    SLL_NEUMANN_MODE_0, &
    SLL_DIRICHLET)
  
  call poisson%compute_phi_from_rho( phi, rho )
  print *,'#bounds on proc', sll_get_collective_rank(sll_world_collective),'=',maxval(phi),minval(phi)

  if(sll_get_collective_rank(sll_world_collective)==0)then
    print *,'#bounds on proc 0 =',maxval(phi),minval(phi)
    print *,'#bounds on proc 0.1 =',maxval(phi(:,:,1)),minval(phi(:,:,1))
    print *,'#bounds on proc 0.2 =',maxval(phi(:,:,Nc_x3)),minval(phi(:,:,Nc_x3))

    if(err==0)then    
      print *, '#PASSED'
    endif
  endif
  
  call sll_halt_collective()
  
end program unit_test_qn_solver_3d_polar_parallel_x1_wrapper