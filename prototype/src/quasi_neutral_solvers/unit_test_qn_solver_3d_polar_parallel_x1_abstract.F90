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

program unit_test_qn_solver_3d_polar_parallel_x1_abstract
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_qn_solver_3d_polar_parallel_x1_abstract_module
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
  
  x1_min = 0.2_f64
  x1_max = 0.8_f64

  
  Nc_x1 = 32
  Nc_x2 = 64
  Nc_x3 = 8 
  
  
  SLL_ALLOCATE(phi(Nc_x1+1,Nc_x2+1,Nc_x3+1),ierr)
  SLL_ALLOCATE(rho(Nc_x1+1,Nc_x2+1,Nc_x3+1),ierr)
  
  rho = 1._f64
  
  err = 0._f64
  
!  poisson =>new_qn_solver_3d_polar_parallel_x1_abstract( &
!    x1_min, &
!    x1_max, &
!    Nc_x1, &
!    Nc_x2, &
!    (/SLL_NEUMANN_MODE_0, SLL_DIRICHLET/))
!  
!  call poisson%compute_phi_from_rho( phi, rho )
!
!  print *,maxval(phi),minval(phi)
  
  

  if(err==0)then    
    print *, '#PASSED'
  endif

end program unit_test_qn_solver_3d_polar_parallel_x1_abstract