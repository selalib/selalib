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

program unit_test_advection_2d_CSL
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_module_advection_2d_CSL
use sll_module_characteristics_2d_explicit_euler_conservative

implicit none
  
  class(sll_advection_2d_base), pointer :: adv
  class(sll_characteristics_2d_base), pointer :: charac
  sll_real64 :: x1_min
  sll_real64 :: x1_max
  sll_real64 :: x2_min
  sll_real64 :: x2_max
  sll_int32 :: num_cells_x1
  sll_int32 :: num_cells_x2
  sll_real64, dimension(:,:), allocatable :: input
  sll_real64, dimension(:,:), allocatable :: output
  sll_real64, dimension(:), pointer :: x1_mesh
  sll_real64, dimension(:), pointer :: x2_mesh
  sll_real64 :: dt
  sll_real64,dimension(:,:), allocatable :: A1
  sll_real64,dimension(:,:), allocatable :: A2
  sll_real64 :: err
  sll_int32 :: ierr
  sll_int32 :: i
  sll_real64 :: delta_x1
  sll_real64 :: delta_x2
  
  x1_min = 0._f64
  x1_max = 1._f64
  x2_min = 0._f64
  x2_max = 1._f64
  num_cells_x1 = 32
  num_cells_x2 = 32
  dt = 0.1_f64 !0.1_f64
  
  delta_x1 = (x1_max-x1_min)/real(num_cells_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(num_cells_x2,f64)
  SLL_ALLOCATE(x1_mesh(num_cells_x1+1),ierr)
  SLL_ALLOCATE(x2_mesh(num_cells_x2+1),ierr)
  SLL_ALLOCATE(input(num_cells_x1+1,num_cells_x2+1),ierr)
  SLL_ALLOCATE(output(num_cells_x1+1,num_cells_x2+1),ierr)
  SLL_ALLOCATE(A1(num_cells_x1+1,num_cells_x2+1),ierr)
  SLL_ALLOCATE(A2(num_cells_x1+1,num_cells_x2+1),ierr)

  do i=1,num_cells_x1+1
    x1_mesh(i) = x1_min+real(i-1,f64)*delta_x1
  enddo

  do i=1,num_cells_x2+1
    x2_mesh(i) = x2_min+real(i-1,f64)*delta_x2
  enddo

  
  input = 1._f64

  A1 = 1._f64
  A2 = 1._f64

  err=0._f64




  charac => new_explicit_euler_conservative_2d_charac(&
      num_cells_x1+1, &
      num_cells_x2+1, &
      SLL_PERIODIC, &
      SLL_PERIODIC)
  
  adv => new_CSL_2d_advector(&
    charac, &
    num_cells_x1+1, &
    num_cells_x2+1, &
    eta1_coords = x1_mesh, &
    eta2_coords = x2_mesh)
  
  call adv%advect_2d(A1, A2, dt, input, output)
  
  err=maxval(abs(input-output))
  
  print *,'#err=',err
  if(err<1.e-15_f64)then  
    print *,'#PASSED' 
  endif

end program