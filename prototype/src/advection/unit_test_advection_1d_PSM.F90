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

program unit_test_advection_1d_PSM
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_module_advection_1d_PSM

implicit none
  
  class(sll_advection_1d_base), pointer :: adv
  sll_real64 :: x_min
  sll_real64 :: x_max
  sll_real64 :: x_min_bis
  sll_real64 :: x_max_bis
  sll_int32 :: num_cells
  sll_real64, dimension(:), allocatable :: input
  sll_real64, dimension(:), allocatable :: output
  !sll_real64, dimension(:), pointer :: mesh
  sll_real64 :: dt
  sll_real64,dimension(:), allocatable :: A
  sll_real64 :: err
  sll_int32 :: ierr
  sll_int32 :: i
  sll_real64 :: delta
  
  x_min = 0._f64
  x_max = 1._f64
  num_cells = 100
  dt = 0._f64 !0.1_f64
  
  delta = (x_max-x_min)/real(num_cells,f64)
  !SLL_ALLOCATE(mesh(num_cells+1),ierr)
  SLL_ALLOCATE(input(num_cells+1),ierr)
  SLL_ALLOCATE(output(num_cells+1),ierr)
  SLL_ALLOCATE(A(num_cells+1),ierr)

  !do i=1,num_cells+1
  !  mesh(i) = x_min+real(i-1,f64)*delta
  !enddo

  x_min_bis = x_min -0.5_f64*delta
  x_max_bis = x_max -0.5_f64*delta
  
  input = 1._f64

  A = 1._f64

  err=0._f64

  
  adv => new_PSM_1d_advector(&
    num_cells+1, &
    eta_min = x_min, &
    eta_max = x_max)


  call adv%advect_1d(A, dt, input, output)
  
  do i=1,num_cells+1
    print *,i,input(i),output(i)
  enddo
  
  err=maxval(abs(input-output))
  
  print *,'#err=',err
  if(err<1.e-15_f64)then  
    print *,'#PASSED' 
  endif

end program