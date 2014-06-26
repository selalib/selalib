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

program unit_test_derivative_2d_oblic
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_module_derivative_2d_oblic
use sll_module_advection_1d_periodic
implicit none

  type(oblic_2d_derivative), pointer :: deriv
  sll_int32 :: Nc_x1
  sll_real64 :: x1_min
  sll_real64 :: x1_max
  class(sll_advection_1d_base), pointer :: adv_x1
  sll_int32 :: Nc_x2
  sll_real64 :: x2_min
  sll_real64 :: x2_max
  sll_int32 :: stencil_r
  sll_int32 :: stencil_s
  sll_real64, dimension(:,:), allocatable :: input
  sll_real64, dimension(:,:), allocatable :: output
  sll_int32 :: ierr
  sll_real64 :: err
  sll_real64 :: iota
  sll_real64 :: A1  
  sll_real64 :: A2
  sll_real64 :: dt  
  
  Nc_x1 = 256
  x1_min = 0._f64
  x1_max = 1._f64
  !(r,s) = (-1,2) for LAG3 interpolation
  Nc_x2 = 32
  x2_min = 0._f64
  x2_max = 1._f64
  stencil_r = -1
  stencil_s = 2
  
  !iota = 0.43 ! !A1/A2
  dt = 0.1
  A1 = 1._f64
  A2 = 2._f64
  
  adv_x1 => new_periodic_1d_advector( &
    Nc_x1, &
    x1_min, &
    x1_max, &
    SPLINE, & 
    4) 
    
  deriv => new_oblic_2d_derivative( &
    Nc_x1, &
    adv_x1, &
    Nc_x2, &
    x2_min, &
    x2_max, &
    stencil_r, &
    stencil_s )
  
  print *,'#oblic derivative is initialized'
  
  SLL_ALLOCATE(input(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(output(Nc_x1+1,Nc_x2+1),ierr)
  
  err = 0._f64
  
  input = 1._f64
  call compute_oblic_derivative_2d( &
    deriv, &
    !iota, &
    A1, &
    A2, &
    input, &
    output)
  
  err = maxval(abs(output))

  print *,'#err=',err
  if(err<1.e-15_f64)then  
    print *,'#PASSED' 
  endif


end program