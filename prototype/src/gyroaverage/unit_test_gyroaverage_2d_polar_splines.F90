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

program unit_test_gyroaverage_2d_polar_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_module_gyroaverage_2d_polar_splines_solver

implicit none
  
  class(sll_gyroaverage_2d_base), pointer :: gyroaverage 
  sll_real64 :: err
  sll_real64 :: eta_min(2)
  sll_real64 :: eta_max(2)
  sll_int32  :: Nc(2)
  sll_int32  :: N_points
  sll_int32  :: splines_case
  sll_real64 :: larmor_rad
  sll_real64,dimension(:,:),allocatable :: f
  sll_int32  :: ierr
  
  eta_min(1) = 0.1_f64
  eta_max(1) = 0.9_f64
  eta_min(2) = 0._f64
  eta_max(2) = 2._f64*sll_pi  
  
  Nc(1)=16
  Nc(2)=16
  
  SLL_ALLOCATE(f(Nc(1)+1,Nc(2)),ierr)
  
  f = 1._f64
  err = 0._f64
  larmor_rad = 0.01_f64

  N_points = 4
  
  splines_case = 1

  gyroaverage => new_gyroaverage_2d_polar_splines_solver( &
    eta_min, &
    eta_max, &
    Nc, &
    N_points, &
    splines_case)
  
  call gyroaverage%compute_gyroaverage( larmor_rad, f)

  print *,maxval(f),minval(f)

  if(err==0)then    
    print *, '#PASSED'
  endif

end program