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

program unit_test_characteristics_2d_verlet
#include "sll_working_precision.h"
use sll_module_characteristics_1d_trapezoid_conservative
use sll_boundary_condition_descriptors
use sll_cubic_spline_interpolator_1d

implicit none
  
  class(sll_characteristics_1d_base),pointer :: trap

  
  sll_int32 :: Npts
  sll_real64, dimension(:), allocatable :: input
  sll_real64, dimension(:), allocatable :: output
  sll_real64, dimension(:), allocatable :: A
  !sll_int32 :: ierr
  sll_int32 :: i
  sll_real64 :: dt
  sll_real64 :: err
  class(sll_interpolator_1d_base), pointer   :: A_interp

  
  
  
  Npts = 32
  dt = 0.1_f64
  
  

  
  !initialization for verlet
  A_interp => new_cubic_spline_1d_interpolator( &
    Npts, &
    0._f64, &
    1._f64, &
    SLL_PERIODIC)





  trap => &
    new_trapezoid_conservative_1d_charac(&
      Npts, &
      A_interp, &
      bc_type=SLL_PERIODIC)
                  

  allocate(input(Npts))
  allocate(output(Npts))
  allocate(A(Npts))
  

  do i=1,Npts
    input(i) = real(i-1,f64)/real(Npts-1,f64)
  enddo
  
  do i=1,Npts
      A(i) = 1._f64 !-input(i)+0.5_f64
  enddo
      
      
  err = 0._f64  
  
  call trap%compute_characteristics( &
      A, &
      dt, &
      input, &
      output)

  

  if(err==0)then    
    print *, '#PASSED'
  endif

end program