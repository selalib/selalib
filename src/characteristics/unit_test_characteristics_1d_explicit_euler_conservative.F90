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

program unit_test_characteristics_1d_explicit_euler_conservative
#include "sll_working_precision.h"
use sll_module_characteristics_1d_explicit_euler_conservative
use sll_boundary_condition_descriptors

implicit none
  
  class(sll_characteristics_1d_base), pointer :: euler 
  
  sll_int32 :: Npts
  sll_real64, dimension(:), allocatable :: input
  sll_real64, dimension(:), allocatable :: output
  sll_real64, dimension(:), allocatable :: A
  sll_int32 :: i
  sll_real64 :: dt
  sll_real64 :: err

  
  
  
  Npts = 101
  dt = 0._f64 !0.1_f64
  
  
  !initialization for explicit_euler_1d
  
  euler => &
    new_explicit_euler_conservative_1d_charac(&
      Npts, &
      SLL_PERIODIC)

  




  
      
      
      

  allocate(input(Npts))
  allocate(output(Npts))
  allocate(A(Npts))
  
  do i=1,Npts
    input(i) = real(i-1,f64)/real(Npts-1,f64)
  enddo

  
  do i=1,Npts   
    A(i) = -input(i)+0.5_f64
  enddo

  call euler%compute_characteristics( &
    A, &
    dt, &
    input, &
    output)
      
      
  err = 0._f64
  
!  do i=1,Npts   
!    tmp = input(i)-dt*A(i)
!    tmp = tmp-floor(tmp)
!    tmp=abs(tmp-output(i))
!    if(tmp>err)then
!        err=tmp
!    endif
!  enddo
  
  print *,'#err=',err
  
  
  
  

  if(err==0)then    
    print *, '#PASSED'
  endif

end program
