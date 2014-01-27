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

program unit_test_characteristics_2d_explicit_euler_conservative
#include "sll_working_precision.h"
use sll_module_characteristics_2d_explicit_euler_conservative
use sll_boundary_condition_descriptors

implicit none
  
  class(sll_characteristics_2d_base), pointer :: euler 
  
  sll_int32 :: Npts1
  sll_int32 :: Npts2
  sll_real64, dimension(:), allocatable :: input1
  sll_real64, dimension(:), allocatable :: input2
  sll_real64, dimension(:,:), allocatable :: output1
  sll_real64, dimension(:,:), allocatable :: output2
  sll_real64, dimension(:,:), allocatable :: A1
  sll_real64, dimension(:,:), allocatable :: A2
  !sll_int32 :: ierr
  sll_int32 :: i
  sll_int32 :: j
  sll_real64 :: dt
  sll_real64 :: err
  sll_real64 :: tmp

  
  
  
  Npts1 = 28
  Npts2 = 32
  dt = 0.1_f64
  
  
  !initialization for explicit_euler_conservative_2d
  
  euler => &
    new_explicit_euler_conservative_2d_charac(&
      Npts1, &
      Npts2, &
      SLL_SET_TO_LIMIT, &
      SLL_PERIODIC)

  




  
      
      
      

  allocate(input1(Npts1))
  allocate(input2(Npts2))
  allocate(output1(Npts1,Npts2))
  allocate(output2(Npts1,Npts2))
  allocate(A1(Npts1,Npts2))
  allocate(A2(Npts1,Npts2))
  
  do i=1,Npts1
    input1(i) = real(i-1,f64)/real(Npts1-1,f64)
  enddo

  do i=1,Npts2
    input2(i) = real(i-1,f64)/real(Npts2-1,f64)
  enddo
  
  do j=1,Npts2
    do i=1,Npts1   
      A1(i,j) = -input2(j)+0.5_f64
      A2(i,j) = input1(i)-0.5_f64
    enddo
  enddo
  !call compute_explicit_euler_conservative_2d_charac( &
  call euler%compute_characteristics( &
      A1, &
      A2, &
      dt, &
      input1, &
      input2, &
      output1, &
      output2)
      
      
  err = 0._f64
  
  do j=1,Npts2
    do i=1,Npts1   
      tmp = input1(i)-dt*A1(i,j)
      if(tmp>1)then
        tmp = 1._f64
      endif
      if(tmp<0)then
        tmp = 0._f64
      endif
      tmp=abs(tmp-output1(i,j))
      if(tmp>err)then
        err=tmp
      endif

      tmp = input2(j)-dt*A2(i,j)
      tmp = tmp-floor(tmp)
      tmp=abs(tmp-output2(i,j))
      if(tmp>err)then
        err=tmp
      endif      
    enddo
  enddo
  
  print *,'#err=',err
  
  
  
  

  if(err==0)then    
    print *, '#PASSED'
  endif

end program