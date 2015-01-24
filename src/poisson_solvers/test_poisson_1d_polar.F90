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

program unit_test_poisson_1d_polar_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_module_poisson_1d_polar_solver
!use sll_boundary_condition_descriptors

implicit none
  
  class(sll_poisson_1d_base), pointer :: poisson 
  sll_real64 :: err
  sll_real64 :: x1_min
  sll_real64 :: x1_max
  sll_int32 :: Nc_x1
  sll_real64,dimension(:),allocatable :: E
  sll_real64,dimension(:),allocatable :: rho
  sll_real64,dimension(:),allocatable :: check
  sll_real64 :: tmp
  sll_real64 :: tmp2
  sll_real64 :: val
  sll_real64 :: r
  sll_real64 :: dr
  sll_int32 :: ierr
  sll_int32 :: i
  
  x1_min = -1_f64
  x1_max = 1._f64

  
  Nc_x1 = 32
  
  dr = (x1_max-x1_min)/real(Nc_x1,f64)
  
  SLL_ALLOCATE(E(Nc_x1+1),ierr)
  SLL_ALLOCATE(rho(Nc_x1+1),ierr)
  SLL_ALLOCATE(check(Nc_x1+1),ierr)
  
  rho = 1._f64
  
  
  poisson =>new_poisson_1d_polar_solver( &
    x1_min, &
    x1_max, &
    Nc_x1)  
  
  call poisson%compute_E_from_rho( E, rho )
  
  !we check that the result is correct
  !-rE(r) = \int_0^rsrho(s)ds
  
  err = 0._f64
  
  val = abs(x1_min*E(1))
  if(val>err)then
    err = val
  endif
  
  
  !check(i) = int_{r(1)}^{r(i+1)}srho(s)ds
  
  check(1) = 0._f64
  tmp = 0._f64
  do i=1,Nc_x1
    r = x1_min+real(i-1,f64)*dr
    tmp = tmp+0.5_f64*(r*rho(i)+(r+dr)*rho(i+1))*dr
    check(i+1) = -tmp 
    if(abs(r+dr)<1e-12)then
      print *,'#val=',tmp
      tmp2 = -tmp
    endif   
  enddo
  
  check = check - tmp2
  
!      val = abs((r+dr)*E(i+1)-check(i))
!    if(val>err)then
!      err = val      
!    endif


  print *,'#err=',err
  
  if(err>-1e-12)then
    do i=1,Nc_x1+1
      r = x1_min+real(i-1,f64)*dr
      print *,r,r*E(i),check(i)
    enddo
  endif
  

  if(err==0)then    
    print *, '#PASSED'
  endif

end program