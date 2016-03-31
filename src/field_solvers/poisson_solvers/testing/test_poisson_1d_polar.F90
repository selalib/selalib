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

program test_poisson_1d_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_poisson_1d_base, only: &
  sll_c_poisson_1d_base

use sll_m_poisson_1d_polar, only: &
  sll_f_new_poisson_1d_polar

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
class(sll_c_poisson_1d_base), pointer :: poisson 
sll_int32                             :: Nc_x1
sll_int32                             :: ierr
sll_int32                             :: i
sll_real64                            :: err
sll_real64                            :: x1_min
sll_real64                            :: x1_max
sll_real64                            :: tmp
sll_real64                            :: tmp2
sll_real64                            :: dr
sll_real64, dimension(:), allocatable :: r
sll_real64, dimension(:), allocatable :: E
sll_real64, dimension(:), allocatable :: rho
sll_real64, dimension(:), allocatable :: check

x1_min = -1._f64
x1_max = +1._f64

Nc_x1  = 32

dr = (x1_max-x1_min)/real(Nc_x1,f64)

SLL_ALLOCATE(r(Nc_x1+1),ierr)
SLL_ALLOCATE(E(Nc_x1+1),ierr)
SLL_ALLOCATE(rho(Nc_x1+1),ierr)
SLL_ALLOCATE(check(Nc_x1+1),ierr)

E     = 0.0_f64
rho   = 1._f64
check = 0.0_f64

poisson => sll_f_new_poisson_1d_polar( x1_min, x1_max, Nc_x1)  

call poisson%compute_E_from_rho( E, rho )

!we check that the result is correct
!-rE(r) = \int_0^rsrho(s)ds
!check(i) = int_{r(1)}^{r(i+1)}srho(s)ds

check(1) = 0._f64
tmp = 0._f64
do i=1,Nc_x1
  r(i) = x1_min+real(i-1,f64)*dr
  tmp = tmp+0.5_f64*(r(i)*rho(i)+(r(i)+dr)*rho(i+1))*dr
  check(i+1) = -tmp 
  if(abs(r(i)+dr)<1e-12)then
    print *,'#val=',tmp
    tmp2 = -tmp
  endif   
enddo

check = check - tmp2

err = maxval(abs(r(1:Nc_x1)*E(1:Nc_x1)-check(1:Nc_x1)))
  
print*, 'error = ', err
if (err>1e-12) then
  print *, '#FAILED'
else  
  print *, '#PASSED'
endif

end program
