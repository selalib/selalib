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

program unit_test_advection_1d_spectral
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_file_io.h"
use sll_module_advection_1d_base
use sll_module_advection_1d_spectral
use sll_boundary_condition_descriptors

!$ use omp_lib

implicit none
  
type(sll_advection_1d_base_ptr), pointer  :: adv(:)

sll_real64                            :: xmin
sll_real64                            :: xmax
sll_int32                             :: num_cells
sll_real64, dimension(:), allocatable :: input
sll_real64, dimension(:), allocatable :: output
sll_real64, dimension(:), allocatable :: solution
sll_real64                            :: dt
sll_real64                            :: a
sll_real64                            :: err
sll_int32                             :: ierr
sll_int32                             :: prank = 1
sll_int32                             :: psize = 1

xmin      = 0.0_f64
xmax      = 1.0_f64
num_cells = 32
a         = 1._f64
dt        = 0.1_f64
SLL_ALLOCATE(input(num_cells+1),    ierr)
SLL_ALLOCATE(output(num_cells+1),   ierr)
SLL_ALLOCATE(solution(num_cells+1), ierr)

!$OMP PARALLEL
!$ prank = OMP_GET_THREAD_NUM()
!$ print*, ' prank = ', prank
!$ psize = OMP_GET_NUM_THREADS()
!$ print*, ' psize = ', psize, xmin, xmax, num_cells

SLL_ALLOCATE(adv(psize),            ierr)

input = 1.0_f64 

adv(prank+1)%ptr => new_spectral_1d_advector(num_cells, xmin, xmax, SLL_PERIODIC) 

call adv(prank+1)%ptr%advect_1d_constant( a, dt, input, output)
err = maxval(abs(output-input))

print *,'# err=',err
if(err == 0.0)then  
  print *,'#PASSED' 
endif

!$OMP END PARALLEL

end program unit_test_advection_1d_spectral
