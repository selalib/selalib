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

program test_advection_1d_ampere
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_advection_1d_ampere, only: &
    sll_f_new_ampere_1d_advector

  use sll_m_advection_1d_base, only: &
    sll_t_advection_1d_base_ptr

  use sll_m_gnuplot, only: &
    sll_o_gnuplot_1d

#ifdef _OPENMP
  use omp_lib, only: &
    omp_get_num_threads, &
    omp_get_thread_num

#endif
  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
type(sll_t_advection_1d_base_ptr), pointer  :: adv(:)

sll_real64                            :: xmin, vmin
sll_real64                            :: xmax, vmax
sll_int32                             :: nc_x, nc_v
sll_real64, dimension(:), allocatable :: x
sll_real64, dimension(:), allocatable :: input
sll_real64, dimension(:), allocatable :: output
sll_real64, dimension(:), allocatable :: solution
sll_real64                            :: dt
sll_real64                            :: a
sll_real64                            :: err
sll_int32                             :: ierr
sll_int32                             :: prank = 0
sll_int32                             :: psize = 1
sll_int32                             :: istep
sll_int32                             :: i
sll_int32                             :: nstep = 100

xmin = 0.0_f64
xmax = 1.0_f64
nc_x = 32
vmin = 0.0_f64
vmax = 1.0_f64
nc_v = 32
a    = 1._f64
dt   = 0.01_f64
SLL_ALLOCATE(x(nc_x+1),    ierr)
SLL_ALLOCATE(solution(nc_x+1),    ierr)
SLL_ALLOCATE(input(nc_x+1),    ierr)
SLL_ALLOCATE(output(nc_x+1),   ierr)

do i = 1, nc_x+1
  x(i) = xmin + (i-1)*(xmax-xmin)/nc_x - 0.5
end do

!$OMP PARALLEL
!$ prank = OMP_GET_THREAD_NUM()
!$ print*, ' prank = ', prank
!$ psize = OMP_GET_NUM_THREADS()
!$ print*, ' psize = ', psize, xmin, xmax, nc_x

SLL_ALLOCATE(adv(psize),            ierr)

solution = exp(-(x*x)/0.01)
input = solution

adv(prank+1)%ptr => sll_f_new_ampere_1d_advector(nc_x, xmin, xmax  )

do istep = 1, nstep
   call adv(prank+1)%ptr%advect_1d_constant( a, dt, input, output)
   input = output
   call sll_o_gnuplot_1d(output, x, 'f_ampere', istep)
end do

err = maxval(abs(solution-input))

print *,'# err=',err
if(err <= 1e-3)then  
  print *,'PASSED' 
endif

!$OMP END PARALLEL

end program test_advection_1d_ampere
