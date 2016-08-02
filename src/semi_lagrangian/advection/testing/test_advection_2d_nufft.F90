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

program test_advection_2d_nufft
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_advection_2d_base
use sll_m_advection_2d_bsl
use sll_m_boundary_condition_descriptors
use sll_m_characteristics_2d_base
use sll_m_characteristics_2d_explicit_euler
use sll_m_nufft_interpolator_2d
use sll_m_interpolators_2d_base

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
class(sll_c_advection_2d_base),       pointer :: adv
class(sll_c_interpolator_2d),         pointer :: interp
type(sll_t_nufft_interpolator_2d),    target  :: nufft2d
class(sll_c_characteristics_2d_base), pointer :: charac
sll_real64                                    :: x1_min
sll_real64                                    :: x1_max
sll_real64                                    :: x2_min
sll_real64                                    :: x2_max
sll_int32                                     :: num_cells_x1
sll_int32                                     :: num_cells_x2
sll_real64, dimension(:,:),       allocatable :: input
sll_real64, dimension(:,:),       allocatable :: output
sll_real64, dimension(:),         pointer     :: x1_mesh
sll_real64, dimension(:),         pointer     :: x2_mesh
sll_real64                                    :: dt
sll_real64, dimension(:,:),       allocatable :: A1
sll_real64, dimension(:,:),       allocatable :: A2
sll_real64                                    :: error
sll_int32                                     :: ierr
sll_int32                                     :: i
sll_real64                                    :: delta_x1
sll_real64                                    :: delta_x2

x1_min       = 0._f64
x1_max       = 1._f64
x2_min       = 0._f64
x2_max       = 1._f64
num_cells_x1 = 32
num_cells_x2 = 32
dt           = 0.1_f64

delta_x1 = (x1_max-x1_min)/real(num_cells_x1,f64)
delta_x2 = (x2_max-x2_min)/real(num_cells_x2,f64)

SLL_ALLOCATE(x1_mesh(num_cells_x1+1),ierr)
SLL_ALLOCATE(x2_mesh(num_cells_x2+1),ierr)
SLL_ALLOCATE(input(  num_cells_x1+1,num_cells_x2+1),ierr)
SLL_ALLOCATE(output( num_cells_x1+1,num_cells_x2+1),ierr)
SLL_ALLOCATE(A1(     num_cells_x1+1,num_cells_x2+1),ierr)
SLL_ALLOCATE(A2(     num_cells_x1+1,num_cells_x2+1),ierr)

do i=1,num_cells_x1+1
  x1_mesh(i) = x1_min+real(i-1,f64)*delta_x1
enddo

do i=1,num_cells_x2+1
  x2_mesh(i) = x2_min+real(i-1,f64)*delta_x2
enddo

input = 1._f64

A1 = 1._f64
A2 = 1._f64

error = 0._f64

call sll_s_nufft_interpolator_2d_init( nufft2d,        &
                                       num_cells_x1+1, &
                                       num_cells_x2+1, &
                                       x1_min,         &
                                       x1_max,         &
                                       x2_min,         &
                                       x2_max)

interp => nufft2d


charac => sll_f_new_explicit_euler_2d_charac( num_cells_x1+1, &
                                              num_cells_x2+1, &
                                              sll_p_periodic, &
                                              sll_p_periodic)

adv => sll_f_new_bsl_2d_advector( interp,                &
                                  charac,                &
                                  num_cells_x1+1,        &
                                  num_cells_x2+1,        &
                                  eta1_coords = x1_mesh, &
                                  eta2_coords = x2_mesh)

call adv%advect_2d(A1, A2, dt, input, output)

error = maxval(abs(input-output))

print*, '#err=', error

if (error < 1.e-15_f64) print *, '#PASSED' 

end program test_advection_2d_nufft
