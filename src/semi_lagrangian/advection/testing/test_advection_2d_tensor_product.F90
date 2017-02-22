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

program test_advection_2d_tensor_product
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_advection_1d_base, only: sll_c_advector_1d
use sll_m_advection_1d_bsl,  only: sll_t_advector_1d_bsl
use sll_m_advection_2d_base, only: sll_c_advector_2d
use sll_m_advection_2d_tensor_product, only: sll_t_advector_2d_tensor_product
use sll_m_boundary_condition_descriptors, only: sll_p_periodic
use sll_m_characteristics_1d_base, only: sll_c_characteristics_1d_base
use sll_m_characteristics_1d_explicit_euler, only: sll_f_new_charac_1d_explicit_euler
use sll_m_characteristics_2d_base, only: sll_c_characteristics_2d_base
use sll_m_characteristics_2d_explicit_euler, only: sll_t_charac_2d_explicit_euler
use sll_m_cubic_spline_interpolator_1d, only: sll_f_new_cubic_spline_interpolator_1d
use sll_m_cubic_spline_interpolator_2d, only: sll_t_cubic_spline_interpolator_2d
use sll_m_interpolators_1d_base, only: sll_c_interpolator_1d
use sll_m_interpolators_2d_base, only: sll_c_interpolator_2d

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
class(sll_c_advector_2d),                  pointer :: adv
class(sll_c_advector_1d),                  pointer :: adv_x1
class(sll_c_advector_1d),                  pointer :: adv_x2
class(sll_c_interpolator_1d),              pointer :: interp_x1
class(sll_c_interpolator_1d),              pointer :: interp_x2
class(sll_c_characteristics_1d_base),      pointer :: charac_x1
class(sll_c_characteristics_1d_base),      pointer :: charac_x2
class(sll_c_interpolator_2d),              pointer :: interp
type(sll_t_cubic_spline_interpolator_2d),  target  :: interp_cs2d
class(sll_c_characteristics_2d_base),      pointer :: charac

sll_real64                                         :: x1_min
sll_real64                                         :: x1_max
sll_real64                                         :: x2_min
sll_real64                                         :: x2_max
sll_int32                                          :: num_cells_x1
sll_int32                                          :: num_cells_x2
sll_real64, dimension(:,:), allocatable            :: input
sll_real64, dimension(:,:), allocatable            :: output
sll_real64, dimension(:),   pointer                :: x1_mesh
sll_real64, dimension(:),   pointer                :: x2_mesh
sll_real64                                         :: dt
sll_real64, dimension(:,:), allocatable            :: A1
sll_real64, dimension(:,:), allocatable            :: A2
sll_real64                                         :: error
sll_int32                                          :: ierr
sll_int32                                          :: i
sll_real64                                         :: delta_x1
sll_real64                                         :: delta_x2

x1_min = 0._f64
x1_max = 1._f64
x2_min = 0._f64
x2_max = 1._f64
num_cells_x1 = 32
num_cells_x2 = 32
dt = 0.1_f64

delta_x1 = (x1_max-x1_min)/real(num_cells_x1,f64)
delta_x2 = (x2_max-x2_min)/real(num_cells_x2,f64)
SLL_ALLOCATE(x1_mesh(num_cells_x1+1),ierr)
SLL_ALLOCATE(x2_mesh(num_cells_x2+1),ierr)
SLL_ALLOCATE(input(num_cells_x1+1,num_cells_x2+1),ierr)
SLL_ALLOCATE(output(num_cells_x1+1,num_cells_x2+1),ierr)
SLL_ALLOCATE(A1(num_cells_x1+1,num_cells_x2+1),ierr)
SLL_ALLOCATE(A2(num_cells_x1+1,num_cells_x2+1),ierr)

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

interp_x1 => sll_f_new_cubic_spline_interpolator_1d( num_cells_x1+1, &
  x1_min, x1_max, sll_p_periodic)
interp_x2 => sll_f_new_cubic_spline_interpolator_1d( num_cells_x2+1, &
  x2_min, x2_max, sll_p_periodic)

charac_x1 => sll_f_new_charac_1d_explicit_euler( num_cells_x1+1, sll_p_periodic)
charac_x2 => sll_f_new_charac_1d_explicit_euler( num_cells_x2+1, sll_p_periodic)

allocate( sll_t_advector_1d_bsl :: adv_x1 )
select type( a => adv_x1 )
type is (sll_t_advector_1d_bsl)
  call a%init( interp_x1, charac_x1, num_cells_x1+1, eta_coords = x1_mesh)
end select

allocate( sll_t_advector_1d_bsl :: adv_x2 )
select type( a => adv_x2 )
type is (sll_t_advector_1d_bsl)
  call a%init( interp_x2, charac_x2, num_cells_x2+1, eta_coords = x2_mesh)
end select

call interp_cs2d%init( num_cells_x1+1, num_cells_x2+1, &
  x1_min, x1_max, x2_min, x2_max, sll_p_periodic, sll_p_periodic)

interp => interp_cs2d

allocate( sll_t_charac_2d_explicit_euler :: charac )

select type ( c => charac )

type is ( sll_t_charac_2d_explicit_euler )

  call c%init( num_cells_x1+1, num_cells_x2+1, sll_p_periodic, sll_p_periodic)

end select

allocate( sll_t_advector_2d_tensor_product :: adv )

select type (a => adv )
type is( sll_t_advector_2d_tensor_product )

  call a%init( adv_x1, adv_x2, num_cells_x1+1, num_cells_x2+1)
  
end select

call adv%advect_2d(a1, a2, dt, input, output)
error = maxval(abs(input-output))

print *,'#error = ', error
if (error < 1.e-15_f64) print *,'#PASSED' 

end program test_advection_2d_tensor_product
