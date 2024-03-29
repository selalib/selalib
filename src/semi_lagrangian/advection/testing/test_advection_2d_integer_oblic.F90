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

program test_advection_2d_integer_oblic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_advection_1d_base, only: &
      sll_c_advector_1d

   use sll_m_advection_1d_periodic, only: &
      sll_f_new_periodic_1d_advector

   use sll_m_advection_2d_integer_oblic, only: &
      sll_t_integer_oblic_2d_advector, &
      sll_s_integer_oblic_advect_2d, &
      sll_f_new_integer_oblic_2d_advector

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic

   use sll_m_characteristics_2d_base, only: &
      sll_c_characteristics_2d_base

   use sll_m_characteristics_2d_explicit_euler, only: &
      sll_f_new_explicit_euler_2d_charac

   use sll_m_cubic_spline_interpolator_2d, only: &
      sll_t_cubic_spline_interpolator_2d

   use sll_m_interpolators_2d_base, only: &
      sll_c_interpolator_2d

   use sll_m_periodic_interp, only: &
      sll_p_spline

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type(sll_t_integer_oblic_2d_advector), pointer :: adv
   class(sll_c_advector_1d), pointer :: adv_x1
   class(sll_c_advector_1d), pointer :: adv_aligned
   class(sll_c_interpolator_2d), pointer :: interp
   type(sll_t_cubic_spline_interpolator_2d), target :: interp_cs2d
   class(sll_c_characteristics_2d_base), pointer :: charac
   sll_real64 :: x1_min
   sll_real64 :: x1_max
   sll_real64 :: x2_min
   sll_real64 :: x2_max
   sll_int32 :: num_cells_x1
   sll_int32 :: num_cells_x2
   sll_real64, dimension(:, :), allocatable :: input
   sll_real64, dimension(:, :), allocatable :: output
   sll_real64, dimension(:), pointer :: x1_mesh
   sll_real64, dimension(:), pointer :: x2_mesh
   sll_real64 :: dt
   sll_real64, dimension(:, :), allocatable :: phi
   sll_real64 :: err
   sll_int32 :: ierr
   sll_int32 :: i
   sll_real64 :: delta_x1
   sll_real64 :: delta_x2
   sll_int32 :: shift

   x1_min = 0._f64
   x1_max = 1._f64
   x2_min = 0._f64
   x2_max = 1._f64
   num_cells_x1 = 32
   num_cells_x2 = 32
   dt = 0.1_f64
   shift = 0

   delta_x1 = (x1_max - x1_min)/real(num_cells_x1, f64)
   delta_x2 = (x2_max - x2_min)/real(num_cells_x2, f64)
   SLL_ALLOCATE(x1_mesh(num_cells_x1 + 1), ierr)
   SLL_ALLOCATE(x2_mesh(num_cells_x2 + 1), ierr)
   SLL_ALLOCATE(input(num_cells_x1 + 1, num_cells_x2 + 1), ierr)
   SLL_ALLOCATE(output(num_cells_x1 + 1, num_cells_x2 + 1), ierr)
   SLL_ALLOCATE(phi(num_cells_x1 + 1, num_cells_x2 + 1), ierr)

   do i = 1, num_cells_x1 + 1
      x1_mesh(i) = x1_min + real(i - 1, f64)*delta_x1
   end do

   do i = 1, num_cells_x2 + 1
      x2_mesh(i) = x2_min + real(i - 1, f64)*delta_x2
   end do

   input = 1._f64

   phi = 1._f64

   err = 0._f64

   adv_x1 => sll_f_new_periodic_1d_advector( &
             num_cells_x1, &
             x1_min, &
             x1_max, &
             sll_p_spline, &
             4)

   adv_aligned => sll_f_new_periodic_1d_advector( &
                  num_cells_x1, &
                  x1_min, &
                  x1_max, &
                  sll_p_spline, &
                  4)

   call interp_cs2d%init( &
      num_cells_x1 + 1, &
      num_cells_x2 + 1, &
      x1_min, &
      x1_max, &
      x2_min, &
      x2_max, &
      sll_p_periodic, &
      sll_p_periodic)

   interp => interp_cs2d

   charac => sll_f_new_explicit_euler_2d_charac( &
             num_cells_x1 + 1, &
             num_cells_x2 + 1, &
             sll_p_periodic, &
             sll_p_periodic)

   adv => sll_f_new_integer_oblic_2d_advector( &
          adv_x1, &
          adv_aligned, &
          interp, &
          charac, &
          num_cells_x1 + 1, &
          num_cells_x2 + 1, &
          eta1_coords=x1_mesh, &
          eta2_coords=x2_mesh)

   call sll_s_integer_oblic_advect_2d(adv, phi, shift, dt, input, output)

   err = maxval(abs(input - output))

   print *, '#err=', err
   if (err < 1.e-15_f64) then
      print *, '#PASSED'
   end if

end program test_advection_2d_integer_oblic
