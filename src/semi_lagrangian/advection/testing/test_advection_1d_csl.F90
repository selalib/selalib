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

program test_advection_1d_CSL
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_advection_1d_base, only: &
      sll_c_advector_1d

   use sll_m_advection_1d_csl, only: &
      sll_f_new_csl_1d_advector

   use sll_m_advection_1d_psm, only: &
      sll_f_new_psm_1d_advector

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic

   use sll_m_characteristics_1d_base, only: &
      sll_c_characteristics_1d_base

   use sll_m_characteristics_1d_trapezoid_conservative, only: &
      sll_f_new_trapezoid_conservative_1d_charac

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_cubic_spline_interpolator_1d, only: &
      sll_f_new_cubic_spline_interpolator_1d

   use sll_m_interpolators_1d_base, only: &
      sll_c_interpolator_1d

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   class(sll_c_advector_1d), pointer :: adv
   class(sll_c_advector_1d), pointer :: adv_ref
   class(sll_c_interpolator_1d), pointer :: interp
   class(sll_c_interpolator_1d), pointer :: A_interp
   class(sll_c_characteristics_1d_base), pointer :: charac
   sll_real64 :: x_min
   sll_real64 :: x_max
   sll_real64 :: x_min_bis
   sll_real64 :: x_max_bis
   sll_int32 :: num_cells
   sll_real64, dimension(:), allocatable :: input
   sll_real64, dimension(:), allocatable :: output
   sll_real64, dimension(:), allocatable :: output_ref
   !sll_real64, dimension(:), pointer :: mesh
   sll_real64 :: dt
   sll_real64, dimension(:), allocatable :: A
   sll_real64 :: err
   sll_int32 :: ierr
   sll_int32 :: i
   sll_real64 :: delta
   sll_real64 :: x

   x_min = 0._f64
   x_max = 1._f64
   num_cells = 100
   dt = 0.01_f64 !0.1_f64

   delta = (x_max - x_min)/real(num_cells, f64)
   !SLL_ALLOCATE(mesh(num_cells+1),ierr)
   SLL_ALLOCATE(input(num_cells + 1), ierr)
   SLL_ALLOCATE(output(num_cells + 1), ierr)
   SLL_ALLOCATE(output_ref(num_cells + 1), ierr)
   SLL_ALLOCATE(A(num_cells + 1), ierr)

   !do i=1,num_cells+1
   !  mesh(i) = x_min+real(i-1,f64)*delta
   !enddo

   x_min_bis = x_min - 0.5_f64*delta
   x_max_bis = x_max - 0.5_f64*delta

   input = 1._f64

   A = 1._f64

   do i = 1, num_cells + 1
      x = x_min + real(i, f64)*delta
      A(i) = sin(2._f64*sll_p_pi*x)
      input(i) = sin(2._f64*sll_p_pi*x)
   end do

   err = 0._f64

   interp => sll_f_new_cubic_spline_interpolator_1d( &
             num_cells + 1, &
             x_min_bis, &
             x_max_bis, &
             sll_p_periodic)

   A_interp => sll_f_new_cubic_spline_interpolator_1d( &
               num_cells + 1, &
               x_min, &
               x_max, &
               sll_p_periodic)

   charac => sll_f_new_trapezoid_conservative_1d_charac( &
             num_cells + 1, &
             A_interp, &
             eta_min=x_min_bis, &
             eta_max=x_max_bis, &
             bc_type=sll_p_periodic)

   adv => sll_f_new_csl_1d_advector( &
          interp, &
          charac, &
          num_cells + 1, &
          eta_min=x_min_bis, &
          eta_max=x_max_bis, &
          bc_type=sll_p_periodic)

   adv_ref => sll_f_new_psm_1d_advector( &
              num_cells + 1, &
              eta_min=x_min, &
              eta_max=x_max)

   call adv%advect_1d(A, dt, input, output)

   call adv_ref%advect_1d(A, dt, input, output_ref)

   do i = 1, num_cells + 1
      print *, i, input(i), output(i), output_ref(i)
   end do

   err = maxval(abs(output_ref - output))

   print *, '#err=', err
   if (err < 1.e-15_f64) then
      print *, '#PASSED'
   end if

end program test_advection_1d_CSL
