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

!Here is an example of rotation test case

program rotation_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_advection_2d_base, only: &
      sll_c_advector_2d

   use sll_m_advection_2d_bsl, only: &
      sll_f_new_advector_2d_bsl

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_hermite, &
      sll_p_set_to_limit

   use sll_m_characteristics_2d_base, only: &
      sll_c_characteristics_2d_base

   use sll_m_characteristics_2d_verlet, only: &
      sll_f_new_verlet_2d_charac

   use sll_m_cubic_spline_interpolator_1d, only: &
      sll_f_new_cubic_spline_interpolator_1d

   use sll_m_cubic_spline_interpolator_2d, only: &
      sll_t_cubic_spline_interpolator_2d

   use sll_m_interpolators_1d_base, only: &
      sll_c_interpolator_1d

   use sll_m_interpolators_2d_base, only: &
      sll_c_interpolator_2d

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   sll_int32 :: Nc_x1
   sll_int32 :: Nc_x2
   sll_real64 :: err
   sll_int32 :: nb_step
   sll_int32 :: step
   sll_real64, dimension(:, :), allocatable :: f
   sll_real64, dimension(:, :), allocatable :: f_exact
   sll_real64, dimension(:, :), allocatable :: f_init
   sll_real64, dimension(:, :), allocatable :: f_old
   sll_real64, dimension(:, :), allocatable :: A1
   sll_real64, dimension(:, :), allocatable :: A2
   sll_int32 ::ierr
   sll_int32 ::i1
   sll_int32 ::i2
   sll_real64 :: x1
   sll_real64 :: x2
   sll_real64 :: x1_min
   sll_real64 :: x1_max
   sll_real64 :: x2_min
   sll_real64 :: x2_max
   sll_real64 :: delta_x1
   sll_real64 :: delta_x2
   sll_real64 :: dt

   class(sll_c_advector_2d), pointer :: adv
   class(sll_c_interpolator_2d), pointer :: interp
   type(sll_t_cubic_spline_interpolator_2d), target   :: interp_cs2d
   class(sll_c_characteristics_2d_base), pointer :: charac
   class(sll_c_interpolator_2d), pointer   :: A1_interp_x1x2
   class(sll_c_interpolator_2d), pointer   :: A2_interp_x1x2
   class(sll_c_interpolator_1d), pointer   :: A1_interp_x1
   class(sll_c_interpolator_1d), pointer   :: A2_interp_x1
   type(sll_t_cubic_spline_interpolator_2d), target   :: A1_cs2d
   type(sll_t_cubic_spline_interpolator_2d), target   :: A2_cs2d

   nb_step = 1000
   dt = 0.01_f64
   Nc_x1 = 64
   Nc_x2 = 64

   x1_max = 10._f64
   x1_min = -x1_max
   x2_min = x1_min
   x2_max = x1_max

   delta_x1 = (x1_max - x1_min)/real(Nc_x1, f64)
   delta_x2 = (x2_max - x2_min)/real(Nc_x2, f64)

   SLL_ALLOCATE(f(Nc_x1 + 1, Nc_x2 + 1), ierr)
   SLL_ALLOCATE(f_old(Nc_x1 + 1, Nc_x2 + 1), ierr)
   SLL_ALLOCATE(f_init(Nc_x1 + 1, Nc_x2 + 1), ierr)
   SLL_ALLOCATE(f_exact(Nc_x1 + 1, Nc_x2 + 1), ierr)
   SLL_ALLOCATE(A1(Nc_x1 + 1, Nc_x2 + 1), ierr)
   SLL_ALLOCATE(A2(Nc_x1 + 1, Nc_x2 + 1), ierr)

   !we initialize the distribution function

   do i2 = 1, Nc_x2
      do i1 = 1, Nc_x1
         x1 = x1_min + real(i1 - 1, f64)*delta_x1
         x2 = x2_min + real(i2 - 1, f64)*delta_x2
         f_init(i1, i2) = exp(-0.5_f64*(x1**2 + x2**2))
         f_exact(i1, i2) = exp(-0.5_f64*(x1**2 + x2**2))
         A1(i1, i2) = -x2
         A2(i1, i2) = x1
      end do
   end do

   f = f_init

   !we initialize the interpolator
   call interp_cs2d%init( &
      Nc_x1 + 1, &
      Nc_x2 + 1, &
      x1_min, &
      x1_max, &
      x2_min, &
      x2_max, &
      sll_p_hermite, &
      sll_p_hermite)
   interp => interp_cs2d

   !we initialize the characteristics
   A1_interp_x1 => sll_f_new_cubic_spline_interpolator_1d( &
                   Nc_x1 + 1, &
                   x1_min, &
                   x1_max, &
                   sll_p_hermite)
   A2_interp_x1 => sll_f_new_cubic_spline_interpolator_1d( &
                   Nc_x1 + 1, &
                   x1_min, &
                   x1_max, &
                   sll_p_hermite)

   call A1_cs2d%init( &
      Nc_x1 + 1, &
      Nc_x2 + 1, &
      x1_min, &
      x1_max, &
      x2_min, &
      x2_max, &
      sll_p_hermite, &
      sll_p_hermite)

   A1_interp_x1x2 => A1_cs2d

   call A2_cs2d%init( &
      Nc_x1 + 1, &
      Nc_x2 + 1, &
      x1_min, &
      x1_max, &
      x2_min, &
      x2_max, &
      sll_p_hermite, &
      sll_p_hermite)

   A2_interp_x1x2 => A2_cs2d

   charac => sll_f_new_verlet_2d_charac( &
             Nc_x1 + 1, &
             Nc_x2 + 1, &
             A1_interp_x1x2, &
             A2_interp_x1x2, &
             A1_interp_x1, &
             A2_interp_x1, &
             bc_type_1=sll_p_set_to_limit, &
             bc_type_2=sll_p_set_to_limit)

   adv => sll_f_new_advector_2d_bsl( &
          interp, &
          charac, &
          Nc_x1 + 1, &
          Nc_x2 + 1, &
          eta1_min=x1_min, &
          eta1_max=x1_max, &
          eta2_min=x2_min, &
          eta2_max=x2_max)

   err = 0._f64

   do step = 1, nb_step
      f_old = f
      call adv%advect_2d(A1, A2, dt, f_old, f)
   end do

   err = maxval(abs(f - f_old))

   print *, '#err=', err

   print *, dt, Nc_x1, Nc_x2, err

end program

