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

program test_characteristics_2d_verlet
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_hermite, &
      sll_p_periodic, &
      sll_p_set_to_limit

   use sll_m_characteristics_2d_base, only: &
      sll_c_characteristics_2d_base

   use sll_m_characteristics_2d_verlet, only: &
      sll_t_charac_2d_verlet

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

   class(sll_c_characteristics_2d_base), pointer :: charac
   type(sll_t_charac_2d_verlet), target  :: verlet

   sll_int32 :: Npts1
   sll_int32 :: Npts2
   sll_real64, dimension(:), allocatable :: input1
   sll_real64, dimension(:), allocatable :: input2
   sll_real64, dimension(:, :), allocatable :: output1
   sll_real64, dimension(:, :), allocatable :: output2
   sll_real64, dimension(:, :), allocatable :: A1
   sll_real64, dimension(:, :), allocatable :: A2
   !sll_int32 :: ierr
   sll_int32 :: i
   sll_int32 :: j
   sll_real64 :: dt
   sll_real64 :: err
   !sll_real64 :: tmp
   class(sll_c_interpolator_2d), pointer :: A1_interp_x1x2
   class(sll_c_interpolator_2d), pointer :: A2_interp_x1x2
   type(sll_t_cubic_spline_interpolator_2d), target  :: A1_cs2d
   type(sll_t_cubic_spline_interpolator_2d), target  :: A2_cs2d
   class(sll_c_interpolator_1d), pointer :: A1_interp_x1
   class(sll_c_interpolator_1d), pointer :: A2_interp_x1

   Npts1 = 28
   Npts2 = 32
   dt = 0.1_f64

   !initialization for verlet
   A1_interp_x1 => sll_f_new_cubic_spline_interpolator_1d( &
                   Npts1, &
                   0._f64, &
                   1._f64, &
                   sll_p_hermite)

   A2_interp_x1 => sll_f_new_cubic_spline_interpolator_1d( &
                   Npts1, &
                   0._f64, &
                   1._f64, &
                   sll_p_hermite)

   call A1_cs2d%init( &
      Npts1, &
      Npts2, &
      0.0_f64, &
      1.0_f64, &
      0.0_f64, &
      1.0_f64, &
      sll_p_hermite, &
      sll_p_periodic)

   A1_interp_x1x2 => A1_cs2d

   call A2_cs2d%init( &
      Npts1, &
      Npts2, &
      0.0_f64, &
      1.0_f64, &
      0.0_f64, &
      1.0_f64, &
      sll_p_hermite, &
      sll_p_periodic)

   A2_interp_x1x2 => A2_cs2d

   call verlet%init( &
      Npts1, &
      Npts2, &
      A1_interp_x1x2, &
      A2_interp_x1x2, &
      A1_interp_x1, &
      A2_interp_x1, &
      bc_type_1=sll_p_set_to_limit, &
      bc_type_2=sll_p_periodic)

   charac => verlet

   allocate (input1(Npts1))
   allocate (input2(Npts2))
   allocate (output1(Npts1, Npts2))
   allocate (output2(Npts1, Npts2))
   allocate (A1(Npts1, Npts2))
   allocate (A2(Npts1, Npts2))

   do i = 1, Npts1
      input1(i) = real(i - 1, f64)/real(Npts1 - 1, f64)
   end do

   do i = 1, Npts2
      input2(i) = real(i - 1, f64)/real(Npts2 - 1, f64)
   end do

   do j = 1, Npts2
      do i = 1, Npts1
         A1(i, j) = -input2(j) + 0.5_f64
         A2(i, j) = input1(i) - 0.5_f64
      end do
   end do

   err = 0._f64

   call charac%compute_characteristics( &
      A1, &
      A2, &
      dt, &
      input1, &
      input2, &
      output1, &
      output2)

   if (err == 0) then
      print *, '#PASSED'
   end if

end program test_characteristics_2d_verlet
