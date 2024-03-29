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

program test_characteristics_1d_trapezoid_conservative
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic

   use sll_m_characteristics_1d_base, only: &
      sll_c_characteristics_1d_base

   use sll_m_characteristics_1d_trapezoid_conservative, only: &
      sll_t_trapezoid_conservative_1d_charac

   use sll_m_cubic_spline_interpolator_1d, only: &
      sll_f_new_cubic_spline_interpolator_1d

   use sll_m_interpolators_1d_base, only: &
      sll_c_interpolator_1d

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   class(sll_c_characteristics_1d_base), pointer :: charac
   type(sll_t_trapezoid_conservative_1d_charac), target  :: trap

   sll_int32 :: Npts
   sll_real64, dimension(:), allocatable :: input
   sll_real64, dimension(:), allocatable :: output
   sll_real64, dimension(:), allocatable :: A
   !sll_int32 :: ierr
   sll_int32 :: i
   sll_real64 :: dt
   sll_real64 :: err
   class(sll_c_interpolator_1d), pointer   :: A_interp

   Npts = 32
   dt = 0.1_f64

   !initialization for verlet
   A_interp => sll_f_new_cubic_spline_interpolator_1d( &
               Npts, &
               0._f64, &
               1._f64, &
               sll_p_periodic)

   call trap%init( &
      Npts, &
      A_interp, &
      bc_type=sll_p_periodic)

   charac => trap

   allocate (input(Npts))
   allocate (output(Npts))
   allocate (A(Npts))

   do i = 1, Npts
      input(i) = real(i - 1, f64)/real(Npts - 1, f64)
   end do

   do i = 1, Npts
      A(i) = 1._f64 !-input(i)+0.5_f64
   end do

   err = 0._f64

   call charac%compute_characteristics( &
      A, &
      dt, &
      input, &
      output)

   if (err == 0) then
      print *, '#PASSED'
   end if

end program test_characteristics_1d_trapezoid_conservative
