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

program test_gyroaverage_2d_polar_hermite
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_gyroaverage_2d_base, only: &
      sll_c_gyroaverage_2d_base

   use sll_m_gyroaverage_2d_polar_hermite_solver, only: &
      sll_f_new_gyroaverage_2d_polar_hermite_solver

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   class(sll_c_gyroaverage_2d_base), pointer :: gyroaverage
   sll_real64 :: err
   sll_real64 :: eta_min(2)
   sll_real64 :: eta_max(2)
   sll_int32  :: Nc(2)
   sll_int32  :: N_points
   sll_int32  :: hermite_case
   sll_int32  :: interp_degree(2)
   sll_real64 :: larmor_rad
   sll_real64, dimension(:, :), allocatable :: f
   sll_int32  :: ierr
   sll_int32  :: i

   eta_min(1) = 0.1_f64
   eta_max(1) = 0.9_f64
   eta_min(2) = 0._f64
   eta_max(2) = 2._f64*sll_p_pi

   Nc(1) = 16
   Nc(2) = 16

   SLL_ALLOCATE(f(Nc(1) + 1, Nc(2) + 1), ierr)

   !call compute_init_f_polar(f,(/1,1/),Nc,eta_min,eta_max)

   f = 1._f64
   err = 0._f64
   larmor_rad = 0.01_f64

   N_points = 4

   interp_degree(1) = 3
   interp_degree(2) = 3

   hermite_case = 2

   gyroaverage => sll_f_new_gyroaverage_2d_polar_hermite_solver( &
                  eta_min, &
                  eta_max, &
                  Nc, &
                  N_points, &
                  interp_degree, &
                  larmor_rad, &
                  hermite_case)

   do i = 1, 10
      call gyroaverage%compute_gyroaverage(larmor_rad, f)
      print *, minval(f), maxval(f)
   end do

   if (err == 0) then
      print *, '#PASSED'
   end if

end program test_gyroaverage_2d_polar_hermite

