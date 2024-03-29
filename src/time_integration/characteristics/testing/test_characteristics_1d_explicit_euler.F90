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

program test_characteristics_1d_explicit_euler
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic

   use sll_m_characteristics_1d_base, only: &
      sll_c_characteristics_1d_base

   use sll_m_characteristics_1d_explicit_euler, only: &
      sll_f_new_charac_1d_explicit_euler

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   class(sll_c_characteristics_1d_base), pointer :: euler

   sll_int32 :: Npts
   sll_real64, dimension(:), allocatable :: input
   sll_real64, dimension(:), allocatable :: output
   sll_real64, dimension(:), allocatable :: A
   sll_int32 :: i
   sll_real64 :: dt
   sll_real64 :: err
   sll_real64 :: tmp

   Npts = 28
   dt = 0.1_f64

   !initialization for explicit_euler_1d

   euler => &
      sll_f_new_charac_1d_explicit_euler( &
      Npts, &
      sll_p_periodic)

   allocate (input(Npts))
   allocate (output(Npts))
   allocate (A(Npts))

   do i = 1, Npts
      input(i) = real(i - 1, f64)/real(Npts - 1, f64)
   end do

   do i = 1, Npts
      A(i) = -input(i) + 0.5_f64
   end do

   call euler%compute_characteristics( &
      A, &
      dt, &
      input, &
      output)

   err = 0._f64

   do i = 1, Npts
      tmp = input(i) - dt*A(i)
      tmp = tmp - real(floor(tmp), f64)
      tmp = abs(tmp - output(i))
      if (tmp > err) then
         err = tmp
      end if
   end do

   print *, '#err=', err

   if (err == 0) then
      print *, '#PASSED'
   end if

end program test_characteristics_1d_explicit_euler
