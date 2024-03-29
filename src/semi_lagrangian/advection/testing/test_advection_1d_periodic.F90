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

program test_advection_1d_periodic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use sll_m_advection_1d_base, only: &
      sll_c_advector_1d

   use sll_m_advection_1d_periodic, only: &
      sll_f_new_periodic_1d_advector

   use sll_m_periodic_interp, only: &
      sll_p_spline

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   class(sll_c_advector_1d), pointer :: adv
   sll_real64 :: xmin
   sll_real64 :: xmax
   sll_int32 :: num_cells
   sll_real64, dimension(:), allocatable :: input
   sll_real64, dimension(:), allocatable :: output
   sll_real64 :: dt
   sll_real64 :: A
   sll_int32 :: order
   sll_real64 :: err

   xmin = 0._f64
   xmax = 1._f64

   num_cells = 32
   A = 1._f64
   dt = 0.1_f64
   order = 4
!  order = 18

   allocate (input(num_cells + 1))
   allocate (output(num_cells + 1))

   input = 1._f64

   adv => sll_f_new_periodic_1d_advector( &
          num_cells, &
          xmin, &
          xmax, &
          sll_p_spline, &
          order)

   call adv%advect_1d_constant( &
      A, &
      dt, &
      input, &
      output)

   err = maxval(abs(input - output))

   print *, '#err=', err
   if (err == 0) then
      print *, '#PASSED'
   end if

end program test_advection_1d_periodic
