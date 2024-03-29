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

!> @ingroup characteristics
!> @brief Abstract class to compute the characteristic in two dimensions
module sll_m_characteristics_2d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

   implicit none

   public :: &
      sll_f_process_outside_point_periodic, &
      sll_f_process_outside_point_set_to_limit, &
      sll_i_signature_process_outside_point, &
      sll_c_characteristics_2d_base

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! For computing the characteristics in 2d
   type, abstract :: sll_c_characteristics_2d_base
   contains
      procedure(signature_compute_characteristics_2d), deferred, pass(charac) :: &
         compute_characteristics
!    procedure(sll_i_signature_process_outside_point), deferred, pass    :: &
!      process_outside_point1
!    procedure(sll_i_signature_process_outside_point), deferred, pass    :: &
!      process_outside_point2

   end type sll_c_characteristics_2d_base

   abstract interface
      !solves eta1'(t) = A1(eta1(t),eta2(t)), eta2'(t)=A2(eta1(t),eta2(t))
      !A1 <=> A1
      !A2 <=> A2
      !dt <=> dt
      !eta1(dt) <=> input1
      !eta2(dt) <=> input2
      !eta1(0) <=> output1
      !eta2(0) <=> output2
      subroutine signature_compute_characteristics_2d( &
         charac, &
         A1, &
         A2, &
         dt, &
         input1, &
         input2, &
         output1, &
         output2)
         use sll_m_working_precision
         import sll_c_characteristics_2d_base
         class(sll_c_characteristics_2d_base) :: charac
         sll_real64, dimension(:, :), intent(in) :: A1
         sll_real64, dimension(:, :), intent(in) :: A2
         sll_real64, intent(in) :: dt
         sll_real64, dimension(:), intent(in) :: input1
         sll_real64, dimension(:), intent(in) :: input2
         sll_real64, dimension(:, :), intent(out) :: output1
         sll_real64, dimension(:, :), intent(out) :: output2

      end subroutine signature_compute_characteristics_2d
   end interface

   abstract interface
      ! change the value of eta when eta<=eta_min or eta>=eta_max
      ! depending on boundary conditions
      function sll_i_signature_process_outside_point(eta, eta_min, eta_max) result(eta_out)
         use sll_m_working_precision
         sll_real64, intent(in)  :: eta
         sll_real64, intent(in) :: eta_min
         sll_real64, intent(in) :: eta_max
         sll_real64 :: eta_out
      end function sll_i_signature_process_outside_point
   end interface

contains

   ! periodic case
   ! called when bc_type = sll_p_periodic
   function sll_f_process_outside_point_periodic(eta, eta_min, eta_max) result(eta_out)
      use sll_m_working_precision
      sll_real64, intent(in)  :: eta
      sll_real64, intent(in) :: eta_min
      sll_real64, intent(in) :: eta_max
      sll_real64 :: eta_out

      eta_out = (eta - eta_min)/(eta_max - eta_min)
      eta_out = eta_out - real(floor(eta_out), kind=f64)
      if (eta_out == 1._f64) then
         eta_out = 0._f64
      end if
      if (.not. ((eta_out >= 0) .and. (eta_out < 1))) then
         print *, '#eta=', eta
         print *, '#eta_min=', eta_min
         print *, '#eta_max=', eta_max
         print *, '#(eta-eta_min)/(eta_max-eta_min)=', (eta - eta_min)/(eta_max - eta_min)
         print *, '#floor(-1e-19)', floor(-1e-19)
         print *, '#eta_out=', eta_out
      end if
      SLL_ASSERT((eta_out >= 0) .and. (eta_out < 1))
      eta_out = eta_min + eta_out*(eta_max - eta_min)
      SLL_ASSERT((eta_out >= eta_min) .and. (eta_out < eta_max))

   end function sll_f_process_outside_point_periodic

   ! set to limit case
   ! called when bc_type = SLL_SET_TO_LIMIT

   function sll_f_process_outside_point_set_to_limit(eta, eta_min, eta_max) result(eta_out)
      use sll_m_working_precision
      sll_real64, intent(in)  :: eta
      sll_real64, intent(in) :: eta_min
      sll_real64, intent(in) :: eta_max
      sll_real64 :: eta_out

      eta_out = (eta - eta_min)/(eta_max - eta_min)
      if (eta_out > 1) then
         eta_out = 1._f64
      end if
      if (eta_out < 0) then
         eta_out = 0._f64
      end if
      SLL_ASSERT((eta_out >= 0) .and. (eta_out <= 1))
      eta_out = eta_min + eta_out*(eta_max - eta_min)
      SLL_ASSERT((eta_out >= eta_min) .and. (eta_out <= eta_max))

   end function sll_f_process_outside_point_set_to_limit

end module sll_m_characteristics_2d_base
