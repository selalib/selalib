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

!> @ingroup advection
module sll_m_advection_2d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   implicit none

   public :: sll_c_advector_2d

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!solves \partial_t f +A1\partial_x1 f+A2\partial_x2 f = 0
! A1 <=> A1
! A2 <=> A2
! dt <=> dt
! f(dt) <=> input
! f(0) <=> output

   type, abstract :: sll_c_advector_2d

   contains

      procedure(advect_2d_executor), deferred, pass(adv) :: advect_2d

   end type sll_c_advector_2d

   abstract interface

      subroutine advect_2d_executor(adv, A1, A2, dt, input, output)

         use sll_m_working_precision
         import sll_c_advector_2d

         class(sll_c_advector_2d)                :: adv
         sll_real64, dimension(:, :), intent(in)  :: A1
         sll_real64, dimension(:, :), intent(in)  :: A2
         sll_real64, intent(in)  :: dt
         sll_real64, dimension(:, :), intent(in)  :: input
         sll_real64, dimension(:, :), intent(out) :: output

      end subroutine advect_2d_executor

   end interface

end module sll_m_advection_2d_base
