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
!> computes the characteristic with explicit euler scheme
module sll_m_characteristics_1d_explicit_euler
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic, &
      sll_p_set_to_limit

   use sll_m_characteristics_1d_base, only: &
      sll_f_process_outside_point_periodic, &
      sll_f_process_outside_point_set_to_limit, &
      sll_i_signature_process_outside_point_1d, &
      sll_c_characteristics_1d_base

   implicit none

   public :: &
      sll_t_charac_1d_explicit_euler, &
      sll_f_new_charac_1d_explicit_euler

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, extends(sll_c_characteristics_1d_base) :: sll_t_charac_1d_explicit_euler
      sll_int32                               :: Npts
      sll_real64                              :: eta_min
      sll_real64                              :: eta_max
      procedure(sll_i_signature_process_outside_point_1d), pointer, nopass    :: &
         process_outside_point
      logical :: feet_inside

   contains
      procedure, pass(charac) :: init => initialize_charac_1d_explicit_euler
      procedure, pass(charac) :: compute_characteristics => &
         compute_charac_1d_explicit_euler
   end type sll_t_charac_1d_explicit_euler

contains
   function sll_f_new_charac_1d_explicit_euler( &
      Npts, &
      bc_type, &
      eta_min, &
      eta_max, &
      process_outside_point, &
      feet_inside) &
      result(charac)

      type(sll_t_charac_1d_explicit_euler), pointer :: charac
      sll_int32, intent(in) :: Npts
      sll_int32, intent(in), optional :: bc_type
      sll_real64, intent(in), optional  :: eta_min
      sll_real64, intent(in), optional  :: eta_max
      procedure(sll_i_signature_process_outside_point_1d), optional    :: &
         process_outside_point
      logical, optional :: feet_inside
      sll_int32 :: ierr

      SLL_ALLOCATE(charac, ierr)
      call initialize_charac_1d_explicit_euler( &
         charac, &
         Npts, &
         bc_type, &
         eta_min, &
         eta_max, &
         process_outside_point, &
         feet_inside)

   end function sll_f_new_charac_1d_explicit_euler

   subroutine initialize_charac_1d_explicit_euler( &
      charac, &
      Npts, &
      bc_type, &
      eta_min, &
      eta_max, &
      process_outside_point, &
      feet_inside)

      class(sll_t_charac_1d_explicit_euler) :: charac
      sll_int32, intent(in) :: Npts
      sll_int32, intent(in), optional :: bc_type
      sll_real64, intent(in), optional  :: eta_min
      sll_real64, intent(in), optional  :: eta_max
      procedure(sll_i_signature_process_outside_point_1d), optional    :: &
         process_outside_point
      logical, optional :: feet_inside

      charac%Npts = Npts

      if (present(eta_min)) then
         charac%eta_min = eta_min
      else
         charac%eta_min = 0._f64
      end if
      if (present(eta_max)) then
         charac%eta_max = eta_max
      else
         charac%eta_max = 1._f64
      end if

      if (present(process_outside_point)) then
         charac%process_outside_point => process_outside_point
      else if (.not. (present(bc_type))) then
         print *, '#provide boundary condition'
         print *, '#bc_type or process_outside_point function'
         print *, '#in initialize_charac_1d_explicit_euler'
         stop
      else
         select case (bc_type)
         case (sll_p_periodic)
            charac%process_outside_point => sll_f_process_outside_point_periodic
         case (sll_p_set_to_limit)
            charac%process_outside_point => sll_f_process_outside_point_set_to_limit
         case default
            print *, '#bad value of boundary condition'
            print *, '#in initialize_charac_1d_explicit_euler'
            stop
         end select
      end if

      if ((present(process_outside_point)) .and. (present(bc_type))) then
         print *, '#provide either process_outside_point or bc_type'
         print *, '#and not both'
         print *, '#in initialize_explicit_euler_2d_charac'
         stop
      end if

      if (present(feet_inside)) then
         charac%feet_inside = feet_inside
      else
         charac%feet_inside = .true.
      end if

   end subroutine initialize_charac_1d_explicit_euler

   subroutine compute_charac_1d_explicit_euler( &
      charac, &
      A, &
      dt, &
      input, &
      output)

      class(sll_t_charac_1d_explicit_euler) :: charac
      sll_real64, dimension(:), intent(in) :: A
      sll_real64, intent(in) :: dt
      sll_real64, dimension(:), intent(in) ::  input
      sll_real64, dimension(:), intent(out) :: output
      sll_int32 :: i
      sll_int32 :: Npts
      sll_real64 :: eta_min
      sll_real64 :: eta_max

      Npts = charac%Npts
      eta_min = charac%eta_min
      eta_max = charac%eta_max

      SLL_ASSERT(size(A) >= charac%Npts)
      SLL_ASSERT(size(input) >= charac%Npts)
      SLL_ASSERT(size(output) >= charac%Npts)

      do i = 1, Npts
         output(i) = input(i) - dt*A(i)
         if (((output(i) <= eta_min) .or. (output(i) >= eta_max)) .and. (charac%feet_inside)) then
            output(i) = charac%process_outside_point(output(i), eta_min, eta_max)
         end if
      end do

   end subroutine compute_charac_1d_explicit_euler

end module sll_m_characteristics_1d_explicit_euler
