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
!> @brief computes the characteristic with explicit euler conservative scheme
module sll_m_characteristics_1d_explicit_euler_conservative
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_periodic, &
      sll_p_set_to_limit, &
      sll_p_user_defined

   use sll_m_characteristics_1d_base, only: &
      sll_f_process_outside_point_periodic, &
      sll_f_process_outside_point_set_to_limit, &
      sll_i_signature_process_outside_point_1d, &
      sll_c_characteristics_1d_base

   implicit none

   public :: &
      sll_f_new_explicit_euler_conservative_1d_charac

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, extends(sll_c_characteristics_1d_base) :: explicit_euler_conservative_1d_charac_computer
      sll_int32                               :: Npts
      sll_real64                              :: eta_min
      sll_real64                              :: eta_max
      sll_int32 :: bc_type
      procedure(sll_i_signature_process_outside_point_1d), pointer, nopass    :: &
         process_outside_point

   contains
      procedure, pass(charac) :: initialize => &
         initialize_explicit_euler_conservative_1d_charac
      procedure, pass(charac) :: compute_characteristics => &
         compute_explicit_euler_conservative_1d_charac
   end type explicit_euler_conservative_1d_charac_computer

contains
   function sll_f_new_explicit_euler_conservative_1d_charac( &
      Npts, &
      bc_type, &
      eta_min, &
      eta_max, &
      process_outside_point) &
      result(charac)

      type(explicit_euler_conservative_1d_charac_computer), pointer :: charac
      sll_int32, intent(in) :: Npts
      sll_int32, intent(in), optional :: bc_type
      sll_real64, intent(in), optional  :: eta_min
      sll_real64, intent(in), optional  :: eta_max
      procedure(sll_i_signature_process_outside_point_1d), optional    :: &
         process_outside_point
      sll_int32 :: ierr

      SLL_ALLOCATE(charac, ierr)
      call initialize_explicit_euler_conservative_1d_charac( &
         charac, &
         Npts, &
         bc_type, &
         eta_min, &
         eta_max, &
         process_outside_point)

   end function sll_f_new_explicit_euler_conservative_1d_charac

   subroutine initialize_explicit_euler_conservative_1d_charac( &
      charac, &
      Npts, &
      bc_type, &
      eta_min, &
      eta_max, &
      process_outside_point)

      class(explicit_euler_conservative_1d_charac_computer) :: charac
      sll_int32, intent(in) :: Npts
      sll_int32, intent(in), optional :: bc_type
      sll_real64, intent(in), optional  :: eta_min
      sll_real64, intent(in), optional  :: eta_max
      procedure(sll_i_signature_process_outside_point_1d), optional    :: &
         process_outside_point

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
         charac%bc_type = sll_p_user_defined
      else if (.not. (present(bc_type))) then
         print *, '#provide boundary condition'
         print *, '#bc_type or process_outside_point function'
         print *, '#in initialize_charac_1d_explicit_euler'
         stop
      else
         charac%bc_type = bc_type
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

   end subroutine initialize_explicit_euler_conservative_1d_charac

   subroutine compute_explicit_euler_conservative_1d_charac( &
      charac, &
      A, &
      dt, &
      input, &
      output)

      class(explicit_euler_conservative_1d_charac_computer) :: charac
      sll_real64, dimension(:), intent(in) :: A
      sll_real64, intent(in) :: dt
      sll_real64, dimension(:), intent(in) ::  input
      sll_real64, dimension(:), intent(out) :: output
      sll_int32 :: i
      sll_int32 :: Npts
      sll_real64 :: eta_min
      sll_real64 :: eta_max
      sll_real64 :: output_min
      sll_real64 :: output_max

      Npts = charac%Npts
      eta_min = charac%eta_min
      eta_max = charac%eta_max

      SLL_ASSERT(size(A) >= charac%Npts - 1)
      SLL_ASSERT(size(input) >= charac%Npts)
      SLL_ASSERT(size(output) >= charac%Npts)

      do i = 1, Npts - 1
         output(i) = 0.5_f64*(input(i) + input(i + 1)) - dt*A(i)
      end do
      select case (charac%bc_type)
      case (sll_p_periodic)
         output_min = output(Npts - 1) - (eta_max - eta_min)
         output_max = output(1) + (eta_max - eta_min)
         !print *,"output_min=",output_min
         !print *,"output_max=",output_max
         !stop
      case (sll_p_set_to_limit)
         output_min = 2._f64*eta_min - output(1)
         output_max = 2._f64*eta_max - output(Npts - 1)
      case default
         print *, '#bad value for charac%bc_type'
         stop
      end select

      output(Npts) = 0.5_f64*(output(Npts - 1) + output_max)

      do i = Npts - 1, 2, -1
         output(i) = 0.5_f64*(output(i) + output(i - 1))
      end do
      output(1) = 0.5_f64*(output(1) + output_min)

!print *,eta_min,eta_max
!print *,output_min,output_max
!
!do i=1,Npts
!  print *,i, input(i),output(i)
!enddo

!    do i=1,Npts
!      if((output(i)<=eta_min).or.(output(i)>=eta_max))then
!        output(i)=charac%process_outside_point(output(i),eta_min,eta_max)
!      endif
!    enddo

   end subroutine compute_explicit_euler_conservative_1d_charac

end module sll_m_characteristics_1d_explicit_euler_conservative
