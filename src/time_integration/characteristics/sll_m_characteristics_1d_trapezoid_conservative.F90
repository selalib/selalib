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

!http://en.wikipedia.org/wiki/Trapezoidal_rule_%28differential_equations%29

!> @ingroup characteristics
!> @brief conservative version of trapezoid
module sll_m_characteristics_1d_trapezoid_conservative
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

   use sll_m_interpolators_1d_base, only: &
      sll_c_interpolator_1d

   implicit none

   public :: &
      sll_t_trapezoid_conservative_1d_charac, &
      sll_f_new_trapezoid_conservative_1d_charac

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, extends(sll_c_characteristics_1d_base) :: sll_t_trapezoid_conservative_1d_charac
      sll_int32                               :: Npts
      sll_real64                              :: eta_min
      sll_real64                              :: eta_max
      procedure(sll_i_signature_process_outside_point_1d), pointer, nopass    :: &
         process_outside_point
      class(sll_c_interpolator_1d), pointer               :: A_interp
      sll_int32 :: maxiter
      sll_real64 :: tol
      sll_int32 :: bc_type
   contains
      procedure, pass(charac) :: init => &
         initialize_trapezoid_conservative_1d_charac
      procedure, pass(charac) :: compute_characteristics => &
         compute_trapezoid_conservative_1d_charac
   end type sll_t_trapezoid_conservative_1d_charac

contains
   function sll_f_new_trapezoid_conservative_1d_charac( &
      Npts, &
      A_interp, &
      bc_type, &
      eta_min, &
      eta_max, &
      process_outside_point, &
      maxiter, &
      tol) &
      result(charac)

      type(sll_t_trapezoid_conservative_1d_charac), pointer :: charac
      sll_int32, intent(in) :: Npts
      sll_int32, intent(in), optional :: bc_type
      sll_real64, intent(in), optional  :: eta_min
      sll_real64, intent(in), optional  :: eta_max
      procedure(sll_i_signature_process_outside_point_1d), optional    :: &
         process_outside_point
      class(sll_c_interpolator_1d), target :: A_interp
      sll_int32, intent(in), optional :: maxiter
      sll_real64, intent(in), optional :: tol
      sll_int32 :: ierr

      SLL_ALLOCATE(charac, ierr)
      call initialize_trapezoid_conservative_1d_charac( &
         charac, &
         Npts, &
         A_interp, &
         bc_type, &
         eta_min, &
         eta_max, &
         process_outside_point, &
         maxiter, &
         tol)

   end function sll_f_new_trapezoid_conservative_1d_charac
   subroutine initialize_trapezoid_conservative_1d_charac( &
      charac, &
      Npts, &
      A_interp, &
      bc_type, &
      eta_min, &
      eta_max, &
      process_outside_point, &
      maxiter, &
      tol)

      class(sll_t_trapezoid_conservative_1d_charac) :: charac
      sll_int32, intent(in) :: Npts
      sll_int32, intent(in), optional :: bc_type
      sll_real64, intent(in), optional  :: eta_min
      sll_real64, intent(in), optional  :: eta_max
      procedure(sll_i_signature_process_outside_point_1d), optional    :: &
         process_outside_point
      class(sll_c_interpolator_1d), target :: A_interp
      sll_int32, intent(in), optional :: maxiter
      sll_real64, intent(in), optional :: tol

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
         print *, '#in initialize_trapezoid_conservative_1d_charac'
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
            print *, '#in initialize_trapezoid_conservative_1d_charac'
            stop
         end select
      end if

      if ((present(process_outside_point)) .and. (present(bc_type))) then
         print *, '#provide either process_outside_point or bc_type'
         print *, '#and not both'
         print *, '#in initialize_trapezoid_conservative_1d_charac'
         stop
      end if

      charac%A_interp => A_interp

      if (present(maxiter)) then
         charac%maxiter = maxiter
      else
         charac%maxiter = 1000
      end if

      if (present(tol)) then
         charac%tol = tol
      else
         charac%tol = 1.e-12_f64
      end if

   end subroutine initialize_trapezoid_conservative_1d_charac

   subroutine compute_trapezoid_conservative_1d_charac( &
      charac, &
      A, &
      dt, &
      input, &
      output)

      class(sll_t_trapezoid_conservative_1d_charac) :: charac
      sll_real64, dimension(:), intent(in) :: A
      sll_real64, intent(in) :: dt
      sll_real64, dimension(:), intent(in) ::  input
      sll_real64, dimension(:), intent(out) :: output
      sll_int32 :: j
      sll_real64 :: x2
      sll_real64 :: x2_old
      sll_real64 :: x2_i !i for inside, so that interpolation is possible
      sll_int32 :: iter
      sll_int32 :: i
      sll_int32 :: Npts
      sll_real64 :: eta_min
      sll_real64 :: eta_max
      sll_real64 :: output_min
      sll_real64 :: output_max

      Npts = charac%Npts
      eta_min = charac%eta_min
      eta_max = charac%eta_max

      SLL_ASSERT(size(A) >= Npts)
      SLL_ASSERT(size(input) >= Npts)
      SLL_ASSERT(size(output) >= Npts)

      ! [YG: 8 Aug 2016]
      ! Optional arguments 'eta_coords=input' and 'size_eta_coords=Npts'
      ! not passed because:
      !  1. most interpolators do not implement such an option
      !  2. 'input' is (usually) the same mesh used to initialize the interpolator
      call charac%A_interp%compute_interpolants(A)
!    call charac%A_interp%compute_interpolants( &
!      A, &
!      input, &
!      Npts)

      do j = 1, Npts - 1
         !We start from Y(t_{n+1})=y_j
         !and look for Y(t_n) = Yn
         !Yn = y_j-(A(y_j)+A(Yn))*dt/2
         x2 = 0.5_f64*(input(j) + input(j + 1)) - dt*A(j)
         x2_old = 0._f64
         iter = 0
         do while (iter < charac%maxiter .and. abs(x2_old - x2) > charac%tol)
            x2_old = x2
            x2_i = x2
            if ((x2 <= eta_min) .or. (x2 >= eta_max)) then
               x2_i = charac%process_outside_point(x2, eta_min, eta_max)
            else
               x2_i = x2
            end if
            x2 = 0.5_f64*(input(j) + input(j + 1)) &
                 - 0.5_f64*dt*(charac%A_interp%interpolate_from_interpolant_value(x2_i) + A(j))
         end do
         if (iter == charac%maxiter .and. abs(x2_old - x2) > charac%tol) then
            print *, '#not enough iterations for compute_trapezoid_conservative_1d_charac', &
               iter, abs(x2_old - x2)
            stop
         end if
         !if((x2<=eta_min).or.(x2>=eta_max))then
         !  x2 =  charac%process_outside_point(x2,eta_min,eta_max)
         !endif
         output(j) = x2
      end do

      select case (charac%bc_type)
      case (sll_p_periodic)
         output_min = output(Npts - 1) - (eta_max - eta_min)
         output_max = output(1) + (eta_max - eta_min)
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

   end subroutine compute_trapezoid_conservative_1d_charac

end module sll_m_characteristics_1d_trapezoid_conservative
