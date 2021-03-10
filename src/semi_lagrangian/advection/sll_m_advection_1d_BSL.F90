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

! in development; should be at least cubic splines
! attached with computation of characteristics

module sll_m_advection_1d_bsl
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_errors.h"

   use sll_m_advection_1d_base, only: sll_c_advector_1d
   use sll_m_characteristics_1d_base, only: sll_c_characteristics_1d_base
   use sll_m_interpolators_1d_base, only: sll_c_interpolator_1d

   implicit none

   public :: sll_f_new_advector_1d_bsl, sll_t_advector_1d_bsl

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, extends(sll_c_advector_1d) :: sll_t_advector_1d_bsl

      class(sll_c_interpolator_1d), pointer  :: interp
      class(sll_c_characteristics_1d_base), pointer  :: charac
      sll_real64, dimension(:), pointer :: eta_coords
      sll_real64, dimension(:), pointer :: charac_feet
      sll_int32 :: npts

   contains

      procedure, pass(adv) :: init => initialize_advector_1d_bsl
      procedure, pass(adv) :: advect_1d => bsl_advect_1d
      procedure, pass(adv) :: advect_1d_constant => bsl_advect_1d_constant
      procedure, pass(adv) :: delete => delete_bsl_1d_adv

   end type sll_t_advector_1d_bsl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   function sll_f_new_advector_1d_bsl(interp, charac, npts, eta_min, eta_max, &
                                      eta_coords) result(adv)

      type(sll_t_advector_1d_bsl), pointer :: adv
      class(sll_c_interpolator_1d), pointer :: interp
      class(sll_c_characteristics_1d_base), pointer :: charac

      sll_int32, intent(in)           :: npts
      sll_real64, intent(in), optional :: eta_min
      sll_real64, intent(in), optional :: eta_max
      sll_real64, pointer, optional :: eta_coords(:)

      sll_int32 :: ierr

      SLL_ALLOCATE(adv, ierr)

      call adv%init(interp, charac, npts, eta_min, eta_max, eta_coords)

   end function sll_f_new_advector_1d_bsl

   subroutine initialize_advector_1d_bsl(adv, interp, charac, npts, &
                                         eta_min, eta_max, eta_coords)

      class(sll_t_advector_1d_bsl), intent(inout) :: adv
      class(sll_c_interpolator_1d), target        :: interp
      class(sll_c_characteristics_1d_base), target        :: charac
      sll_int32, intent(in)    :: npts

      sll_real64, intent(in), optional :: eta_min
      sll_real64, intent(in), optional :: eta_max
      sll_real64, dimension(:), pointer, optional :: eta_coords

      sll_int32  :: ierr
      sll_int32  :: i
      sll_real64 :: delta_eta

      adv%npts = npts
      adv%interp => interp
      adv%charac => charac

      SLL_ALLOCATE(adv%eta_coords(npts), ierr)
      SLL_ALLOCATE(adv%charac_feet(npts), ierr)

      if (present(eta_min) .and. present(eta_max)) then
         if (present(eta_coords)) then
            SLL_ERROR('initialize_advector_1d_bsl', 'provide either eta_coords or eta_min and eta_max and not both')
         else
            delta_eta = (eta_max - eta_min)/real(npts - 1, f64)
            do i = 1, npts
               adv%eta_coords(i) = eta_min + real(i - 1, f64)*delta_eta
            end do
         end if
      else if (present(eta_coords)) then
         if (size(eta_coords) < npts) then
            SLL_ERROR('initialize_advector_1d_bsl', 'bad size for eta_coords in initialize_sll_t_advector_1d_bsl')
         else
            adv%eta_coords(1:npts) = eta_coords(1:npts)
         end if
      else
         SLL_WARNING('initialize_advector_1d_bsl', 'we assume eta_min = 0._f64 eta_max = 1._f64')
         delta_eta = 1._f64/real(npts - 1, f64)
         do i = 1, npts
            adv%eta_coords(i) = real(i - 1, f64)*delta_eta
         end do
      end if

   end subroutine initialize_advector_1d_bsl

   subroutine bsl_advect_1d(adv, A, dt, input, output)

      class(sll_t_advector_1d_bsl)          :: adv
      sll_real64, dimension(:), intent(in)  :: A
      sll_real64, intent(in)                :: dt
      sll_real64, dimension(:), intent(in)  :: input
      sll_real64, dimension(:), intent(out) :: output

      call adv%charac%compute_characteristics(A, dt, adv%eta_coords, adv%charac_feet)

      call adv%interp%interpolate_array(adv%npts, input, adv%charac_feet, output)

   end subroutine bsl_advect_1d

   subroutine bsl_advect_1d_constant(adv, A, dt, input, output)

      class(sll_t_advector_1d_bsl)          :: adv
      sll_real64, intent(in)  :: A
      sll_real64, intent(in)  :: dt
      sll_real64, dimension(:), intent(in)  :: input
      sll_real64, dimension(:), intent(out) :: output
      sll_real64, dimension(:), allocatable :: A1
      sll_int32 :: ierr

      SLL_ALLOCATE(A1(adv%npts), ierr)

      A1 = A

      call adv%charac%compute_characteristics(A1, dt, adv%eta_coords, adv%charac_feet)

      call adv%interp%interpolate_array(adv%npts, input, adv%charac_feet, output)

      SLL_DEALLOCATE_ARRAY(A1, ierr)

   end subroutine bsl_advect_1d_constant

   subroutine delete_bsl_1d_adv(adv)

      class(sll_t_advector_1d_bsl), intent(inout) :: adv

      sll_int32 :: ierr

      SLL_DEALLOCATE(adv%eta_coords, ierr)
      SLL_DEALLOCATE(adv%charac_feet, ierr)

   end subroutine delete_bsl_1d_adv

end module sll_m_advection_1d_bsl
