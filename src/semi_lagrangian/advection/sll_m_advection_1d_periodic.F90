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

! for the moment mimic of sll_m_periodic_interpolator_1d.F90

!> @ingroup advection
module sll_m_advection_1d_periodic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_advection_1d_base, only: sll_c_advector_1d

   use sll_m_periodic_interp, only: &
      sll_s_periodic_interp_init, &
      sll_s_periodic_interp, &
      sll_t_periodic_interp_work

   implicit none

   public :: sll_f_new_periodic_1d_advector
   public :: sll_t_advector_1d_periodic

   private

   type, extends(sll_c_advector_1d) :: sll_t_advector_1d_periodic

      sll_int32                         :: num_cells
      sll_real64                        :: xmin
      sll_real64                        :: xmax
      type(sll_t_periodic_interp_work)  :: per_interp

   contains

      procedure, pass(adv) :: init => initialize_periodic_1d_advector
      procedure, pass(adv) :: advect_1d_constant => periodic_advect_1d_constant
      procedure, pass(adv) :: advect_1d => periodic_advect_1d_fake
      procedure, pass(adv) :: delete => delete_periodic_1d_advector
   end type sll_t_advector_1d_periodic

contains

   function sll_f_new_periodic_1d_advector( &
      num_cells, &
      xmin, &
      xmax, &
      type, &
      order) &
      result(adv)
      type(sll_t_advector_1d_periodic), pointer :: adv
      sll_int32, intent(in)               :: num_cells
      sll_real64, intent(in)               :: xmin
      sll_real64, intent(in)               :: xmax
      sll_int32, intent(in)               :: type
      sll_int32, intent(in)               :: order
      sll_int32 :: ierr

      SLL_ALLOCATE(adv, ierr)
      call initialize_periodic_1d_advector(adv, num_cells, xmin, xmax, type, order)

   end function sll_f_new_periodic_1d_advector

   subroutine initialize_periodic_1d_advector( &
      adv, &
      num_cells, &
      xmin, &
      xmax, &
      type, &
      order)

      class(sll_t_advector_1d_periodic) :: adv
      sll_int32, intent(in)            :: num_cells
      sll_real64, intent(in)            :: xmin
      sll_real64, intent(in)            :: xmax
      sll_int32, intent(in)            :: type
      sll_int32, intent(in)            :: order

      call sll_s_periodic_interp_init(adv%per_interp, num_cells, type, order)

      adv%num_cells = num_cells
      adv%xmin = xmin
      adv%xmax = xmax

   end subroutine initialize_periodic_1d_advector

   subroutine periodic_advect_1d_constant(adv, A, dt, input, output)

      class(sll_t_advector_1d_periodic)     :: adv
      sll_real64, intent(in)                :: A
      sll_real64, intent(in)                :: dt
      sll_real64, dimension(:), intent(in)  :: input
      sll_real64, dimension(:), intent(out) :: output

      sll_real64 :: shift
      sll_real64 :: xmin
      sll_real64 :: xmax
      sll_int32  :: num_cells

      num_cells = adv%num_cells
      xmin = adv%xmin
      xmax = adv%xmax
      shift = A*dt/(xmax - xmin)*real(num_cells, f64)

      call sll_s_periodic_interp( &
         adv%per_interp, &
         output(1:num_cells), &
         input(1:num_cells), &
         shift)

      ! complete by periodicity
      if (size(output) > num_cells) then
         output(num_cells + 1) = output(1)
      end if

   end subroutine periodic_advect_1d_constant

   subroutine periodic_advect_1d_fake(adv, A, dt, input, output)

      class(sll_t_advector_1d_periodic)     :: adv
      sll_real64, dimension(:), intent(in)  :: A
      sll_real64, intent(in)                :: dt
      sll_real64, dimension(:), intent(in)  :: input
      sll_real64, dimension(:), intent(out) :: output

      print *, '#periodic_advect_1d_fake'
      print *, '#not implemented'
      print *, maxval(A)
      print *, dt
      print *, maxval(input)
      output = 0._f64
      print *, adv%num_cells
      stop

   end subroutine periodic_advect_1d_fake

   subroutine delete_periodic_1d_advector(adv)
      class(sll_t_advector_1d_periodic), intent(inout) :: adv
      SLL_ASSERT(storage_size(adv) > 0)
   end subroutine delete_periodic_1d_advector

end module sll_m_advection_1d_periodic
