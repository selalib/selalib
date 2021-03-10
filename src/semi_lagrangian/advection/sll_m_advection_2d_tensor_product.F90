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

! in development
! use of CSL1D in 2D
! begin with a specific example (for characteristics and interpolation)

module sll_m_advection_2d_tensor_product
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_advection_1d_base, only: &
      sll_c_advector_1d

   use sll_m_advection_2d_base, only: &
      sll_c_advector_2d

   implicit none

   public :: &
      sll_f_new_advector_2d_tensor_product, &
      sll_t_advector_2d_tensor_product

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, extends(sll_c_advector_2d) :: sll_t_advector_2d_tensor_product

      class(sll_c_advector_1d), pointer :: advect_x1
      class(sll_c_advector_1d), pointer :: advect_x2
      sll_int32                         :: npts1
      sll_int32                         :: npts2
      sll_real64, dimension(:), pointer :: buf1d

   contains

      procedure, pass(adv) :: init => initialize_advector_2d_tensor_product
      procedure, pass(adv) :: advect_2d => tensor_product_advect_2d

   end type sll_t_advector_2d_tensor_product

contains

   function sll_f_new_advector_2d_tensor_product(advect_x1, advect_x2, &
                                                 npts1, npts2) result(adv)

      type(sll_t_advector_2d_tensor_product), pointer :: adv
      class(sll_c_advector_1d), pointer :: advect_x1
      class(sll_c_advector_1d), pointer :: advect_x2
      sll_int32, intent(in) :: npts1
      sll_int32, intent(in) :: npts2

      sll_int32 :: ierr

      SLL_ALLOCATE(adv, ierr)

      call adv%init(advect_x1, advect_x2, npts1, npts2)

   end function sll_f_new_advector_2d_tensor_product

   subroutine initialize_advector_2d_tensor_product( &
      adv, &
      advect_x1, &
      advect_x2, &
      npts1, &
      npts2)

      class(sll_t_advector_2d_tensor_product), intent(inout) :: adv
      class(sll_c_advector_1d), pointer :: advect_x1
      class(sll_c_advector_1d), pointer :: advect_x2
      sll_int32, intent(in) :: npts1
      sll_int32, intent(in) :: npts2

      sll_int32 :: ierr

      adv%advect_x1 => advect_x1
      adv%advect_x2 => advect_x2

      adv%npts1 = npts1
      adv%npts2 = npts2

      SLL_ALLOCATE(adv%buf1d(max(npts1, npts2)), ierr)

   end subroutine initialize_advector_2d_tensor_product

   subroutine tensor_product_advect_2d(adv, A1, A2, dt, input, output)

      class(sll_t_advector_2d_tensor_product) :: adv
      sll_real64, dimension(:, :), intent(in)  :: A1
      sll_real64, dimension(:, :), intent(in)  :: A2
      sll_real64, intent(in)  :: dt
      sll_real64, dimension(:, :), intent(in)  :: input
      sll_real64, dimension(:, :), intent(out) :: output

      sll_int32 :: i1
      sll_int32 :: i2
      sll_int32 :: npts1
      sll_int32 :: npts2

      npts1 = adv%npts1
      npts2 = adv%npts2

      do i2 = 1, npts2
         adv%buf1d(1:npts1) = input(1:npts1, i2)
         call adv%advect_x1%advect_1d( &
            A1(1:npts1, i2), &
            0.5_f64*dt, &
            adv%buf1d(1:npts1), &
            output(1:npts1, i2))
      end do

      do i1 = 1, npts1
         adv%buf1d(1:npts2) = output(i1, 1:npts2)
         call adv%advect_x2%advect_1d( &
            A2(i1, 1:npts2), &
            dt, &
            adv%buf1d(1:npts2), &
            output(i1, 1:npts2))
      end do

      do i2 = 1, npts2
         adv%buf1d(1:npts1) = output(1:npts1, i2)
         call adv%advect_x1%advect_1d( &
            A1(1:npts1, i2), &
            0.5_f64*dt, &
            adv%buf1d(1:npts1), &
            output(1:npts1, i2))
      end do

   end subroutine tensor_product_advect_2d

end module sll_m_advection_2d_tensor_product
