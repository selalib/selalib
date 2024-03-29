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
!> @brief Spectral 1d advection
!> @details Boundary conditions must be periodic.
module sll_m_advection_1d_spectral
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_advection_1d_base, only: sll_c_advector_1d

   use sll_m_constants, only: sll_p_pi

   use sll_m_fft, only: &
      sll_s_fft_init_c2r_1d, &
      sll_s_fft_init_r2c_1d, &
      sll_s_fft_exec_c2r_1d, &
      sll_s_fft_exec_r2c_1d, &
      sll_s_fft_free, &
      sll_t_fft

   implicit none

   public :: sll_f_new_advector_1d_spectral
   public :: sll_t_advector_1d_spectral

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, extends(sll_c_advector_1d) :: sll_t_advector_1d_spectral

      sll_int32                         :: num_cells
      sll_real64                        :: eta_min
      sll_real64                        :: eta_max
      sll_real64                        :: delta_eta

      sll_real64, dimension(:), pointer :: d_dx
      sll_real64, dimension(:), pointer :: kx

      type(sll_t_fft)                   :: fwx
      type(sll_t_fft)                   :: bwx
      sll_comp64, dimension(:), pointer :: fk

   contains

      procedure, pass(adv) :: init => advector_1d_spectral_init
      procedure, pass(adv) :: advect_1d
      procedure, pass(adv) :: advect_1d_constant
      procedure, pass(adv) :: delete

   end type sll_t_advector_1d_spectral

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function sll_f_new_advector_1d_spectral(num_cells, &
                                           eta_min, &
                                           eta_max &
                                           ) result(adv)

      type(sll_t_advector_1d_spectral), pointer :: adv

      sll_int32, intent(in)              :: num_cells
      sll_real64, intent(in)              :: eta_min
      sll_real64, intent(in)              :: eta_max

      sll_int32                           :: error

      SLL_ALLOCATE(adv, error)

      call adv%init(num_cells, eta_min, eta_max)

   end function sll_f_new_advector_1d_spectral

   subroutine advector_1d_spectral_init(adv, num_cells, eta_min, eta_max)

      class(sll_t_advector_1d_spectral), intent(inout) :: adv
      sll_int32, intent(in)    :: num_cells
      sll_real64, intent(in)    :: eta_min
      sll_real64, intent(in)    :: eta_max

      sll_real64     :: kx0
      sll_int32      :: i, error

      adv%eta_min = eta_min
      adv%eta_max = eta_max
      adv%num_cells = num_cells

      SLL_CLEAR_ALLOCATE(adv%d_dx(1:num_cells), error)
      SLL_ALLOCATE(adv%fk(1:num_cells/2 + 1), error)
      adv%fk(1:num_cells/2 + 1) = cmplx(0.0, 0.0, kind=f64)
!$OMP CRITICAL
      call sll_s_fft_init_r2c_1d(adv%fwx, num_cells, adv%d_dx, adv%fk)
      call sll_s_fft_init_c2r_1d(adv%bwx, num_cells, adv%fk, adv%d_dx)
!$OMP END CRITICAL

      SLL_CLEAR_ALLOCATE(adv%kx(1:num_cells/2 + 1), error)

      kx0 = 2._f64*sll_p_pi/(eta_max - eta_min)

      adv%kx(1) = 1.0_f64
      do i = 2, num_cells/2 + 1
         adv%kx(i) = real(i - 1, f64)*kx0
      end do

   end subroutine advector_1d_spectral_init

   subroutine advect_1d(adv, a, dt, input, output)

      class(sll_t_advector_1d_spectral)           :: adv
      sll_real64, dimension(:), intent(in)  :: a
      sll_real64, intent(in)  :: dt
      sll_real64, dimension(:), intent(in)  :: input
      sll_real64, dimension(:), intent(out) :: output

      character(len=*), parameter :: this_sub_name = 'advect_1d'
      sll_int32                   :: num_cells

      num_cells = adv%num_cells

      output = input
      print *, size(a), dt
      SLL_ERROR(this_sub_name, "Not implemented.")

   end subroutine advect_1d

   subroutine advect_1d_constant(adv, a, dt, input, output)

      class(sll_t_advector_1d_spectral)              :: adv
      sll_real64, intent(in)  :: a
      sll_real64, intent(in)  :: dt
      sll_real64, dimension(:), intent(in)  :: input
      sll_real64, dimension(:), intent(out) :: output

      sll_int32  :: i, num_cells

      num_cells = adv%num_cells

      adv%d_dx = input(1:num_cells)
      call sll_s_fft_exec_r2c_1d(adv%fwx, adv%d_dx, adv%fk)

      !f = f^n exp(-i kx vx dt)
      do i = 2, num_cells/2 + 1
         adv%fk(i) = adv%fk(i)*cmplx(cos(adv%kx(i)*a*dt), &
                                     -sin(adv%kx(i)*a*dt), kind=f64)
      end do

      call sll_s_fft_exec_c2r_1d(adv%bwx, adv%fk, adv%d_dx)

      output(1:num_cells) = adv%d_dx/real(num_cells, f64)

      output(num_cells + 1) = output(1)

   end subroutine advect_1d_constant

   subroutine delete(adv)

      class(sll_t_advector_1d_spectral), intent(inout) :: adv

      call sll_s_fft_free(adv%fwx)
      call sll_s_fft_free(adv%bwx)

   end subroutine delete

end module sll_m_advection_1d_spectral
