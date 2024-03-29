#ifndef DOXYGEN_SHOULD_SKIP_THIS
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

!> @ingroup poisson_solvers
!> solves axisymmetric poisson
module sll_m_poisson_1d_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_errors.h"

   use sll_m_poisson_1d_base, only: &
      sll_c_poisson_1d_base

   implicit none

   public :: &
      sll_f_new_poisson_1d_polar, &
      sll_t_poisson_1d_polar

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, extends(sll_c_poisson_1d_base) :: sll_t_poisson_1d_polar
      sll_real64 :: length
      sll_int32 :: nc_eta1
      !type(sll_plan_poisson_polar), pointer                   :: poiss

   contains
      procedure, pass(poisson) :: init => &
         initialize_poisson_1d_polar
      procedure, pass(poisson) :: compute_phi_from_rho => &
         compute_phi_from_rho_1d_polar
      procedure, pass(poisson) :: compute_E_from_rho => &
         compute_E_from_rho_1d_polar
!    procedure, pass(poisson) :: compute_E_from_phi => &
!      compute_E_from_phi_2d_polar

   end type sll_t_poisson_1d_polar

contains
   function sll_f_new_poisson_1d_polar( &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      bc) &
      result(poisson)

      type(sll_t_poisson_1d_polar), pointer :: poisson
      sll_real64, intent(in) :: eta1_min
      sll_real64, intent(in) :: eta1_max
      sll_int32, intent(in) :: nc_eta1
      sll_int32, intent(in), optional :: bc
      sll_int32 :: ierr

      SLL_ALLOCATE(poisson, ierr)
      call initialize_poisson_1d_polar( &
         poisson, &
         eta1_min, &
         eta1_max, &
         nc_eta1, &
         bc)

   end function sll_f_new_poisson_1d_polar

   subroutine initialize_poisson_1d_polar( &
      poisson, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      bc)
      class(sll_t_poisson_1d_polar) :: poisson
      sll_real64, intent(in) :: eta1_min
      sll_real64, intent(in) :: eta1_max
      sll_int32, intent(in) :: nc_eta1
      sll_int32, intent(in), optional :: bc
      !sll_int32 :: ierr

      if (present(bc)) then
         print *, '#Warning bc=', bc, 'present but not used'
      end if
      poisson%length = eta1_max - eta1_min
      poisson%nc_eta1 = nc_eta1

   end subroutine initialize_poisson_1d_polar

   ! solves -\Delta phi = rho in 1d
   subroutine compute_phi_from_rho_1d_polar(poisson, phi, rho)
      class(sll_t_poisson_1d_polar)       :: poisson
      sll_real64, dimension(:), intent(in)  :: rho
      sll_real64, dimension(:), intent(out) :: phi

      SLL_ERROR('compute_phi_from_rho_1d_polar', '#not implemented yet')

   end subroutine compute_phi_from_rho_1d_polar

   ! solves E = -\nabla Phi with -\Delta phi = rho in 1d
   subroutine compute_E_from_rho_1d_polar(poisson, E, rho)
      class(sll_t_poisson_1d_polar) :: poisson
      sll_real64, dimension(:), intent(in) :: rho
      sll_real64, dimension(:), intent(out) :: E
      sll_int32 :: N
      sll_real64 :: L

      N = poisson%nc_eta1
      L = poisson%length

      E(1:N + 1) = rho(1:N + 1)

      call poisson1dpolar(E, L, N)

   end subroutine compute_E_from_rho_1d_polar

   subroutine poisson1dpolar(E, L, N)
      integer, intent(in)::N
      !sll_real64,dimension(0:N),intent(inout)::E
      sll_real64, dimension(:), intent(inout) :: E
      sll_real64, intent(in) :: L
      integer :: i
      sll_real64 :: eold
      sll_real64 :: enew
      sll_real64 :: dx
      sll_real64 :: tmp

      !dx = L/real(N,f64)
      dx = L/(2._f64*real(N, f64))

      eold = E(1 + N/2)*dx
      tmp = 0._f64
      E(1 + N/2) = 0._f64
      do i = 1, N/2
         enew = E(1 + N/2 + i)*dx
         tmp = (tmp - eold)*(1._f64 - 1._f64/real(i, f64)) - enew
         eold = enew
         E(1 + N/2 + i) = tmp
         E(1 + N/2 - i) = -tmp
      end do

   end subroutine poisson1dpolar

end module sll_m_poisson_1d_polar
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
