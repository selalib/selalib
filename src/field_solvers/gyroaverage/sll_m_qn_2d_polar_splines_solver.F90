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

module sll_m_qn_2d_polar_splines_solver
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_qn_2d_base, only: &
      sll_c_qn_2d_base

   use sll_m_qn_2d_polar, only: &
      sll_s_qn_2d_polar_init, &
      sll_s_precompute_gyroaverage_index, &
      sll_s_precompute_inverse_qn_matrix_polar_splines, &
      sll_t_qn_2d_polar, &
      sll_s_qn_2d_polar_solve

   implicit none

   public :: &
      sll_f_new_qn_2d_polar_splines_solver, &
      sll_s_qn_2d_polar_splines_solver_init, &
      sll_t_qn_2d_polar_splines_solver

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, extends(sll_c_qn_2d_base) :: sll_t_qn_2d_polar_splines_solver

      type(sll_t_qn_2d_polar) :: quasineutral

   contains
      procedure, pass(qn) :: init => &
         sll_s_qn_2d_polar_splines_solver_init
      procedure, pass(qn) :: precompute_qn => &
         precompute_qn_2d_polar_splines
      procedure, pass(qn) :: solve_qn => &
         solve_qn_2d_polar_splines

   end type sll_t_qn_2d_polar_splines_solver

contains

   function sll_f_new_qn_2d_polar_splines_solver( &
      eta_min, &
      eta_max, &
      Nc, &
      N_points, &
      lambda, &
      T_i) &
      result(qn)

      type(sll_t_qn_2d_polar_splines_solver), pointer :: qn
      sll_real64, intent(in) :: eta_min(2)
      sll_real64, intent(in) :: eta_max(2)
      sll_int32, intent(in)  :: Nc(2)
      sll_int32, intent(in)  :: N_points
      sll_real64, dimension(:), intent(in)  :: lambda
      sll_real64, dimension(:), intent(in)  :: T_i
      sll_int32 :: ierr

      SLL_ALLOCATE(qn, ierr)
      call qn%init( &
         eta_min, &
         eta_max, &
         Nc, &
         N_points, &
         lambda, &
         T_i)

   end function sll_f_new_qn_2d_polar_splines_solver

   subroutine sll_s_qn_2d_polar_splines_solver_init( &
      qn, &
      eta_min, &
      eta_max, &
      Nc, &
      N_points, &
      lambda, &
      T_i)

      class(sll_t_qn_2d_polar_splines_solver), intent(out) :: qn
      sll_real64, intent(in) :: eta_min(2)
      sll_real64, intent(in) :: eta_max(2)
      sll_int32, intent(in)  :: Nc(2)
      sll_int32, intent(in)  :: N_points
      sll_real64, dimension(:), intent(in)  :: lambda
      sll_real64, dimension(:), intent(in)  :: T_i

      call sll_s_qn_2d_polar_init(qn%quasineutral, eta_min, eta_max, Nc, N_points, lambda, T_i)

   end subroutine sll_s_qn_2d_polar_splines_solver_init

   subroutine precompute_qn_2d_polar_splines(qn, mu_points, mu_weights, N_mu)
      class(sll_t_qn_2d_polar_splines_solver), target :: qn
      sll_int32, intent(in) :: N_mu
      sll_real64, dimension(1:N_mu), intent(in) :: mu_points
      sll_real64, dimension(1:N_mu), intent(in) :: mu_weights
      sll_real64, dimension(1:N_mu) :: rho_points
      sll_int32 :: i

      do i = 1, N_mu
         rho_points(i) = sqrt(2._f64*mu_points(i))
      end do

      call sll_s_precompute_gyroaverage_index(qn%quasineutral, rho_points, N_mu)
      call sll_s_precompute_inverse_qn_matrix_polar_splines(qn%quasineutral, mu_points, mu_weights, N_mu)

   end subroutine precompute_qn_2d_polar_splines

   subroutine solve_qn_2d_polar_splines(qn, phi)
      class(sll_t_qn_2d_polar_splines_solver), target :: qn
      sll_real64, dimension(:, :), intent(inout) :: phi

      call sll_s_qn_2d_polar_solve(qn%quasineutral, phi)

   end subroutine solve_qn_2d_polar_splines

end module sll_m_qn_2d_polar_splines_solver
