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

module sll_m_gyroaverage_2d_polar_splines_solver
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_gyroaverage_2d_base, only: &
      sll_c_gyroaverage_2d_base

   use sll_m_gyroaverage_2d_polar, only: &
      sll_s_compute_gyroaverage_points_polar_spl, &
      sll_s_compute_gyroaverage_points_polar_with_invar_spl, &
      sll_s_compute_gyroaverage_pre_compute_polar_spl, &
      sll_s_compute_gyroaverage_pre_compute_polar_spl_fft, &
      sll_f_new_plan_gyroaverage_polar_splines, &
      sll_s_pre_compute_gyroaverage_polar_spl, &
      sll_s_pre_compute_gyroaverage_polar_spl_fft, &
      sll_t_plan_gyroaverage_polar

   implicit none

   public :: &
      sll_f_new_gyroaverage_2d_polar_splines_solver

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, extends(sll_c_gyroaverage_2d_base) :: gyroaverage_2d_polar_splines_solver

      type(sll_t_plan_gyroaverage_polar), pointer                   :: gyro
      sll_int32 :: splines_case
      ! splines_case
      ! 1 : splines
      ! 2 : splines with invariance
      ! 3 : splines pre-compute
      ! 4 : splines pre-compute with FFT

   contains
      procedure, pass(gyroaverage) :: initialize => &
         initialize_gyroaverage_2d_polar_splines_solver
      procedure, pass(gyroaverage) :: compute_gyroaverage => &
         compute_gyroaverage_2d_polar_splines

   end type gyroaverage_2d_polar_splines_solver

contains
   function sll_f_new_gyroaverage_2d_polar_splines_solver( &
      eta_min, &
      eta_max, &
      Nc, &
      N_points, &
      splines_case) &
      result(gyroaverage)

      type(gyroaverage_2d_polar_splines_solver), pointer :: gyroaverage
      sll_real64, intent(in) :: eta_min(2)
      sll_real64, intent(in) :: eta_max(2)
      sll_int32, intent(in)  :: Nc(2)
      sll_int32              :: N_points
      sll_int32, optional    :: splines_case
      sll_int32 :: ierr

      SLL_ALLOCATE(gyroaverage, ierr)
      call initialize_gyroaverage_2d_polar_splines_solver( &
         gyroaverage, &
         eta_min, &
         eta_max, &
         Nc, &
         N_points, &
         splines_case)

   end function sll_f_new_gyroaverage_2d_polar_splines_solver

   subroutine initialize_gyroaverage_2d_polar_splines_solver( &
      gyroaverage, &
      eta_min, &
      eta_max, &
      Nc, &
      N_points, &
      splines_case)
      class(gyroaverage_2d_polar_splines_solver), intent(out) :: gyroaverage
      sll_real64, intent(in) :: eta_min(2)
      sll_real64, intent(in) :: eta_max(2)
      sll_int32, intent(in)  :: Nc(2)
      sll_int32, intent(in)  :: N_points
      sll_int32, optional    :: splines_case

      if (.not. (present(splines_case))) then
         gyroaverage%splines_case = 1
      else
         gyroaverage%splines_case = splines_case
      end if

      select case (gyroaverage%splines_case)
      case (1, 2, 3, 4)
         gyroaverage%gyro => sll_f_new_plan_gyroaverage_polar_splines( &
                             eta_min, &
                             eta_max, &
                             Nc, &
                             N_points)
      case default
         print *, '#bad value of splines_case=', gyroaverage%splines_case
         print *, '#not implemented'
         print *, '#in initialize_gyroaverage_2d_polar_splines_solver'
         stop
      end select

   end subroutine initialize_gyroaverage_2d_polar_splines_solver

   subroutine compute_gyroaverage_2d_polar_splines(gyroaverage, larmor_rad, f)
      class(gyroaverage_2d_polar_splines_solver), target :: gyroaverage
      sll_real64, intent(in) :: larmor_rad
      sll_real64, dimension(:, :), intent(inout) :: f

      select case (gyroaverage%splines_case)
      case (1)
         call sll_s_compute_gyroaverage_points_polar_spl(gyroaverage%gyro, f, larmor_rad)
      case (2)
         call sll_s_compute_gyroaverage_points_polar_with_invar_spl(gyroaverage%gyro, f, larmor_rad)
      case (3)
         call sll_s_pre_compute_gyroaverage_polar_spl(gyroaverage%gyro, larmor_rad)
         call sll_s_compute_gyroaverage_pre_compute_polar_spl(gyroaverage%gyro, f)
      case (4)
         call sll_s_pre_compute_gyroaverage_polar_spl_fft(gyroaverage%gyro, larmor_rad)
         call sll_s_compute_gyroaverage_pre_compute_polar_spl_fft(gyroaverage%gyro, f)
      case default
         print *, '#bad value of splines_case=', gyroaverage%splines_case
         print *, '#not implemented'
         print *, 'compute_gyroaverage_2d_polar_splines'
         stop
      end select

   end subroutine compute_gyroaverage_2d_polar_splines

end module sll_m_gyroaverage_2d_polar_splines_solver
