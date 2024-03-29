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

!> @ingroup interpolators
!> @brief
!> Class for the cubic spline sll_c_interpolator_2d
!> @details
!> Implements the sll_c_interpolator_2d interface
module sll_m_cubic_spline_interpolator_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_hermite, &
      sll_p_periodic

   use sll_m_cubic_splines, only: &
      sll_s_cubic_spline_2d_compute_interpolant, &
      sll_f_cubic_spline_2d_get_x1_delta, &
      sll_f_cubic_spline_2d_get_x1_max, &
      sll_f_cubic_spline_2d_get_x1_min, &
      sll_f_cubic_spline_2d_get_x2_delta, &
      sll_f_cubic_spline_2d_get_x2_max, &
      sll_f_cubic_spline_2d_get_x2_min, &
      sll_f_cubic_spline_2d_eval, &
      sll_f_cubic_spline_2d_eval_deriv_x1, &
      sll_f_cubic_spline_2d_eval_deriv_x2, &
      sll_s_cubic_spline_2d_init, &
      sll_t_cubic_spline_2d, &
      sll_s_cubic_spline_2d_free

   use sll_m_interpolators_2d_base, only: &
      sll_c_interpolator_2d

   implicit none

   public :: &
      sll_f_new_cubic_spline_interpolator_2d, &
      sll_t_cubic_spline_interpolator_2d, &
      sll_o_delete

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> @brief
   !> The spline-based interpolator is only a wrapper around the capabilities
   !> of the cubic splines.
   !> @details
   !> All interpolators share a common interface with
   !> respect to their use, as described by the interpolator_2d_base class.
   !> Where the diverse interpolators diverge is in the way to initialize them.
   type, extends(sll_c_interpolator_2d) :: sll_t_cubic_spline_interpolator_2d
      !> Number of points along first direction
      sll_int32                           :: npts1
      !> Number of points along second direction
      sll_int32                           :: npts2
      !> Cubic spline object in two dimensions
      type(sll_t_cubic_spline_2d)         :: spline
      !> Boundary condition type in first direction
      sll_int32                           :: bc_type1
      !> Boundary condition type in second direction
      sll_int32                           :: bc_type2
      !> Interpolated values
      sll_real64, dimension(:, :), pointer :: interpolation_points
   contains
      !> Initialization
      procedure, pass(interpolator) :: init => initialize_cs2d_interpolator
      !> Compute interpolants
      procedure :: compute_interpolants => compute_interpolants_cs2d
      !> Interpolate values after compute interpolants
      procedure :: interpolate_from_interpolant_value => interpolate_value_cs2d
      !> Interpolate values of first direction derivative after compute interpolants
      procedure :: interpolate_from_interpolant_derivative_eta1 => interpolate_deriv1_cs2d
      !> Interpolate values of second direction derivative after compute interpolants
      procedure :: interpolate_from_interpolant_derivative_eta2 => interpolate_deriv2_cs2d
      !> Compute interpolants and interpolate values of a 2d array
      procedure, pass :: interpolate_array => spline_interpolate2d
      !> Compute interpolants and interpolate values of a 2d array with given displacement
      procedure, pass :: interpolate_array_disp => spline_interpolate2d_disp
      !> Set spline coefficients
      procedure, pass :: set_coefficients => set_coefficients_cs2d
      !> Get spline coefficients
      procedure, pass :: get_coefficients => get_coefficients_cs2d
      !> Check if spline coefficents are set
      procedure, pass :: coefficients_are_set => coefficients_are_set_cs2d
      !> Free memory
      procedure, pass :: delete => delete_sll_cubic_spline_interpolator_2d
   end type sll_t_cubic_spline_interpolator_2d

   !> Pointer to this interpolator derived type
   type :: sll_cubic_spline_interpolator_2d_ptr
      type(sll_t_cubic_spline_interpolator_2d), pointer :: interp
   end type sll_cubic_spline_interpolator_2d_ptr

   !> Deallocate the interpolator object
   interface sll_o_delete
      module procedure delete_sll_cubic_spline_interpolator_2d
   end interface sll_o_delete

contains

   subroutine delete_sll_cubic_spline_interpolator_2d(interpolator)
      class(sll_t_cubic_spline_interpolator_2d), intent(inout) :: interpolator
      call sll_s_cubic_spline_2d_free(interpolator%spline)
   end subroutine delete_sll_cubic_spline_interpolator_2d

   !> Function that return a pointer to a cubic spline interpolator 2d object.
   !> The result can be the target of a interpolator 2d base class
   !> (sll_c_interpolator_2d)
   function sll_f_new_cubic_spline_interpolator_2d( &
      npts1, &
      npts2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      eta1_bc_type, &
      eta2_bc_type, &
      const_eta1_min_slope, &
      const_eta1_max_slope, &
      const_eta2_min_slope, &
      const_eta2_max_slope, &
      eta1_min_slopes, &
      eta1_max_slopes, &
      eta2_min_slopes, &
      eta2_max_slopes) &
      result(interpolator)

      type(sll_t_cubic_spline_interpolator_2d), pointer :: interpolator

      sll_int32, intent(in)           :: npts1
      sll_int32, intent(in)           :: npts2
      sll_real64, intent(in)           :: eta1_min
      sll_real64, intent(in)           :: eta1_max
      sll_real64, intent(in)           :: eta2_min
      sll_real64, intent(in)           :: eta2_max
      sll_int32, intent(in)           :: eta1_bc_type
      sll_int32, intent(in)           :: eta2_bc_type
      sll_real64, intent(in), optional :: const_eta1_min_slope
      sll_real64, intent(in), optional :: const_eta1_max_slope
      sll_real64, intent(in), optional :: const_eta2_min_slope
      sll_real64, intent(in), optional :: const_eta2_max_slope
      sll_real64, dimension(:), intent(in), optional :: eta1_min_slopes
      sll_real64, dimension(:), intent(in), optional :: eta1_max_slopes
      sll_real64, dimension(:), intent(in), optional :: eta2_min_slopes
      sll_real64, dimension(:), intent(in), optional :: eta2_max_slopes
      sll_int32 :: ierr

      SLL_ALLOCATE(interpolator, ierr)

      call interpolator%init( &
         npts1, &
         npts2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         eta1_bc_type, &
         eta2_bc_type, &
         const_eta1_min_slope, &
         const_eta1_max_slope, &
         const_eta2_min_slope, &
         const_eta2_max_slope, &
         eta1_min_slopes, &
         eta1_max_slopes, &
         eta2_min_slopes, &
         eta2_max_slopes)

   end function sll_f_new_cubic_spline_interpolator_2d

   ! We allow to use the enumerators of the splines module in this interpolator
   ! because:
   ! a. This is just a wrapper and is intimately related to the underlying
   !    cubic splines module.
   ! b. There is no uniform interface for the initialization anyway.
   ! The underlying implementation with the splines module could be hidden but
   ! I can't see a compelling reason why.
   subroutine initialize_cs2d_interpolator( &
      interpolator, &
      npts1, &
      npts2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      eta1_bc_type, &
      eta2_bc_type, &
      const_eta1_min_slope, &
      const_eta1_max_slope, &
      const_eta2_min_slope, &
      const_eta2_max_slope, &
      eta1_min_slopes, &
      eta1_max_slopes, &
      eta2_min_slopes, &
      eta2_max_slopes)

      class(sll_t_cubic_spline_interpolator_2d), intent(inout) :: interpolator
      sll_int32, intent(in)                         :: npts1
      sll_int32, intent(in)                         :: npts2
      sll_real64, intent(in)                        :: eta1_min
      sll_real64, intent(in)                        :: eta1_max
      sll_real64, intent(in)                        :: eta2_min
      sll_real64, intent(in)                        :: eta2_max
      sll_int32, intent(in)                          :: eta1_bc_type
      sll_int32, intent(in)                         :: eta2_bc_type
      sll_real64, intent(in), optional              :: const_eta1_min_slope
      sll_real64, intent(in), optional              :: const_eta1_max_slope
      sll_real64, intent(in), optional              :: const_eta2_min_slope
      sll_real64, intent(in), optional              :: const_eta2_max_slope
      sll_real64, dimension(:), intent(in), optional :: eta1_min_slopes
      sll_real64, dimension(:), intent(in), optional :: eta1_max_slopes
      sll_real64, dimension(:), intent(in), optional :: eta2_min_slopes
      sll_real64, dimension(:), intent(in), optional :: eta2_max_slopes

      interpolator%npts1 = npts1
      interpolator%npts2 = npts2
      interpolator%bc_type1 = eta1_bc_type
      interpolator%bc_type2 = eta2_bc_type

      call sll_s_cubic_spline_2d_init( &
         interpolator%spline, &
         npts1, &
         npts2, &
         eta1_min, &
         eta1_max, &
         eta2_min, &
         eta2_max, &
         eta1_bc_type, &
         eta2_bc_type, &
         const_slope_x1_min=const_eta1_min_slope, &
         const_slope_x1_max=const_eta1_max_slope, &
         const_slope_x2_min=const_eta2_min_slope, &
         const_slope_x2_max=const_eta2_max_slope, &
         x1_min_slopes=eta1_min_slopes, &
         x1_max_slopes=eta1_max_slopes, &
         x2_min_slopes=eta2_min_slopes, &
         x2_max_slopes=eta2_max_slopes)

   end subroutine

   subroutine compute_interpolants_cs2d( &
      interpolator, &
      data_array, &
      eta1_coords, &
      size_eta1_coords, &
      eta2_coords, &
      size_eta2_coords)
      class(sll_t_cubic_spline_interpolator_2d), intent(inout) :: interpolator
      sll_real64, dimension(:, :), intent(in)           :: data_array
      sll_real64, dimension(:), intent(in), optional :: eta1_coords
      sll_real64, dimension(:), intent(in), optional :: eta2_coords
      sll_int32, intent(in), optional :: size_eta1_coords
      sll_int32, intent(in), optional :: size_eta2_coords

      if (present(eta1_coords) .or. present(eta2_coords) &
          .or. present(size_eta1_coords) .or. present(size_eta2_coords)) then
         SLL_ERROR('compute_interpolants_cs2d', 'This case is not yet implemented')
      end if

      call sll_s_cubic_spline_2d_compute_interpolant(data_array, interpolator%spline)

   end subroutine

   function interpolate_value_cs2d(interpolator, eta1, eta2) result(val)
      class(sll_t_cubic_spline_interpolator_2d), intent(in) :: interpolator
      sll_real64 :: val
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2
      val = sll_f_cubic_spline_2d_eval(eta1, eta2, interpolator%spline)
   end function

   function interpolate_deriv1_cs2d(interpolator, eta1, eta2) result(val)
      class(sll_t_cubic_spline_interpolator_2d), intent(in) :: interpolator
      sll_real64 :: val
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2
      val = sll_f_cubic_spline_2d_eval_deriv_x1(eta1, eta2, interpolator%spline)
   end function

   function interpolate_deriv2_cs2d(interpolator, eta1, eta2) result(val)
      class(sll_t_cubic_spline_interpolator_2d), intent(in) :: interpolator
      sll_real64 :: val
      sll_real64, intent(in) :: eta1
      sll_real64, intent(in) :: eta2

      val = sll_f_cubic_spline_2d_eval_deriv_x2(eta1, eta2, interpolator%spline)

   end function

   subroutine spline_interpolate2d(this, num_points1, num_points2, data_in, &
                                   eta1, eta2, data_out)

      class(sll_t_cubic_spline_interpolator_2d), intent(inout) :: this
      sll_int32, intent(in)                           :: num_points1
      sll_int32, intent(in)                           :: num_points2
      sll_real64, dimension(:, :), intent(in)           :: eta1
      sll_real64, dimension(:, :), intent(in)           :: eta2
      sll_real64, dimension(:, :), intent(in)           :: data_in
      sll_real64, intent(out)          :: data_out(num_points1, num_points2)
      ! local variables
      sll_int32 :: i, j
      ! compute the interpolating spline coefficients
      call sll_s_cubic_spline_2d_compute_interpolant(data_in, this%spline)
      do j = 1, num_points2
      do i = 1, num_points1
         data_out(i, j) = this%interpolate_from_interpolant_value(eta1(i, j), eta2(i, j))
      end do
      end do

   end subroutine spline_interpolate2d

   subroutine spline_interpolate2d_disp(this, &
                                        num_points1, &
                                        num_points2, &
                                        data_in, &
                                        alpha1, &
                                        alpha2, &
                                        data_out)

      class(sll_t_cubic_spline_interpolator_2d), intent(inout) :: this

      sll_int32, intent(in)                         :: num_points1
      sll_int32, intent(in)                         :: num_points2
      sll_real64, dimension(:, :), intent(in)         :: alpha1
      sll_real64, dimension(:, :), intent(in)         :: alpha2
      sll_real64, dimension(:, :), intent(in)         :: data_in
      sll_real64, intent(out)        :: data_out(num_points1, num_points2)
      sll_real64                                     :: eta1
      sll_real64                                     :: eta1_min
      sll_real64                                     :: eta1_max
      sll_real64                                     :: delta_eta1
      sll_real64                                     :: eta2
      sll_real64                                     :: eta2_min
      sll_real64                                     :: eta2_max
      sll_real64                                     :: delta_eta2
      sll_int32                                      :: i
      sll_int32                                      :: j

      eta1_min = sll_f_cubic_spline_2d_get_x1_min(this%spline) !this%spline%x1_min
      eta1_max = sll_f_cubic_spline_2d_get_x1_max(this%spline) !this%spline%x1_max
      eta2_min = sll_f_cubic_spline_2d_get_x2_min(this%spline) !this%spline%x2_min
      eta2_max = sll_f_cubic_spline_2d_get_x2_max(this%spline) !this%spline%x2_max
      delta_eta1 = sll_f_cubic_spline_2d_get_x1_delta(this%spline) !this%spline%x1_delta
      delta_eta2 = sll_f_cubic_spline_2d_get_x2_delta(this%spline) !this%spline%x2_delta

      call sll_s_cubic_spline_2d_compute_interpolant(data_in, this%spline)

      if (this%bc_type1 == sll_p_periodic .and. &
          this%bc_type2 == sll_p_periodic) then

         do j = 1, num_points2
            do i = 1, num_points1
               eta1 = eta1_min + (i - 1)*delta_eta1
               eta2 = eta2_min + (j - 1)*delta_eta2
               eta1 = eta1_min + &
                      modulo(eta1 - eta1_min + alpha1(i, j), eta1_max - eta1_min)
               eta2 = eta2_min + &
                      modulo(eta2 - eta2_min + alpha2(i, j), eta2_max - eta2_min)
               data_out(i, j) = this%interpolate_from_interpolant_value(eta1, eta2)
            end do
         end do

      else if (this%bc_type1 == sll_p_hermite .and. &
               this%bc_type2 == sll_p_hermite) then

         do j = 1, num_points2
            do i = 1, num_points1
               eta1 = eta1_min + (i - 1)*delta_eta1 + alpha1(i, j)
               eta2 = eta2_min + (j - 1)*delta_eta2 + alpha2(i, j)
               eta1 = min(eta1, eta1_max)
               eta2 = min(eta2, eta2_max)
               eta1 = max(eta1, eta1_min)
               eta2 = max(eta2, eta2_min)
               data_out(i, j) = this%interpolate_from_interpolant_value(eta1, eta2)
            end do
         end do

      else

         do j = 1, num_points2
            do i = 1, num_points1
               eta1 = eta1_min + (i - 1)*delta_eta1 + alpha1(i, j)
               eta2 = eta2_min + (j - 1)*delta_eta2 + alpha2(i, j)
               SLL_ASSERT(eta1_min <= eta1 .and. eta1 <= eta1_max)
               SLL_ASSERT(eta2_min <= eta2 .and. eta2 <= eta2_max)
               data_out(i, j) = this%interpolate_from_interpolant_value(eta1, eta2)
            end do
         end do
      end if
   end subroutine spline_interpolate2d_disp

   subroutine set_coefficients_cs2d( &
      interpolator, &
      coeffs_1d, &
      coeffs_2d, &
      coeff2d_size1, &
      coeff2d_size2, &
      knots1, &
      size_knots1, &
      knots2, &
      size_knots2)
      class(sll_t_cubic_spline_interpolator_2d), intent(inout) :: interpolator
      sll_real64, dimension(:), intent(in), optional :: coeffs_1d
      sll_real64, dimension(:, :), intent(in), optional :: coeffs_2d
      ! size coeffs 2D
      sll_int32, intent(in), optional :: coeff2d_size1
      sll_int32, intent(in), optional :: coeff2d_size2
      sll_real64, dimension(:), intent(in), optional   :: knots1
      sll_real64, dimension(:), intent(in), optional   :: knots2
      sll_int32, intent(in), optional :: size_knots1
      sll_int32, intent(in), optional :: size_knots2
      print *, 'set_coefficients_cs2d(): ERROR: This function has not been ', &
         'implemented yet.'
      print *, interpolator%npts1
      if (present(coeffs_1d)) then
         print *, 'coeffs_1d present but not used'
      end if
      if (present(coeffs_2d)) then
         print *, 'coeffs_2d present but not used'
      end if
      if (present(coeff2d_size1)) then
         print *, 'coeff2d_size1 present but not used'
      end if
      if (present(coeff2d_size2)) then
         print *, 'coeff2d_size2 present but not used'
      end if
      if (present(knots1)) then
         print *, 'knots1 present but not used'
      end if
      if (present(knots2)) then
         print *, 'knots2 present but not used'
      end if
      if (present(size_knots1)) then
         print *, 'size_knots1 present but not used'
      end if
      if (present(size_knots2)) then
         print *, 'size_knots2 present but not used'
      end if

      stop
   end subroutine !set_coefficients_cs2d

!!$  subroutine compute_spl_coeff_cs2d(interpolator, &
!!$       data_array, &
!!$       eta1_coords, &
!!$       size_eta1_coords, &
!!$       eta2_coords, &
!!$       size_eta2_coords )
!!$    class(sll_t_cubic_spline_interpolator_2d), intent(inout)  :: interpolator
!!$    sll_real64, dimension(:,:), intent(in)     :: data_array
!!$    sll_real64, dimension(:), intent(in),optional       :: eta1_coords
!!$    sll_real64, dimension(:), intent(in),optional       :: eta2_coords
!!$    sll_int32, intent(in), optional                     :: size_eta1_coords
!!$    sll_int32, intent(in),optional                      :: size_eta2_coords
!!$
!!$    print *, 'compute_coefficients_cs2d(): ERROR: This function has not been',&
!!$         'implemented yet.'
!!$    stop
!!$  end subroutine compute_spl_coeff_cs2d

   function get_coefficients_cs2d(interpolator)
      class(sll_t_cubic_spline_interpolator_2d), intent(in)    :: interpolator
      sll_real64, dimension(:, :), pointer            :: get_coefficients_cs2d

      print *, 'get_coefficients_cs2d(): ERROR: This function has not been ', &
         'implemented yet.'
      get_coefficients_cs2d => null()
      print *, interpolator%npts1
      stop
   end function get_coefficients_cs2d

   function coefficients_are_set_cs2d(interpolator) result(res)
      class(sll_t_cubic_spline_interpolator_2d), intent(in) :: interpolator
      logical :: res
      res = .false.
      print *, 'coefficients_are_set_cs2d(): this function has not been implemented yet.'
      print *, '#', interpolator%npts1
      !stop
   end function coefficients_are_set_cs2d

end module sll_m_cubic_spline_interpolator_2d
