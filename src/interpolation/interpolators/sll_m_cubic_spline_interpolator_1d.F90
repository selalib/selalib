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
!! Interpolator 1d using cubic splines on regular mesh
!! @details
!! the following provides an implementation for the abstract interface
!! sll_interpolator_1d and define spline interpolation of values in 
!! data define on original grid at points coordinates
!!
module sll_m_cubic_spline_interpolator_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_cubic_splines, only: &
    sll_s_compute_cubic_spline_1d, &
    sll_f_interpolate_derivative, &
    sll_s_interpolate_from_interpolant_array, &
    sll_s_interpolate_from_interpolant_derivatives_eta1, &
    sll_f_interpolate_from_interpolant_value, &
    sll_f_new_cubic_spline_1d, &
    sll_t_cubic_spline_1d, &
    sll_o_delete

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  implicit none

  public :: &
    sll_f_new_cubic_spline_interpolator_1d, &
    sll_t_cubic_spline_interpolator_1d, &
    sll_o_delete

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Cubic spline interpolator 1d regular grid.
type, extends(sll_c_interpolator_1d) :: sll_t_cubic_spline_interpolator_1d

   sll_real64, dimension(:), pointer  :: interpolation_points !< points position
   sll_int32                          :: num_points           !< size
   sll_int32                          :: bc_type              !< boundary condition
   type(sll_t_cubic_spline_1d), pointer :: spline               !< spline object

contains

   procedure :: initialize => initialize_cs1d_interpolator
   procedure :: compute_interpolants => compute_interpolants_cs1d
   procedure :: interpolate_from_interpolant_value => interpolate_value_cs1d
   procedure :: interpolate_from_interpolant_derivative_eta1 => interpolate_deriv1_cs1d
   procedure :: interpolate_from_interpolant_array => interpolate_values_cs1d
   procedure :: interpolate_from_interpolant_derivatives_eta1 => interpolate_derivatives_cs1d
   procedure :: interpolate_array => spline_interpolate1d
   procedure :: interpolate_array_disp => spline_interpolate1d_disp
   procedure :: interpolate_array_disp_inplace => spline_interpolate1d_disp_inplace
   procedure :: set_coefficients => set_coefficients_cs1d
   procedure :: get_coefficients => get_coefficients_cs1d

end type sll_t_cubic_spline_interpolator_1d


!> Pointer to cubic spline interpolator implementation 1D
type :: sll_cubic_spline_interpolator_1d_ptr
   type(sll_t_cubic_spline_interpolator_1d), pointer :: interp
end type sll_cubic_spline_interpolator_1d_ptr


!> Deallocate the interpolator object
interface sll_o_delete
   module procedure delete_cs1d
end interface sll_o_delete

contains  ! ****************************************************************



  subroutine spline_interpolate1d(this, num_pts, data, coordinates, output_array)

    class(sll_t_cubic_spline_interpolator_1d),  intent(in)       :: this
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(num_pts), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_pts), intent(out)  :: output_array
    
    call sll_s_compute_cubic_spline_1d( data, this%spline )
    call sll_s_interpolate_from_interpolant_array( coordinates, output_array, num_pts, &
         this%spline )

  end subroutine spline_interpolate1d

  subroutine spline_interpolate1d_disp(this, num_pts, data, alpha, output_array)
    class(sll_t_cubic_spline_interpolator_1d),  intent(in)       :: this
    !class(sll_t_cubic_spline_1d),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_pts
    sll_real64,  intent(in)   :: alpha
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_pts), intent(out)      :: output_array
    ! local variables
    sll_real64, dimension(num_pts)      :: coordinates
    sll_real64 :: length, delta
    sll_real64 :: xmin, xmax
    sll_int32 :: i
    ! compute the interpolating spline coefficients
    call sll_s_compute_cubic_spline_1d( data, this%spline )
    ! compute array of coordinates where interpolation is performed from displacement
    length = this%interpolation_points(this%num_points) - &
             this%interpolation_points(1)
    delta = this%interpolation_points(2) - this%interpolation_points(1)
    xmin = this%interpolation_points(1)
    xmax = this%interpolation_points(this%num_points)
    if (this%bc_type == sll_p_periodic) then
       ! The case alpha = 0.0 is problematic. We need to further try to make
       ! this computation in general m re efficient, minimize the use of modulo
       ! and even explore a uniform grid representation...
       if( alpha == 0.0_f64 ) then
          coordinates(:) = this%interpolation_points(:)
       else ! alpha != 0.0
          do i = 1, num_pts
             coordinates(i) = xmin + &
                  modulo(this%interpolation_points(i) - xmin + alpha, length)
!!$             write (*,'(a,z,f21.16,a,i,a,z,f21.16)') 'xmin = ', &
!!$                  xmin, xmin, '  coordinates(',i,') = ', coordinates(i), &
!!$                  coordinates(i)
             SLL_ASSERT(coordinates(i) >= xmin)
             SLL_ASSERT(coordinates(i) <= xmax)
          end do
       end if
    else ! any other BC? better a case statement
       if (alpha < 0 ) then
          do i = 1, num_pts
             coordinates(i) = max(this%interpolation_points(i) + alpha, xmin)
             SLL_ASSERT((xmin <=coordinates(i)).and.(coordinates(i) <= xmax))
          end do
       else
          do i = 1, num_pts
             coordinates(i) = min(this%interpolation_points(i) + alpha, xmax)
             SLL_ASSERT((xmin <=coordinates(i)).and.(coordinates(i) <= xmax))
          end do
       endif
    end if
    call sll_s_interpolate_from_interpolant_array( coordinates, output_array, num_pts, &
         this%spline )
  end subroutine spline_interpolate1d_disp

  subroutine spline_interpolate1d_disp_inplace(this, num_pts, data, alpha)
    class(sll_t_cubic_spline_interpolator_1d),  intent(in)       :: this
    !class(sll_t_cubic_spline_1d),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_pts
    sll_real64,  intent(in)   :: alpha
    sll_real64, dimension(num_pts), intent(inout)   :: data

    ! local variables
    sll_real64, dimension(num_pts)      :: coordinates
    sll_real64 :: length, delta
    sll_real64 :: xmin, xmax
    sll_int32 :: i
    ! compute the interpolating spline coefficients
    call sll_s_compute_cubic_spline_1d( data, this%spline )
    ! compute array of coordinates where interpolation is performed from displacement
    length = this%interpolation_points(this%num_points) - &
             this%interpolation_points(1)
    delta = this%interpolation_points(2) - this%interpolation_points(1)
    xmin = this%interpolation_points(1)
    xmax = this%interpolation_points(this%num_points)
    if (this%bc_type == sll_p_periodic) then
       ! The case alpha = 0.0 is problematic. We need to further try to make
       ! this computation in general m re efficient, minimize the use of modulo
       ! and even explore a uniform grid representation...
       if( alpha == 0.0_f64 ) then
          coordinates(:) = this%interpolation_points(:)
       else ! alpha != 0.0
          do i = 1, num_pts
             coordinates(i) = xmin + &
                  modulo(this%interpolation_points(i) - xmin + alpha, length)
!!$             write (*,'(a,z,f21.16,a,i,a,z,f21.16)') 'xmin = ', &
!!$                  xmin, xmin, '  coordinates(',i,') = ', coordinates(i), &
!!$                  coordinates(i)
             SLL_ASSERT(coordinates(i) >= xmin)
             SLL_ASSERT(coordinates(i) <= xmax)
          end do
       end if
    else ! any other BC? better a case statement
       if (alpha < 0 ) then
          do i = 1, num_pts
             coordinates(i) = max(this%interpolation_points(i) + alpha, xmin)
             SLL_ASSERT((xmin <=coordinates(i)).and.(coordinates(i) <= xmax))
          end do
       else
          do i = 1, num_pts
             coordinates(i) = min(this%interpolation_points(i) + alpha, xmax)
             SLL_ASSERT((xmin <=coordinates(i)).and.(coordinates(i) <= xmax))
          end do
       endif
    end if
    call sll_s_interpolate_from_interpolant_array( coordinates, data, num_pts, &
         this%spline )
  end subroutine spline_interpolate1d_disp_inplace


  subroutine compute_interpolants_cs1d( interpolator, data_array,&
       eta_coords, &
       size_eta_coords)
    class(sll_t_cubic_spline_interpolator_1d), intent(inout) :: interpolator
    sll_real64, dimension(:), intent(in)               :: data_array
    sll_real64, dimension(:), intent(in),optional  :: eta_coords
    sll_int32, intent(in),optional                 :: size_eta_coords
    call sll_s_compute_cubic_spline_1d( data_array, interpolator%spline )

    if(present(eta_coords))then
      !print *,'#warning eta_coords not taken into account'
    endif
    if(present(size_eta_coords))then
      !print *,'#warning size_eta_coords not taken into account'
    endif
  end subroutine compute_interpolants_cs1d

  ! Alternative implementation for the function meant to interpolate a
  ! whole array. This implementation fixes some problems in the previous
  ! function. Furthermore, it separates the operation into the more
  ! elementary steps: one is supposed to first compute the interpolants,
  ! then request to interpolate array values.
  subroutine interpolate_values_cs1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )
    class(sll_t_cubic_spline_interpolator_1d),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(num_pts), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(num_pts), intent(out)  :: output_array
    call sll_s_interpolate_from_interpolant_array( vals_to_interpolate, output_array, &
         num_pts, interpolator%spline )
  end subroutine interpolate_values_cs1d


  subroutine interpolate_derivatives_cs1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )
    class(sll_t_cubic_spline_interpolator_1d),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    call sll_s_interpolate_from_interpolant_derivatives_eta1( vals_to_interpolate, output_array, &
         num_pts,  interpolator%spline )
  end subroutine interpolate_derivatives_cs1d


  function interpolate_value_cs1d( interpolator, eta1 ) result(val)
    class(sll_t_cubic_spline_interpolator_1d), intent(in) :: interpolator
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    val = sll_f_interpolate_from_interpolant_value( eta1, interpolator%spline )
  end function

  function interpolate_deriv1_cs1d( interpolator, eta1 ) result(val)
    class(sll_t_cubic_spline_interpolator_1d), intent(in) :: interpolator
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
    val = sll_f_interpolate_derivative(eta1,interpolator%spline)
  end function interpolate_deriv1_cs1d

!PN DEFINED BUT NOT USED
! function interpolate_derivative_f95( interpolator, eta1 ) result(val)
!   class(sll_t_cubic_spline_interpolator_1d), intent(in) :: interpolator
!   sll_real64 :: val
!   sll_real64, intent(in) :: eta1
!   val = sll_f_interpolate_derivative(eta1,interpolator%spline)
! end function interpolate_derivative_f95

    !> PLEASE ADD DOCUMENTATION
  function sll_f_new_cubic_spline_interpolator_1d( &
    num_points, &
    xmin, &
    xmax, &
    bc_type, &
    slope_left, &
    slope_right ) result(res)

    type(sll_t_cubic_spline_interpolator_1d),  pointer :: res
    sll_int32,  intent(in)               :: num_points
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type
    sll_real64, intent(in), optional     :: slope_left
    sll_real64, intent(in), optional     :: slope_right
    sll_int32 :: ierr
    SLL_ALLOCATE(res,ierr)
    call initialize_cs1d_interpolator( &
         res, &
         num_points, &
         xmin, &
         xmax, &
         bc_type, &
         slope_left, &
         slope_right )
  end function sll_f_new_cubic_spline_interpolator_1d

  ! Why is the name of this function changing depending on the standard?
  ! only one will be compiled anyway!!

  !> initialize cubic spline interpolator
  subroutine initialize_cs1d_interpolator( &
    interpolator, &
    num_points, &
    xmin, &
    xmax, &
    bc_type, &
    slope_left, &
    slope_right )

    class(sll_t_cubic_spline_interpolator_1d),  intent(inout) :: interpolator
    sll_int32,  intent(in)               :: num_points
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type
    sll_real64, intent(in), optional     :: slope_left
    sll_real64, intent(in), optional     :: slope_right
    sll_int32                            :: ierr
    sll_int32  :: i
    sll_real64 :: delta

    interpolator%num_points = num_points
    SLL_ALLOCATE(interpolator%interpolation_points(num_points),ierr)
    interpolator%interpolation_points(1) = xmin
    delta = (xmax - xmin) / (num_points - 1)
    do i = 2, num_points
       interpolator%interpolation_points(i) = &
            interpolator%interpolation_points(i-1) + delta
    end do
    interpolator%interpolation_points(num_points) = xmax
    interpolator%bc_type = bc_type
    if (present(slope_left).and.present(slope_right)) then
       interpolator%spline => sll_f_new_cubic_spline_1d( &
            num_points, &
            xmin, xmax, &
            bc_type, &
            slope_left, &
            slope_right )
    else
       interpolator%spline => &
            sll_f_new_cubic_spline_1d(num_points, xmin, xmax, bc_type)
    end if

  end subroutine

  subroutine delete_cs1d( obj )
    class(sll_t_cubic_spline_interpolator_1d) :: obj
    call sll_o_delete(obj%spline)
  end subroutine delete_cs1d

  subroutine set_coefficients_cs1d( interpolator, coeffs )
    class(sll_t_cubic_spline_interpolator_1d),  intent(inout) :: interpolator
    sll_real64, dimension(:), intent(in), optional :: coeffs
    print *, '#set_coefficients_cs1d(): ERROR: This function has not been ', &
         'implemented yet.'
    if(present(coeffs))then
      print *,'#coefs are present'
    endif
    print *,interpolator%num_points
    stop
  end subroutine set_coefficients_cs1d


  function get_coefficients_cs1d(interpolator)
    class(sll_t_cubic_spline_interpolator_1d),  intent(in) :: interpolator
    sll_real64, dimension(:), pointer            :: get_coefficients_cs1d

    print *, 'get_coefficients_cs1d(): ERROR: This function has not been ', &
         'implemented yet.'
    get_coefficients_cs1d => null()
    print *,  interpolator%num_points
    stop
  end function get_coefficients_cs1d

end module sll_m_cubic_spline_interpolator_1d
