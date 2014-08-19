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

module sll_cubic_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#ifndef STDF95
use sll_module_interpolators_1d_base
#endif
use sll_cubic_splines
  implicit none

!> Object for 1d cubic spline interpolation on uniform mesh
#ifdef STDF95

type ::  cubic_spline_1d_interpolator
   sll_real64, dimension(:), pointer :: interpolation_points !< points positions
   sll_int32                         :: num_points           !< size
   sll_int32                         :: bc_type            !< boundary condition
   type(sll_cubic_spline_1D), pointer      :: spline       !< spline object
end type cubic_spline_1d_interpolator

#else

type, extends(sll_interpolator_1d_base) ::  cubic_spline_1d_interpolator

   sll_real64, dimension(:), pointer :: interpolation_points !< points position
   sll_int32                         :: num_points           !< size
   sll_int32                         :: bc_type            !< boundary condition
   type(sll_cubic_spline_1D), pointer      :: spline       !< spline object

contains

procedure, pass(interpolator) :: initialize => initialize_cs1d_interpolator
procedure, pass :: compute_interpolants => compute_interpolants_cs1d
procedure, pass :: interpolate_value => interpolate_value_cs1d
procedure, pass :: interpolate_derivative_eta1 => interpolate_deriv1_cs1d
procedure, pass :: interpolate_array_values => interpolate_values_cs1d
procedure, pass :: interpolate_pointer_values => interpolate_pointer_values_cs1d
procedure, pass :: interpolate_array_derivatives => interpolate_derivatives_cs1d
procedure, pass :: interpolate_pointer_derivatives => interpolate_pointer_derivatives_cs1d
procedure, pass:: interpolate_array => spline_interpolate1d
procedure, pass:: interpolate_array_disp => spline_interpolate1d_disp
procedure, pass:: reconstruct_array => reconstruct_array ! this is suspicious...
procedure, pass :: set_coefficients => set_coefficients_cs1d
procedure, pass :: get_coefficients => get_coefficients_cs1d
!generic :: initialize => initialize_cs1d_interpolator

end type cubic_spline_1d_interpolator


type :: cubic_spline_1d_interpolator_ptr
   type(cubic_spline_1d_interpolator), pointer :: interp
end type cubic_spline_1d_interpolator_ptr


#endif

  interface delete
     module procedure delete_cs1d
  end interface delete

contains  ! ****************************************************************


  ! the following provides an implementation for the abstract interface
  !interpolate1d
  !> Define spline interpolation of values in data define on original grid at
  !> points coordinates
  ! Issues with the following function:
  ! - entities referenced through "this" are modified, violating the declared
  !   intent.
  ! - it is probably better to convert this into a subroutine, since data_out
  !   will be allocated on the stack (too big an array will crash the program),
  !   and some copy operation might be involved when "catching" the results.

#ifdef STDF95
  function cubic_spline_interpolate_array(this, num_points, data, coordinates) &
       result(data_out)
    type(cubic_spline_1d_interpolator),  intent(in)       :: this
#else
  function spline_interpolate1d(this, num_points, data, coordinates) &
       result(data_out)
    class(cubic_spline_1d_interpolator),  intent(in)       :: this
#endif
    !class(sll_cubic_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    ! compute the interpolating spline coefficients
    call compute_cubic_spline_1D( data, this%spline )
    call interpolate_array_values( coordinates, data_out, num_points, &
         this%spline )
  end function

#ifdef STDF95
  function cubic_spline_interpolate_array_at_displacement(this, num_points, &
       data, alpha ) & ! coordinates) &
       result(data_out)
    type(cubic_spline_1d_interpolator),  intent(in)       :: this
#else
  function spline_interpolate1d_disp(this, num_points, data, alpha) &
       result(data_out)
    class(cubic_spline_1d_interpolator),  intent(in)       :: this
#endif
    !class(sll_cubic_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64,  intent(in)   :: alpha
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    ! local variables
    sll_real64, dimension(num_points)      :: coordinates
    sll_real64 :: length, delta
    sll_real64 :: xmin, xmax
    sll_int32 :: i
    ! compute the interpolating spline coefficients
    call compute_cubic_spline_1D( data, this%spline )
    ! compute array of coordinates where interpolation is performed from displacement
    length = this%interpolation_points(num_points) - &
             this%interpolation_points(1)
    delta = this%interpolation_points(2) - this%interpolation_points(1)
    xmin = this%interpolation_points(1)
    xmax = this%interpolation_points(num_points)
    if (this%bc_type == SLL_PERIODIC) then
       ! The case alpha = 0.0 is problematic. We need to further try to make
       ! this computation in general m re efficient, minimize the use of modulo
       ! and even explore a uniform grid representation...
       if( alpha == 0.0_f64 ) then
          coordinates(:) = this%interpolation_points(:)
       else ! alpha != 0.0
          do i = 1, num_points
             coordinates(i) = xmin + &
                  modulo(this%interpolation_points(i) - xmin - alpha, length)
!!$             write (*,'(a,z,f21.16,a,i,a,z,f21.16)') 'xmin = ', &
!!$                  xmin, xmin, '  coordinates(',i,') = ', coordinates(i), &
!!$                  coordinates(i)
             SLL_ASSERT(coordinates(i) >= xmin)
             SLL_ASSERT(coordinates(i) <= xmax)
          end do
       end if
    else ! any other BC? better a case statement
       if (alpha > 0 ) then
          do i = 1, num_points
             coordinates(i) = max(this%interpolation_points(i) - alpha, xmin)
             SLL_ASSERT((xmin <=coordinates(i)).and.(coordinates(i) <= xmax))
          end do
       else
          do i = 1, num_points
             coordinates(i) = min(this%interpolation_points(i) - alpha, xmax)
             SLL_ASSERT((xmin <=coordinates(i)).and.(coordinates(i) <= xmax))
          end do
       endif
    end if
    call interpolate_array_values( coordinates, data_out, num_points, &
         this%spline )
  end function

  ! Both versions F03 and F95 of compute_interpolants_cs1d should have the
  ! same name. In the F95 we should add a generic interface around this
  ! subroutine, selecting on the type of interpolator. In the F03 case the
  ! interface is the compute_interpolants routine which gets assigned to
  ! the cs1d at initialization time.
#ifdef STDF95
  subroutine cubic_spline_compute_interpolants( interpolator, data_array,&
       eta_coords, &
       size_eta_coords)
    type(cubic_spline_1d_interpolator), intent(inout)  :: interpolator
#else
  subroutine compute_interpolants_cs1d( interpolator, data_array,&
       eta_coords, &
       size_eta_coords)
    class(cubic_spline_1d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64, dimension(:), intent(in)               :: data_array
    sll_real64, dimension(:), intent(in),optional  :: eta_coords
    sll_int32, intent(in),optional                 :: size_eta_coords
    call compute_cubic_spline_1D( data_array, interpolator%spline )

    if(present(eta_coords))then
      !print *,'#warning eta_coords not taken into account'
    endif
    if(present(size_eta_coords))then
      !print *,'#warning size_eta_coords not taken into account'
    endif
#ifdef STDF95
  end subroutine cubic_spline_compute_interpolants
#else
  end subroutine compute_interpolants_cs1d
#endif

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
#ifdef STDF95
    type(cubic_spline_1d_interpolator),  intent(in) :: interpolator
#else
    class(cubic_spline_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    call interpolate_array_values( vals_to_interpolate, output_array, &
         num_pts, interpolator%spline )
  end subroutine interpolate_values_cs1d

  subroutine interpolate_pointer_values_cs1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output )
#ifdef STDF95
    type(cubic_spline_1d_interpolator),  intent(in) :: interpolator
#else
    class(cubic_spline_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)            :: num_pts
    sll_real64, dimension(:), pointer :: vals_to_interpolate
    sll_real64, dimension(:), pointer :: output
    call interpolate_pointer_values( vals_to_interpolate, output, &
         num_pts, interpolator%spline )
  end subroutine interpolate_pointer_values_cs1d


  subroutine interpolate_derivatives_cs1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )
#ifdef STDF95
    type(cubic_spline_1d_interpolator),  intent(in) :: interpolator
#else
    class(cubic_spline_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    call interpolate_array_derivatives( vals_to_interpolate, output_array, &
         num_pts,  interpolator%spline )
  end subroutine interpolate_derivatives_cs1d

  subroutine interpolate_pointer_derivatives_cs1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output )
#ifdef STDF95
    type(cubic_spline_1d_interpolator),  intent(in) :: interpolator
#else
    class(cubic_spline_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)              :: num_pts
    sll_real64, dimension(:), pointer   :: vals_to_interpolate
    sll_real64, dimension(:), pointer   :: output
    call interpolate_pointer_derivatives( vals_to_interpolate, output, &
         num_pts, interpolator%spline )
  end subroutine interpolate_pointer_derivatives_cs1d

#ifdef STDF95
  function cubic_spline_interpolate_value( interpolator, eta1 ) result(val)
    type(cubic_spline_1d_interpolator), intent(in) :: interpolator
#else
  function interpolate_value_cs1d( interpolator, eta1 ) result(val)
    class(cubic_spline_1d_interpolator), intent(in) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    val = interpolate_value( eta1, interpolator%spline )
  end function

#ifdef STDF95
  function cubic_spline_interpolate_derivative_eta1( interpolator, eta1 ) &
       result(val)
    type(cubic_spline_1d_interpolator), intent(in)  :: interpolator
#else
  function interpolate_deriv1_cs1d( interpolator, eta1 ) result(val)
    class(cubic_spline_1d_interpolator), intent(in) :: interpolator
#endif
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
    val = interpolate_derivative(eta1,interpolator%spline)
#ifdef STDF95
  end function cubic_spline_interpolate_derivative_eta1
#else
  end function interpolate_deriv1_cs1d
#endif

#ifdef STDF95
  function cubic_spline_interpolate_derivative_f95( interpolator, eta1 ) result(val)
    type(cubic_spline_1d_interpolator), intent(in) :: interpolator
#else
  function interpolate_derivative_f95( interpolator, eta1 ) result(val)
    class(cubic_spline_1d_interpolator), intent(in) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    val = interpolate_derivative(eta1,interpolator%spline)
#ifdef STDF95
  end function cubic_spline_interpolate_derivative_f95
#else
  end function interpolate_derivative_f95
#endif


  function new_cubic_spline_1d_interpolator( &
    num_points, &
    xmin, &
    xmax, &
    bc_type, &
    slope_left, &
    slope_right ) result(res)

    type(cubic_spline_1d_interpolator),  pointer :: res
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
  end function new_cubic_spline_1d_interpolator

  ! Why is the name of this function changing depending on the standard?
  ! only one will be compiled anyway!!

  !> initialize cubic spline interpolator
#ifdef STDF95
  subroutine cubic_spline_1d_interpolator_initialize( &
#else
  subroutine initialize_cs1d_interpolator( &
#endif
    interpolator, &
    num_points, &
    xmin, &
    xmax, &
    bc_type, &
    slope_left, &
    slope_right )

#ifdef STDF95
    type(cubic_spline_1d_interpolator),  intent(inout)  :: interpolator
#else
    class(cubic_spline_1d_interpolator),  intent(inout) :: interpolator
#endif
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
       interpolator%spline => new_cubic_spline_1D( &
            num_points, &
            xmin, xmax, &
            bc_type, &
            slope_left, &
            slope_right )
    else
       interpolator%spline => &
            new_cubic_spline_1D(num_points, xmin, xmax, bc_type)
    end if
  end subroutine

  function reconstruct_array(this, num_points, data) result(res)
    ! dummy procedure
#ifdef STDF95
    type(cubic_spline_1d_interpolator), intent(in)      :: this
#else
    class(cubic_spline_1d_interpolator), intent(in)     :: this
#endif
       sll_int32, intent(in)                :: num_points! size of output array
       sll_real64, dimension(:), intent(in) :: data   ! data to be interpolated
       sll_real64, dimension(num_points)    :: res
       res(:) = 0.0_f64
       print *,'#warning reconstruct_array is dummy'
       print *,'#', this%num_points
       print *,maxval(data)

  end function reconstruct_array

  subroutine delete_cs1d( obj )
#ifdef STDF95
    type(cubic_spline_1d_interpolator) :: obj
#else
    class(cubic_spline_1d_interpolator) :: obj
#endif
    call sll_delete(obj%spline)
  end subroutine delete_cs1d

  subroutine set_coefficients_cs1d( interpolator, coeffs )
#ifdef STDF95
    type(cubic_spline_1d_interpolator),  intent(inout)  :: interpolator
#else
    class(cubic_spline_1d_interpolator),  intent(inout) :: interpolator
#endif
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
#ifdef STDF95
    type(cubic_spline_1d_interpolator),  intent(in)  :: interpolator
#else
    class(cubic_spline_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_real64, dimension(:), pointer            :: get_coefficients_cs1d

    print *, 'get_coefficients_cs1d(): ERROR: This function has not been ', &
         'implemented yet.'
    get_coefficients_cs1d => null()
    print *,  interpolator%num_points
    stop
  end function get_coefficients_cs1d

end module sll_cubic_spline_interpolator_1d
