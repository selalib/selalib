!> @ingroup interpolators
!> @brief
!> Module for 1D interpolation and reconstruction
!> @details
!> This is an abstract class, methods are implemented in other modules
!> @todo
!> delete function for this type
module sll_m_interpolators_1d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_c_interpolator_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !>Abstract class for 1D interpolation and reconstruction
  type, abstract :: sll_c_interpolator_1d

contains

  ! Procedure with precomputation of interpolant(s)
  !> Compute coefficients of the interpolants and stores it in the interpolation object. 
   procedure(interpolator_1d_interpolant), deferred :: &
          compute_interpolants
   !> Set value of coefficients of the interpolant to the given values.
   procedure(interpolator_1d_set_coeffs), deferred :: &
        set_coefficients
   !> Extract the value of the precomputed coefficients of the interpolation from the interpolator object.
   procedure(interpolator_1d_get_coeffs), deferred :: &
        get_coefficients
   !> Compute the value of the interpolant of at a given abscissa \a x from the precomputed coefficients of the interpolator (stored in the interpolation object) 
   procedure(interpolator_one_arg_sub), deferred :: &
          interpolate_from_interpolant_value
   !>  Compute the value of the derivative of the interpolant of at a given abscissa \a x from the precomputed coefficients of the interpolator (stored in the interpolation object).
   procedure(interpolator_one_arg_sub), deferred :: &
          interpolate_from_interpolant_derivative_eta1
   !> Compute the value of the interpolant at several given abscissae given as an array from the precomputed coefficients of the interpolator.
   procedure(interpolator_1d_array_interpolant), deferred :: &
        interpolate_from_interpolant_array


   ! Procedures including whole interpolation process
   !> Compute the value of the interpolant at several given abscissae given as an array from given function values. Does not use a precomputed interpolant.
   procedure(interpolator_1d_array), deferred :: &
        interpolate_array
   !> Compute the value of the interpolant at all grid points shifted by the given displacement from function values. Does not use a precomputed interpolant.
   procedure(interpolator_1d_array_disp), deferred :: &
        interpolate_array_disp
   !> Compute the value of the interpolant at all grid points shifted by the given displacement. Does not use a precomputed interpolant.
   procedure(interpolator_1d_array_disp_inplace), deferred :: &
        interpolate_array_disp_inplace


end type sll_c_interpolator_1d




  !> Signature of compute_interpolants
  abstract interface
     subroutine interpolator_1d_interpolant( &
          interpolator, data_array,&
          eta_coords, &
          size_eta_coords)
       use sll_m_working_precision
       import :: sll_c_interpolator_1d
       class(sll_c_interpolator_1d), intent(inout)        :: interpolator !< interpolator object
       sll_real64,                      intent(in)           :: data_array(:)   !< data at the grid points
       sll_int32,                       intent(in),optional  :: size_eta_coords !< number of coordinates given
       sll_real64,                      intent(in),optional  :: eta_coords(:)   !< coordinates for which to compute the intepolant (optional argument for local interpolants)
     end subroutine interpolator_1d_interpolant
  end interface


  !> Signature of get_coefficients
  abstract interface
     function interpolator_1d_get_coeffs(interpolator)
       use sll_m_working_precision
       import :: sll_c_interpolator_1d
       class(sll_c_interpolator_1d), intent(in) :: interpolator !< interpolator object
       sll_real64, dimension(:), pointer           :: interpolator_1d_get_coeffs !< coefficients
     end function interpolator_1d_get_coeffs
  end interface

  !> Signature of set_coefficients
  abstract interface
     subroutine interpolator_1d_set_coeffs( interpolator, coeffs )
       use sll_m_working_precision
       import :: sll_c_interpolator_1d
       class(sll_c_interpolator_1d), intent(inout) :: interpolator  !< interpolator
       ! We allow the coefficients to be passed as 1d or 2d arrays. This allows
       ! for more flexibility for the children classes.
       sll_real64, dimension(:), intent(in), optional   :: coeffs !< coefficients to be set
     end subroutine interpolator_1d_set_coeffs
  end interface


  !> Signature of interpolate_from_interpolant_value and interpolate_from_interpolant_derivative_eta1
  abstract interface
     function interpolator_one_arg_sub( interpolator, eta1 ) result(val)
       use sll_m_working_precision
       import :: sll_c_interpolator_1d
       sll_real64                                     :: val  !< interpolated value
       class(sll_c_interpolator_1d     ), intent(in)    :: interpolator !< interpolator object
       sll_real64,                      intent(in)    :: eta1 !< abscissa where to interpolate
     end function interpolator_one_arg_sub
  end interface



  !> Signature of interpolate_from_interpolant_value
  abstract interface
     subroutine interpolator_1d_array_interpolant( &
       interpolator, &
       num_pts, &
       vals_to_interpolate, &
       output_array )

       use sll_m_working_precision
       import :: sll_c_interpolator_1d
       class(sll_c_interpolator_1d), intent(in) :: interpolator !< interpolator object
       sll_int32,                       intent(in) :: num_pts      !< size of output array
       sll_real64,                      intent(in) :: vals_to_interpolate(num_pts) !< abscissae where to interpolate (size num_pts)
       sll_real64,                      intent(out):: output_array(num_pts) !< interpolated values at \a vals_to_interpolate
     end subroutine interpolator_1d_array_interpolant
  end interface



  !> Signature of interpolate_array
  abstract interface
     subroutine interpolator_1d_array(this, num_pts, data, coordinates, output_array) 
       use sll_m_working_precision
       import :: sll_c_interpolator_1d
       class(sll_c_interpolator_1d), intent(in)     :: this !< interpolator object
       sll_int32,                       intent(in)     :: num_pts    !< size of output array
       sll_real64,                      intent(in)     :: data(:)  !< function values at grid points
       sll_real64,                      intent(in)     :: coordinates(num_pts) !<  points where output is desired (size num_pts)
       sll_real64,                      intent(out)    :: output_array(num_pts) !< interpolated values at \a coordinates
     end subroutine interpolator_1d_array
  end interface



  !> Signature of interpolate_array_disp
  abstract interface
     subroutine interpolator_1d_array_disp( &
       this, &
       num_pts, &
       data, &
       alpha, &
       output_array)

       use sll_m_working_precision
       import :: sll_c_interpolator_1d
       class(sll_c_interpolator_1d), intent(in)     :: this !< interpolator object
       sll_int32,                       intent(in)     :: num_pts    !< size of output array
       sll_real64,                      intent(in)     :: data(:)  !< data to be interpolated
       sll_real64,                      intent(in)     :: alpha !< displacement
       sll_real64,                      intent(out)    :: output_array(num_pts) !< interpolated values

     end subroutine interpolator_1d_array_disp
  end interface

  !> Signature of interpolate_array_disp
  abstract interface
     subroutine interpolator_1d_array_disp_inplace( &
       this, &
       num_pts, &
       data, &
       alpha)

       use sll_m_working_precision
       import :: sll_c_interpolator_1d
       class(sll_c_interpolator_1d), intent(in)     :: this !< interpolator object
       sll_int32,                       intent(in)     :: num_pts    !< size of output array
       sll_real64,                      intent(inout)  :: data(num_pts)  !< data to be interpolated
       sll_real64,                      intent(in)     :: alpha !< displacement

     end subroutine interpolator_1d_array_disp_inplace
  end interface




end module sll_m_interpolators_1d_base
