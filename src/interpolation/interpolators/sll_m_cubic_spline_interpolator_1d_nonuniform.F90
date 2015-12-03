!> @ingroup interpolators
!> @brief
!> Implements sll_c_interpolator_1d with cubic splines on non uniform mesh
!> @details
!> Define spline interpolation of values in data define on original grid at
!> points coordinates
module sll_m_cubic_spline_interpolator_1d_nonuniform
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_m_interpolators_1d_base
use sll_m_cubic_non_uniform_splines
use sll_m_cubic_splines
implicit none
private

  !> sll_interpolator_1d implemented with cubic splines on non uniform mesh
  type, public, extends(sll_c_interpolator_1d) :: sll_cubic_spline_interpolator_1d_nonuniform
     sll_real64, dimension(:), pointer      :: interpolation_points !< points
     sll_int32                              :: num_points     !< size
     sll_int32                              :: bc_type        !< boundary condition
     type(sll_cubic_spline_1D), pointer     :: spline         !< cubic spline
     type(cubic_nonunif_spline_1D), pointer :: nonunif_spline !< spline
   contains
     !> PLEASE ADD DOCUMENTATION
     procedure, pass(interpolator) :: initialize => initialize_cs1d_interpolator2
     !> PLEASE ADD DOCUMENTATION
     procedure :: compute_interpolants => compute_interpolants_cs1d
     !> PLEASE ADD DOCUMENTATION
     procedure :: interpolate_from_interpolant_value => interpolate_value_cs1d
     !> PLEASE ADD DOCUMENTATION
     procedure :: interpolate_from_interpolant_derivative_eta1 => interpolate_deriv1_cs1d
     !> PLEASE ADD DOCUMENTATION
     procedure :: interpolate_from_interpolant_array => interpolate_values_cs1d
     !> PLEASE ADD DOCUMENTATION
     procedure :: interpolate_from_interpolant_derivatives_eta => interpolate_derivatives_cs1d
     !> PLEASE ADD DOCUMENTATION
     procedure, pass:: interpolate_array => spline_interpolate1d
     !> PLEASE ADD DOCUMENTATION
     procedure, pass:: interpolate_array_disp => spline_interpolate1d_disp
     !generic :: initialize => initialize_cs1d_interpolator
     !> PLEASE ADD DOCUMENTATION
     procedure, pass :: set_coefficients => set_coefficients_cs1d
     !> PLEASE ADD DOCUMENTATION
     procedure, pass :: get_coefficients => get_coefficients_cs1d
  end type sll_cubic_spline_interpolator_1d_nonuniform

  !> Deallocate the interpolator object
  interface sll_delete
     module procedure delete_cs1d
  end interface sll_delete

public sll_delete

contains  ! ****************************************************************



  subroutine spline_interpolate1d(this, num_pts, data, coordinates, output_array)
    class(sll_cubic_spline_interpolator_1d_nonuniform),  intent(in)       :: this
    !class(sll_cubic_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(num_pts), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_pts), intent(out)      :: output_array
    ! compute the interpolating spline coefficients
    call compute_cubic_spline_1D( data, this%spline )
    call interpolate_from_interpolant_array( coordinates, output_array, num_pts, &
         this%spline )
  end subroutine spline_interpolate1d


  subroutine spline_interpolate1d_disp(this, num_pts, data, alpha, output_array)
    class(sll_cubic_spline_interpolator_1d_nonuniform),  intent(in)       :: this
    !class(sll_cubic_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_pts
    sll_real64,  intent(in)   :: alpha
    sll_real64, dimension(num_pts), intent(inout)   :: data
    sll_real64, dimension(num_pts),intent(inout)      :: output_array

    sll_real64, dimension(num_pts)      :: coordinates
    sll_real64 :: length, delta
    sll_real64 :: xmin, xmax
    sll_int32 :: i
    ! compute the interpolating spline coefficients
    call compute_cubic_spline_1D( data, this%spline )
    ! compute array of coordinates where interpolation is performed from displacement
    length = this%interpolation_points(num_pts) - &
             this%interpolation_points(1)
    delta = this%interpolation_points(2) - this%interpolation_points(1)
    xmin = this%interpolation_points(1)
    xmax = this%interpolation_points(num_pts)
    if (this%bc_type == SLL_PERIODIC) then
       do i = 1, num_pts
          coordinates(i) = xmin + modulo(this%interpolation_points(i) - xmin - alpha, length)
          SLL_ASSERT(coordinates(i) >= xmin)
          SLL_ASSERT(coordinates(i) <= xmax)
       end do
    else
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
    call interpolate_from_interpolant_array( coordinates, output_array, num_pts, &
         this%spline )
  end subroutine spline_interpolate1d_disp




  ! Both versions F03 and F95 of compute_interpolants_cs1d should have the
  ! same name. In the F95 we should add a generic interface around this
  ! subroutine, selecting on the type of interpolator. In the F03 case the
  ! interface is the compute_interpolants routine which gets assigned to
  ! the cs1d at initialization time.
    subroutine compute_interpolants_cs1d(interpolator, data_array,&
         eta_coords, &
         size_eta_coords)
      class(sll_cubic_spline_interpolator_1d_nonuniform),intent(inout)::interpolator

      sll_real64, dimension(:), intent(in)           :: data_array
      sll_real64, dimension(:), intent(in),optional  :: eta_coords
      sll_int32, intent(in),optional                 :: size_eta_coords
      if(present(eta_coords))then
        !print *,'#Warning eta_coords present but not used'
      endif
      if(present(size_eta_coords))then
        !print *,'#Warning size_eta_coords present but not used'
      endif
      call compute_cubic_spline_1D( data_array, interpolator%spline )
  end subroutine

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
    class(sll_cubic_spline_interpolator_1d_nonuniform),  intent(in) :: interpolator
    sll_int32,  intent(in)                       :: num_pts
    sll_real64, dimension(num_pts), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(num_pts), intent(out)  :: output_array
    call interpolate_from_interpolant_array( vals_to_interpolate, output_array, &
         num_pts, interpolator%spline )
  end subroutine interpolate_values_cs1d


  subroutine interpolate_derivatives_cs1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )
    class(sll_cubic_spline_interpolator_1d_nonuniform),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    call interpolate_array_derivatives( vals_to_interpolate, output_array, &
         num_pts, interpolator%spline )
  end subroutine interpolate_derivatives_cs1d

  function interpolate_value_cs1d( interpolator, eta1 ) result(val)
    class(sll_cubic_spline_interpolator_1d_nonuniform), intent(in) :: interpolator
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    val = interpolate_from_interpolant_value( eta1, interpolator%spline )
  end function

  function interpolate_deriv1_cs1d( interpolator, eta1 ) result(val)
    class(sll_cubic_spline_interpolator_1d_nonuniform), intent(in) :: interpolator
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
    val = interpolate_derivative(eta1,interpolator%spline)
  end function

!PN DEFINED BUT NOT USED
!  function interpolate_derivative_f95( interpolator, eta1 ) result(val)
!    class(sll_cubic_spline_interpolator_1d_nonuniform), intent(in) :: interpolator
!    sll_real64 :: val
!    sll_real64, intent(in) :: eta1
!    val = interpolate_derivative(eta1,interpolator%spline)
!  end function

  ! Why is the name of this function changing depending on the standard?
  ! only one will be compiled anyway!!

  !> initialize cubic spline interpolator
  subroutine initialize_cs1d_interpolator2( &
    interpolator, &
    num_points, &
    xmin, &
    xmax, &
    bc_type, &
    slope_left, &
    slope_right )

    class(sll_cubic_spline_interpolator_1d_nonuniform),  intent(inout) :: interpolator
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

  subroutine delete_cs1d( obj )
    class(sll_cubic_spline_interpolator_1d_nonuniform) :: obj
    call sll_delete(obj%spline)
  end subroutine delete_cs1d


  subroutine set_coefficients_cs1d( interpolator, coeffs )
    class(sll_cubic_spline_interpolator_1d_nonuniform), intent(inout)  :: interpolator
    sll_real64, dimension(:), intent(in), optional :: coeffs
    print *, 'set_coefficients_cs1d(): ERROR: This function has not been ', &
         'implemented yet.'
    print *,interpolator%num_points
    if(present(coeffs))then
      print *,'coeffs are present'
    endif
    stop
  end subroutine set_coefficients_cs1d


  function get_coefficients_cs1d(interpolator)
    class(sll_cubic_spline_interpolator_1d_nonuniform), intent(in)  :: interpolator
    sll_real64, dimension(:), pointer            :: get_coefficients_cs1d

    print *, 'get_coefficients_cs1d(): ERROR: This function has not been ', &
         'implemented yet.'
    print *,interpolator%num_points
    get_coefficients_cs1d => null()
    stop
  end function get_coefficients_cs1d


end module sll_m_cubic_spline_interpolator_1d_nonuniform
