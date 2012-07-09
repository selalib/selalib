module sll_cubic_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_module_interpolators_1d_base
use sll_splines
  implicit none
  
  type, extends(sll_interpolator_1d_base) ::  cubic_spline_1d_interpolator
     sll_int32                     :: num_points ! size
     sll_int32                     :: bc_type
     type(sll_spline_1D), pointer  :: spline
   contains
     procedure, pass(interpolator) :: initialize => initialize_cs1d_interpolator
     procedure :: compute_interpolants => compute_interpolants_cs1d
     procedure :: interpolate_value => interpolate_value_cs1d
     procedure :: interpolate_derivative_eta1 => interpolate_deriv1_cs1d
     procedure :: interpolate_array_values => interpolate_values_cs1d
     procedure :: interpolate_pointer_values => interpolate_pointer_values_cs1d
     procedure :: interpolate_array_derivatives => interpolate_derivatives_cs1d
     procedure :: interpolate_pointer_derivatives => &
          interpolate_pointer_derivatives_cs1d
     procedure, pass:: interpolate_array => spline_interpolate1d
 !    procedure, pass:: reconstruct_array
     !generic :: initialize => initialize_cs1d_interpolator
  end type cubic_spline_1d_interpolator


  
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
  function spline_interpolate1d(this, num_points, data, coordinates) &
       result(data_out)
    class(cubic_spline_1d_interpolator),  intent(in)       :: this
    !class(sll_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    ! local variables
    sll_int32 :: ierr
    ! compute the interpolating spline coefficients
    call compute_spline_1D( data, this%bc_type, this%spline )
    call interpolate_array_values( coordinates, data_out, num_points, &
         this%spline )
  end function spline_interpolate1d
  
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

    class(cubic_spline_1d_interpolator),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    sll_int32 :: ierr
    call interpolate_array_values( vals_to_interpolate, output_array, &
         num_pts, interpolator%spline )
  end subroutine interpolate_values_cs1d

  subroutine interpolate_pointer_values_cs1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output )

    class(cubic_spline_1d_interpolator),  intent(in) :: interpolator
    sll_int32,  intent(in)            :: num_pts
    sll_real64, dimension(:), pointer :: vals_to_interpolate
    sll_real64, dimension(:), pointer :: output
    sll_int32 :: ierr
    call interpolate_pointer_values( vals_to_interpolate, output, &
         num_pts, interpolator%spline )
  end subroutine interpolate_pointer_values_cs1d


  subroutine interpolate_derivatives_cs1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )

    class(cubic_spline_1d_interpolator),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    sll_int32 :: ierr
    call interpolate_array_derivatives( vals_to_interpolate, num_pts, &
         output_array, interpolator%spline )
  end subroutine interpolate_derivatives_cs1d

  subroutine interpolate_pointer_derivatives_cs1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output )

    class(cubic_spline_1d_interpolator),  intent(in) :: interpolator
    sll_int32,  intent(in)              :: num_pts
    sll_real64, dimension(:), pointer   :: vals_to_interpolate
    sll_real64, dimension(:), pointer   :: output
    sll_int32 :: ierr
    call interpolate_pointer_derivatives( vals_to_interpolate, num_pts, &
         output, interpolator%spline )
  end subroutine interpolate_pointer_derivatives_cs1d


  subroutine compute_interpolants_cs1d( interpolator, data_array )
    class(cubic_spline_1d_interpolator), intent(inout) :: interpolator
    sll_real64, dimension(:), intent(in)               :: data_array
    call compute_spline_1D_bis( data_array, interpolator%spline )
  end subroutine compute_interpolants_cs1d

  function interpolate_value_cs1d( interpolator, eta1 ) result(val)
    sll_real64 :: val
    class(cubic_spline_1d_interpolator), intent(inout) :: interpolator
    sll_real64, intent(in) :: eta1
    val = interpolate_value_1D( eta1, interpolator%spline )
  end function interpolate_value_cs1d
  
  function interpolate_deriv1_cs1d( interpolator, eta1 ) result(val)
    sll_real64 :: val
    class(cubic_spline_1d_interpolator), intent(inout) :: interpolator
    sll_real64, intent(in) :: eta1
    val = interpolate_derivative(eta1,interpolator%spline)
  end function interpolate_deriv1_cs1d

  !> create new spline-based interpolator
  subroutine initialize_cs1d_interpolator( &
    interpolator, &
    num_points, &
    xmin, &
    xmax, &
    bc_type, &
    slope_left, &
    slope_right )

    class(cubic_spline_1d_interpolator),  intent(inout) :: interpolator 
    sll_int32,  intent(in)               :: num_points
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type
    sll_real64, intent(in), optional     :: slope_left
    sll_real64, intent(in), optional     :: slope_right
    sll_int32                            :: ierr
    
    interpolator%num_points = num_points
    interpolator%bc_type = bc_type
    if (present(slope_left).and.present(slope_right)) then
       interpolator%spline => new_spline_1D(num_points, xmin, xmax, bc_type, slope_left, slope_right )
    else
       interpolator%spline => new_spline_1D(num_points, xmin, xmax, bc_type)
    end if
  end subroutine initialize_cs1d_interpolator

  function reconstruct_array(this, num_points, data) result(res)
    ! dummy procedure
    class(cubic_spline_1d_interpolator), intent(in)     :: this
       sll_int32, intent(in)                :: num_points! size of output array
       sll_real64, dimension(:), intent(in) :: data   ! data to be interpolated 
       sll_real64, dimension(num_points)    :: res
       res(:) = 0.0_f64
  end function reconstruct_array
  
end module sll_cubic_spline_interpolator_1d
