module sll_cubic_spline_interpolator_2d
#include "sll_working_precision.h"
  use sll_interpolators_base
  use sll_splines
  implicit none

  ! The spline-based interpolator is only a wrapper around the capabilities
  ! of the cubic splines. All interpolators share a common interface with
  ! respect to their use, as described by the interpolator_2d_base class.
  !
  ! Where the diverse interpolators diverge is in the way to initialize them.
  type, extends(interpolator_2d_base) :: cubic_spline_2d_interpolator
     sll_int32                    :: npts1
     sll_int32                    :: npts2
     type(sll_spline_2D), pointer :: spline
   contains
     procedure, pass(interpolator) :: initialize=>initialize_cs2d_interpolator
     procedure :: compute_interpolants => compute_interpolants_cs2d
     procedure :: interpolate_value => interpolate_value_cs2d
     procedure :: interpolate_derivative_eta1 => interpolate_deriv1_cs2d
     procedure :: interpolate_derivative_eta2 => interpolate_deriv2_cs2d
  end type cubic_spline_2d_interpolator

contains

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
    eta1_bc_type,   &
    eta2_bc_type,   &
    const_eta1_min_slope, &
    const_eta1_max_slope, &
    const_eta2_min_slope, &
    const_eta2_max_slope, &
    eta1_min_slopes, &
    eta1_max_slopes, &
    eta2_min_slopes, &
    eta2_max_slopes )

    class(cubic_spline_2d_interpolator), intent(inout) :: interpolator
    sll_int32, intent(in)                         :: npts1
    sll_int32, intent(in)                         :: npts2
    sll_real64, intent(in)                        :: eta1_min
    sll_real64, intent(in)                        :: eta1_max
    sll_real64, intent(in)                        :: eta2_min
    sll_real64, intent(in)                        :: eta2_max
    sll_int32, intent(in), optional               :: eta1_bc_type
    sll_int32, intent(in), optional               :: eta2_bc_type
    sll_real64, intent(in), optional              :: const_eta1_min_slope
    sll_real64, intent(in), optional              :: const_eta1_max_slope
    sll_real64, intent(in), optional              :: const_eta2_min_slope
    sll_real64, intent(in), optional              :: const_eta2_max_slope
    sll_real64, dimension(:),intent(in), optional :: eta1_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta1_max_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_max_slopes

    interpolator%npts1 = npts1
    interpolator%npts2 = npts2
    interpolator%spline => new_spline_2D( &
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
         x2_max_slopes=eta2_max_slopes )
  end subroutine initialize_cs2d_interpolator

  subroutine compute_interpolants_cs2d( interpolator, data_array )
    class(cubic_spline_2d_interpolator), intent(inout) :: interpolator
    sll_real64, dimension(:,:), intent(in) :: data_array
    call compute_spline_2D( data_array, interpolator%spline )
  end subroutine compute_interpolants_cs2d

  function interpolate_value_cs2d( interpolator, eta1, eta2 ) result(val)
    sll_real64 :: val
    class(cubic_spline_2d_interpolator), intent(inout) :: interpolator
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = interpolate_value_2D( eta1, eta2, interpolator%spline )
  end function interpolate_value_cs2d

  function interpolate_deriv1_cs2d( interpolator, eta1, eta2 ) result(val)
    sll_real64 :: val
    class(cubic_spline_2d_interpolator), intent(inout) :: interpolator
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = interpolate_x1_derivative_2D(eta1,eta2,interpolator%spline)
  end function interpolate_deriv1_cs2d

  function interpolate_deriv2_cs2d( interpolator, eta1, eta2 ) result(val)
    sll_real64 :: val
    class(cubic_spline_2d_interpolator), intent(inout) :: interpolator
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = interpolate_x2_derivative_2D(eta1,eta2,interpolator%spline)
  end function interpolate_deriv2_cs2d


end module sll_cubic_spline_interpolator_2d
