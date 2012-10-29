module sll_cubic_spline_interpolator_2d
#include "sll_working_precision.h"
#ifndef STDF95
  use sll_module_interpolators_2d_base
#endif
  use sll_splines
  implicit none

  ! The spline-based interpolator is only a wrapper around the capabilities
  ! of the cubic splines. All interpolators share a common interface with
  ! respect to their use, as described by the interpolator_2d_base class.
  !
  ! Where the diverse interpolators diverge is in the way to initialize them.
#ifdef STDF95
  type                                :: cubic_spline_2d_interpolator
#else
  type, extends(sll_interpolator_2d_base) :: cubic_spline_2d_interpolator
#endif
     sll_int32                    :: npts1
     sll_int32                    :: npts2
     type(sll_spline_2D), pointer :: spline
#ifdef STDF95
#else
   contains
     procedure, pass(interpolator) :: initialize=>initialize_cs2d_interpolator
     procedure :: compute_interpolants => compute_interpolants_cs2d
     procedure :: interpolate_value => interpolate_value_cs2d
     procedure :: interpolate_derivative_eta1 => interpolate_deriv1_cs2d
     procedure :: interpolate_derivative_eta2 => interpolate_deriv2_cs2d
     procedure, pass:: interpolate_array => spline_interpolate2d
     procedure, pass:: interpolate_array_disp => spline_interpolate2d_disp
#endif
  end type cubic_spline_2d_interpolator

contains

  ! We allow to use the enumerators of the splines module in this interpolator
  ! because:
  ! a. This is just a wrapper and is intimately related to the underlying
  !    cubic splines module.
  ! b. There is no uniform interface for the initialization anyway.
  ! The underlying implementation with the splines module could be hidden but
  ! I can't see a compelling reason why.
#ifdef STDF95
  subroutine cubic_spline_initialize( &
#else
  subroutine initialize_cs2d_interpolator( &
#endif
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

#ifdef STDF95
    type(cubic_spline_2d_interpolator), intent(inout) :: interpolator
#else
    class(cubic_spline_2d_interpolator), intent(inout) :: interpolator
#endif
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
  end subroutine

#ifdef STDF95
  subroutine cubic_spline_compute_interpolants( interpolator, data_array )
    type(cubic_spline_2d_interpolator), intent(inout) :: interpolator
#else
  subroutine compute_interpolants_cs2d( interpolator, data_array )
    class(cubic_spline_2d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64, dimension(:,:), intent(in) :: data_array
    call compute_spline_2D( data_array, interpolator%spline )
  end subroutine

#ifdef STDF95
  function cubic_spline_interpolate_value( interpolator, eta1, eta2 ) result(val)
    type(cubic_spline_2d_interpolator), intent(in) :: interpolator
#else
  function interpolate_value_cs2d( interpolator, eta1, eta2 ) result(val)
    class(cubic_spline_2d_interpolator), intent(in) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = interpolate_value_2D( eta1, eta2, interpolator%spline )
  end function

#ifdef STDF95
  function cubic_spline_interpolate_derivative_eta1( interpolator, eta1, eta2 ) result(val)
    type(cubic_spline_2d_interpolator), intent(in) :: interpolator
#else
  function interpolate_deriv1_cs2d( interpolator, eta1, eta2 ) result(val)
    class(cubic_spline_2d_interpolator), intent(in) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = interpolate_x1_derivative_2D(eta1,eta2,interpolator%spline)
  end function

#ifdef STDF95
  function cubic_spline_interpolate_derivative_eta2( interpolator, eta1, eta2 ) result(val)
    type(cubic_spline_2d_interpolator), intent(in) :: interpolator
#else
  function interpolate_deriv2_cs2d( interpolator, eta1, eta2 ) result(val)
    class(cubic_spline_2d_interpolator), intent(in) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = interpolate_x2_derivative_2D(eta1,eta2,interpolator%spline)
  end function

  function spline_interpolate2d(this, num_points1, num_points2, data_in, &
                                eta1, eta2) &
       result(data_out)
    class(cubic_spline_2d_interpolator),  intent(in)       :: this
    sll_int32,  intent(in)                 :: num_points1
    sll_int32,  intent(in)                 :: num_points2
    sll_real64, dimension(:,:), intent(in)   :: eta1
    sll_real64, dimension(:,:), intent(in)   :: eta2
    sll_real64, dimension(:,:), intent(in)   :: data_in
    sll_real64, dimension(num_points1,num_points2) :: data_out
    ! local variables
    sll_int32 :: i,j, ierr
    ! compute the interpolating spline coefficients
    call compute_spline_2D( data_in, this%spline )
    do j = 1, num_points2
    do i = 1, num_points1
        data_out(i,j) = this%interpolate_value(eta1(i,j),eta2(i,j))
    end do
    end do

  end function 

  function spline_interpolate2d_disp(this, num_points1, num_points2, data_in, &
                                alpha1, alpha2) &
       result(data_out)
    class(cubic_spline_2d_interpolator),  intent(in)       :: this
    sll_int32,  intent(in)                 :: num_points1
    sll_int32,  intent(in)                 :: num_points2
    sll_real64, dimension(:,:), intent(in)   :: alpha1
    sll_real64, dimension(:,:), intent(in)   :: alpha2
    sll_real64, dimension(:,:), intent(in)   :: data_in
    sll_real64, dimension(num_points1,num_points2) :: data_out
    sll_real64 :: eta1, eta1_min, eta1_max, delta_eta1
    sll_real64 :: eta2, eta2_min, eta2_max, delta_eta2
    ! local variables
    sll_int32 :: i,j, ierr
    ! compute the interpolating spline coefficients

    eta1_min  = this%spline%x1_min 
    eta1_max  = this%spline%x1_max 
    eta2_min  = this%spline%x2_min 
    eta2_max  = this%spline%x2_max 
    delta_eta1 =this%spline%x1_delta  
    delta_eta2 =this%spline%x2_delta  
    
    call compute_spline_2D( data_in, this%spline )

    do j = 1, num_points2
        do i = 1, num_points1
           eta1 = eta1_min + modulo(eta1-eta1_min-alpha1(i,j),eta1_max-eta1_min)
           eta2 = eta2_min + modulo(eta2-eta2_min-alpha2(i,j),eta2_max-eta2_min)
           data_out(i,j) = this%interpolate_value(eta1,eta2)
       end do
    end do

  end function 

end module sll_cubic_spline_interpolator_2d
