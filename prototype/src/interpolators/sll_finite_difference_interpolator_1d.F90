module sll_finite_difference_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_interpolators_1d_base_macros.h"
use sll_module_interpolators_1d_base
  implicit none

  ! This is an experimental type meant to solve the following problem:
  ! We are used to fields coming with their own means of interpolation.
  ! The standard interpolator, represented by the base class, is general
  ! in that it can return interpolated values along a continuous coordinate.
  ! In other words, given discrete data on nodes, the general interpolator
  ! can give the illusion of having a continuous function defined on the 
  ! same domain as the discrete data.
  !
  ! However, sometimes the interpolation services needed are very basic.
  ! If we are only interested in having the derivatives on the nodes of an
  ! array, this is a service provided by the basic interpolator, hence we
  ! can 

  type, extends(sll_interpolator_1d_base) :: finite_difference_1d_interpolator
     sll_int32  :: num_points
     sll_real64 :: delta      ! cell size, distance between data points
     sll_real64 :: r_delta    ! reciprocal of the cell size
   contains
     procedure, pass(interpolator) :: initialize => &
          initialize_finite_difference_1d_interp
     procedure :: compute_interpolants => null_fd_1d_array_msg
     procedure :: interpolate_value => null_fd_1d_arg_msg
     ! FIXME: define is we want this to explicitly state that we are taking
     ! derivatives with respect to eta or whatnot.
!     procedure :: interpolate_derivative_eta1 => null_interp_one_arg_msg
     ! yes, we are using the null functions more than once.
     procedure :: interpolate_derivative_eta1 => null_fd_1d_arg_msg
     procedure :: interpolate_array_values => null_fd_1d_array_sub
     procedure :: interpolate_pointer_values => null_fd_1d_ptr_sub
     procedure :: interpolate_array_derivatives => null_fd_1d_array_sub
     procedure :: interpolate_pointer_derivatives => &
          finite_difference_1d_interp_derivatives 
     procedure, pass :: interpolate_array => null_fd_1d_interp_array
     procedure, pass :: reconstruct_array => null_fd_1d_reconstruct
     procedure, pass :: interpolate_array_disp => null_fd_1d_array_disp
  end type finite_difference_1d_interpolator

contains

  subroutine initialize_finite_difference_1d_interp( &
    interpolator, &
    num_points, &
    delta )

    class(finite_difference_1d_interpolator), intent(inout) :: interpolator
    sll_int32, intent(in) :: num_points
    sll_real64, intent(in) :: delta ! cell spacing
    interpolator%num_points = num_points
    interpolator%delta = delta
    interpolator%r_delta = 1.0_f64/delta
  end subroutine initialize_finite_difference_1d_interp

  subroutine finite_difference_1d_interp_derivatives( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output )

    class(finite_difference_1d_interpolator), intent(in) :: interpolator
    sll_int32, intent(in) :: num_pts
    ! the following input array must be a pointer so that this function
    ! is usable with splitting methods in higher dimensional data.
    sll_real64, dimension(:), pointer :: vals_to_interpolate
    sll_real64, dimension(:), pointer :: output
    sll_int32  :: i
    sll_real64 :: r_delta  
    sll_real64, dimension(:), pointer :: dptr

    SLL_ASSERT( size(vals_to_interpolate)  >= num_pts )
    SLL_ASSERT( size(output) >= num_pts )

    r_delta = interpolator%r_delta

    ! Compute derivative in first point with a forward scheme (-3/2, 2, -1/2)
    output(1) = r_delta*( -1.5_f64*vals_to_interpolate(1) + &
                           2.0_f64*vals_to_interpolate(2) - &
                           0.5_f64*vals_to_interpolate(3) )

    ! Compute derivative in the bulk of the array with a centered scheme.
    do i=2,num_pts-1
       output(i) = 0.5_f64*r_delta*( vals_to_interpolate(i+1) - &
                                     vals_to_interpolate(i-1) )
    end do

    ! Compute derivative in the last point with a backward scheme (1/2, -2, 3/2)
    output(num_pts) = r_delta*( 0.5_f64*vals_to_interpolate(num_pts-2) - &
                                2.0_f64*vals_to_interpolate(num_pts-1) + &
                                1.5_f64*vals_to_interpolate(num_pts) )

  end subroutine finite_difference_1d_interp_derivatives

  ! Define the functions that are not needed in this particular interpolator
  ! with the use of the null-function build macros.

  DEFINE_NULL_INTERP_ONE_ARG_MSG(finite_difference_1d_interpolator, null_fd_1d_arg_msg)

  DEFINE_NULL_INTERP_1D_ARRAY(finite_difference_1d_interpolator, null_fd_1d_interp_array)

  DEFINE_NULL_INTERP_1D_ARRAY_MSG(finite_difference_1d_interpolator, null_fd_1d_array_msg)

  DEFINE_NULL_INTERP_1D_ARRAY_SUB(finite_difference_1d_interpolator, null_fd_1d_array_sub)

  DEFINE_NULL_INTERP_1D_POINTER_SUB(finite_difference_1d_interpolator, null_fd_1d_ptr_sub)

  DEFINE_NULL_RECONSTRUCT_1D_ARRAY(finite_difference_1d_interpolator, null_fd_1d_reconstruct)

  DEFINE_NULL_INTERPOLATE_1D_DISP(finite_difference_1d_interpolator, null_fd_1d_array_disp)

end module sll_finite_difference_interpolator_1d
