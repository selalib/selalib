program unit_test
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use sll_constants
  use util_constants

  use sll_module_interpolators_1d_base
  !use WENO_interp
  use sll_cubic_spline_interpolator_1d
  use sll_quintic_spline_interpolator_1d
  use cubic_non_uniform_spline_interpolator_1d
 ! use sll_periodic_interpolator_1d
  implicit none

  class(sll_interpolator_1d_base), pointer     :: interp

  type(cubic_spline_1d_interpolator), target   :: spline
  !type(quintic_spline_1d_interpolator), target :: quintic_spline
  !type(cubic_non_uniform_spline_1d_interpolator), target  :: cubic_nonunif_spline
  !type(WENO_interp_1d), pointer               :: weno

  sll_real64                            :: error
  sll_real64                            :: phase
  sll_real64, allocatable, dimension(:) :: interpolation_points
  sll_real64, allocatable, dimension(:) :: data  
  sll_real64, allocatable, dimension(:) :: out
  sll_real64, allocatable, dimension(:) :: coordinates_d
  sll_real64, allocatable, dimension(:) :: data_interp

  sll_int32 :: ierr, i
  sll_int32, parameter :: n = 512
  sll_real64  :: x_min, x_max, delta

  SLL_ALLOCATE(data(n), ierr)
  SLL_ALLOCATE(out(n), ierr)
  SLL_ALLOCATE(interpolation_points(n), ierr)
  SLL_ALLOCATE(coordinates_d(n), ierr)
  SLL_ALLOCATE(data_interp(n), ierr)

  print *, 'initialize data and interpolation_points array'
  x_min = 0.0_f64
  x_max = 2.0_f64 * sll_pi
  delta = (x_max - x_min ) / real(n-1,f64) 
  phase = 0.0_f64
  do i=1,n
     coordinates_d(i) = (i-1)*delta
     interpolation_points(i) = modulo(coordinates_d(i) - delta/3.0_f64,2.0_f64 * sll_pi)
     data(i)        = 2.0_f64*(sin(coordinates_d(i)) + 2.5_f64 + cos(coordinates_d(i)))
     data_interp(i) = 2.0_f64*(sin(interpolation_points(i)) + 2.5_f64 &
          + cos(interpolation_points(i)))
  end do

  print*, 'Cubic spline interpolation'
  call spline%initialize(n, x_min, x_max, SLL_PERIODIC )
  !call quintic_spline%initialize(n, x_min, x_max, SLL_PERIODIC )
  !call cubic_nonunif_spline%initialize(n, x_min, x_max, SLL_PERIODIC )

  interp =>  spline
!  interp =>  quintic_spline
!  interp =>  cubic_nonunif_spline
  out = interp%interpolate_array(n, data, interpolation_points)


  error = 0.0_f64
  do i=1,n   
     error = max(error, abs(data_interp(i) - out(i)))
    ! print*, i, interpolation_points(i), data_interp(i) - out(i)
  end do
  print*, '    error=',error

!!$  print*, 'WENO interpolation'
!!$  weno = new_WENO_1D( n, x_min, x_max )
!!$  interp => weno
!!$  out = interp%interpolate_array(n, data, interpolation_points)
!!$  !print*, 'delta ', delta , interp%weno%delta, x_min, coordinates_d(1), x_max, coordinates_d(n)
!!$  error = 0.0_f64
!!$  do i=1,n   
!!$     error = max(error, abs(data_interp(i) - out(i)))
!!$     !print*, i, interpolation_points(i), data_interp(i) - out(i)
!!$  end do
!!$  
!!$  print*, '    error=',error

  print *, 'Successful, exiting program.'

end program unit_test
