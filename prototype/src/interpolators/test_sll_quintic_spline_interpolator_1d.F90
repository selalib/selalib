program test_sll_quintic_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use numeric_constants
  use util_constants

#ifndef STDF95
  use sll_module_interpolators_1d_base
#endif
  use WENO_interp
  use sll_quintic_spline_interpolator_1d
    implicit none

#ifdef STDF95
  type(quintic_spline_1d_interpolator), pointer  :: interp
#else
  class(sll_interpolator_1d_base), pointer     :: interp
#endif

  type(quintic_spline_1d_interpolator), target  :: spline
  type(WENO_interp_1d), pointer               :: weno

  sll_real64                            :: error
  sll_real64                            :: phase
  sll_real64, allocatable, dimension(:) :: interpolation_points
  sll_real64, allocatable, dimension(:) :: data  
  sll_real64, allocatable, dimension(:) :: out
  sll_real64, allocatable, dimension(:) :: coordinates_d
  sll_real64, allocatable, dimension(:) :: data_interp
  sll_real64                            :: mu

  sll_int32 :: ierr, i
  sll_int32, parameter :: n = 512
  sll_real64  :: x_min, x_max, delta

  SLL_ALLOCATE(data(n), ierr)
  SLL_ALLOCATE(out(n), ierr)
  SLL_ALLOCATE(interpolation_points(n), ierr)
  SLL_ALLOCATE(coordinates_d(n), ierr)
  SLL_ALLOCATE(data_interp(n), ierr)

  print *, 'initialize data and interpolation_points array'
  x_min = -10.0_f64
  x_max = 10.0_f64 * sll_pi
  mu = (xmin+xmax)/2
  delta = (x_max - x_min ) / real(n-1,f64) 
  phase = 0.0_f64


  do i=1,n
     coordinates_d(i) = (i-1)*delta
     interpolation_points(i) = modulo(coordinates_d(i) - delta/3.0_f64,2.0_f64 * sll_pi)
     data(i)        = exp( - ( coordinates_d(i) - mu )**2  )
     data_interp(i) = exp( - ( interpolation_points(i) - mu )**2  )
  end do

  print*, 'quintic spline interpolation'
#ifdef STDF95
  call quintic_spline_initialize(spline, n, x_min, x_max, 2 )
#else
  call spline%initialize(n, x_min, x_max, 2 )
#endif

  interp =>  spline
#ifdef STDF95
  out = quintic_spline_interpolate_array(interp, n, data, interpolation_points)
#else
  out = interp%interpolate_array(n, data, interpolation_points)
#endif

  error = 0.0_f64
  do i=1,n   
     error = max(error, abs(data_interp(i) - out(i)))
     !print*, i, interpolation_points(i), data_interp(i) - out(i)
  end do
  print*, '    error=',error

  print *, 'Successful, exiting program.'

end program test_sll_quintic_spline_interpolator_1d
