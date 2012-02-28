program unit_test
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use numeric_constants
  use util_constants
  use sll_interpolator_1d
    implicit none

  type(interpolator_1d), pointer        :: interp
  sll_real64                            :: error
  sll_real64                            :: phase
  sll_real64, allocatable, dimension(:) :: coordinates
  sll_real64, allocatable, dimension(:) :: data  
  sll_real64, allocatable, dimension(:) :: out

  sll_int32 :: ierr, i
  sll_int32, parameter :: n = 48
  !sll_real64  :: xmax

  SLL_ALLOCATE(data(n), ierr)
  SLL_ALLOCATE(out(n), ierr)
  SLL_ALLOCATE(coordinates(n), ierr)
  interp =>  new_interpolator_1d('spline', n, 0.0_f64, 2*sll_pi, PERIODIC_SPLINE )


  print *, 'initialize data and coordinates array'
  phase = 0_f64
  do i=1,n
     data(i)        = cos(phase)
     coordinates(i) = modulo(phase + 0.1_f64, 2*sll_pi)
     phase          = phase + interp%spl%delta
  end do

  print*, 'Cubic spline interpolation'
  out = interp%interpolate_1d(interp, n, coordinates, data)

  do i=1,n   
     error = max(error, abs(cos(coordinates(i)) - out(i)))
     !print*, i, coordinates(i), cos(coordinates(i))- out(i)
  end do

  print*, 'error=',error

  print *, 'Successful, exiting program.'

end program unit_test
