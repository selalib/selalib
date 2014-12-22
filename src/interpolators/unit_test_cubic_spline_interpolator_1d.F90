program cubic_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_constants.h"

use util_constants
use sll_module_interpolators_1d_base
use sll_module_cubic_spline_interpolator_1d
use sll_module_cubic_spline_interpolator_1d_nonuniform
use sll_module_periodic_interpolator_1d
implicit none

class(sll_interpolator_1d_base),        pointer :: interp
type(sll_cubic_spline_interpolator_1d), target  :: spline

sll_real64                            :: error
sll_real64, allocatable, dimension(:) :: points
sll_real64, allocatable, dimension(:) :: data  
sll_real64, allocatable, dimension(:) :: out
sll_real64, allocatable, dimension(:) :: coords
sll_real64, allocatable, dimension(:) :: data_interp

sll_int32, parameter :: n = 512
sll_int32            :: ierr, i
sll_real64           :: x_min, x_max, delta

SLL_ALLOCATE(data(n), ierr)
SLL_ALLOCATE(out(n), ierr)
SLL_ALLOCATE(points(n), ierr)
SLL_ALLOCATE(coords(n), ierr)
SLL_ALLOCATE(data_interp(n), ierr)

print *, 'initialize data and points array'
x_min = 0.0_f64
x_max = 2.0_f64 * sll_pi
delta = (x_max - x_min ) / real(n-1,f64) 

do i=1,n
  coords(i) = (i-1)*delta
  points(i) = modulo(coords(i) - delta/3.0_f64,2.0_f64 * sll_pi)
  data(i)        = 2.0_f64*(sin(coords(i)) + 2.5_f64 + cos(coords(i)))
  data_interp(i) = 2.0_f64*(sin(points(i)) + 2.5_f64 + cos(points(i)))
end do

print*, 'Cubic spline interpolation'
call spline%initialize(n, x_min, x_max, SLL_PERIODIC )

interp =>  spline
out = interp%interpolate_array(n, data, points)

error = maxval(abs(data_interp - out))
print*, '    error=',error

print *, 'Successful, exiting program.'
print *, 'PASSED'

end program cubic_spline_interpolator_1d
