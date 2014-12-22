program quintic_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_constants.h"

use util_constants
use sll_module_interpolators_1d_base
use sll_module_quintic_spline_interpolator_1d

implicit none

class(sll_interpolator_1d_base), pointer         :: interp
type(sll_quintic_spline_interpolator_1d), target :: spline

sll_real64                            :: error
sll_real64, allocatable, dimension(:) :: inte_points
sll_real64, allocatable, dimension(:) :: data  
sll_real64, allocatable, dimension(:) :: out
sll_real64, allocatable, dimension(:) :: coordinates
sll_real64, allocatable, dimension(:) :: data_interp

sll_int32 :: ierr, i
sll_int32, parameter :: n = 512
sll_real64  :: x_min, x_max, delta

SLL_ALLOCATE(data(n), ierr)
SLL_ALLOCATE(out(n), ierr)
SLL_ALLOCATE(inte_points(n), ierr)
SLL_ALLOCATE(coordinates(n), ierr)
SLL_ALLOCATE(data_interp(n), ierr)

print *, 'initialize data and interpolation_points array'
x_min = 0.0_f64
x_max = 2.0_f64 * sll_pi
delta = (x_max - x_min ) / real(n-1,f64) 
do i=1,n
  coordinates(i) = (i-1)*delta
  inte_points(i) = modulo(coordinates(i) - delta/3.0_f64,2.0_f64 * sll_pi)
  data(i)        = 2.0_f64*(sin(coordinates(i)) + 2.5_f64 + cos(coordinates(i)))
  data_interp(i) = 2.0_f64*(sin(inte_points(i)) + 2.5_f64 + cos(inte_points(i)))
end do

print*, 'Quintic spline interpolation'
call spline%initialize(n, x_min, x_max, SLL_PERIODIC, SLL_PERIODIC )

interp =>  spline
out = interp%interpolate_array(n, data, inte_points)

error = maxval(abs(data_interp-out))
print*, 'Error=',error
print*, 'Successful, exiting program.'
print*, 'PASSED'

end program quintic_spline_interpolator_1d
