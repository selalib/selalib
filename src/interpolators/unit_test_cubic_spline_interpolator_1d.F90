program cubic_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_constants.h"

use util_constants
use sll_module_interpolators_1d_base
use sll_module_cubic_spline_interpolator_1d

implicit none

class(sll_interpolator_1d_base), pointer       :: interp
type(sll_cubic_spline_interpolator_1d), target :: spline

sll_real64                            :: error
sll_real64                            :: phase
sll_real64, allocatable, dimension(:) :: point
sll_real64, allocatable, dimension(:) :: pdata  
sll_real64, allocatable, dimension(:) :: fdata
sll_real64, allocatable, dimension(:) :: coord
sll_real64, allocatable, dimension(:) :: gdata

sll_int32 :: ierr, i
sll_int32, parameter :: n = 512
sll_real64  :: x_min, x_max, delta

SLL_ALLOCATE(pdata(n), ierr)
SLL_ALLOCATE(fdata(n), ierr)
SLL_ALLOCATE(point(n), ierr)
SLL_ALLOCATE(coord(n), ierr)
SLL_ALLOCATE(gdata(n), ierr)

print*, 'Initialize data and point array'
x_min = 0.0_f64
x_max = 2.0_f64 * sll_pi
delta = (x_max - x_min ) / real(n-1,f64) 
phase = 0.0_f64
call random_number(point)
do i=1,n
   coord(i) = (i-1)*delta
   pdata(i) = f(coord(i))
   point(i) = (x_max - x_min) * point(i)
   gdata(i) = f(point(i))
end do

print*, 'Cubic spline interpolation'
call spline%initialize(n, x_min, x_max, SLL_PERIODIC )

interp =>  spline
fdata = interp%interpolate_array(n, pdata, point)

do i = 1, n
   write(28,*) coord(i), pdata(i), point(i), fdata(i)
end do

error = maxval(abs(gdata-fdata))
if ( error < 2e-10) then
  print*, 'Successful, exiting program.'
  print*, 'PASSED'
else
  print*, 'FAILED'
end if
print*, 'Error=', error

contains

function f(x)
sll_real64 :: x
sll_real64 :: f

f = 2.0_f64*(sin(x) + 2.5_f64 + cos(x))

end function f

end program cubic_spline_interpolator_1d
