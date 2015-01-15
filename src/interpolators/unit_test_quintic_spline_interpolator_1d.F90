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
sll_real64, allocatable, dimension(:) :: point
sll_real64, allocatable, dimension(:) :: pdata  
sll_real64, allocatable, dimension(:) :: fdata
sll_real64, allocatable, dimension(:) :: coord
sll_real64, allocatable, dimension(:) :: gdata

sll_int32 :: ierr, i
sll_int32, parameter :: n = 32
sll_int32, parameter :: m = 512
sll_real64  :: x_min, x_max, delta

SLL_ALLOCATE(coord(n), ierr)
SLL_ALLOCATE(pdata(n), ierr)
SLL_ALLOCATE(fdata(m), ierr)
SLL_ALLOCATE(gdata(m), ierr)
SLL_ALLOCATE(point(m), ierr)

print *, 'initialize data and interpolation_points array'
x_min = 0.0_f64
x_max = 2.0_f64 * sll_pi
delta = (x_max - x_min ) / real(n-1,f64) 
do i=1,n
  coord(i) = (i-1)*delta
  pdata(i) = f(coord(i))
end do

print*, 'Quintic spline interpolation'
call spline%initialize(n, x_min, x_max, SLL_DIRICHLET, SLL_DIRICHLET )

call set_values_at_boundary( spline,       &
                             f(coord(1)),  &
                             f(coord(n)),  &
                             df(coord(1)), &
                             df(coord(n)))
interp => spline

call interp%compute_interpolants(pdata)

delta = (x_max - x_min ) / real(m-1,f64) 
do i=1,m
  point(i) = (i-1)*delta
  gdata(i) = f(point(i))
  fdata(i) = interp%interpolate_value(point(i))
end do

error = maxval(abs(fdata-gdata))
print*, 'Error=',error
print*, 'Successful, exiting program.'
print*, 'PASSED'

contains

function f(x)

  real(8) :: x
  real(8) :: f
  
  f = 2.0_f64*(sin(x) + 2.5_f64 + cos(x))
  
end function

function df(x)

  real(8) :: x
  real(8) :: df
  
  df = 2.0_f64*(cos(x) - sin(x))
  
end function

end program quintic_spline_interpolator_1d
