! N,           the number of data points for the interpolation.
! K,           the order of the spline.
! TAU(N),      the data point abscissas. TAU should be strictly increasing.
! GTAU(N),     the data ordinates.
! T(N+K+2),    the knot sequence.
!
! Q((2*K-1)*(N+2)) is the triangular factorization
! of the coefficient matrix of the linear system for the B-coefficients 
! BCOEF(N+M), the B-spline coefficients of the interpolant.

program test_deboor_hermite
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_assert.h"
#include "sll_deboor_splines.h"
#include "sll_boundary_condition_descriptors.h"

implicit none


type(sll_bspline_1d) :: bsplines

sll_real64, dimension(:), allocatable :: x
sll_real64, dimension(:), allocatable :: y
sll_int32,  parameter                 :: n = 1024
sll_int32,  parameter                 :: k = 4
sll_int32                             :: ierr
sll_real64, dimension(:), allocatable :: gtau
sll_real64, dimension(:), allocatable :: htau

sll_int32                             :: i
sll_int32                             :: j

sll_int32, parameter :: m = 2
sll_real64 :: tau_min = 0.0_f64
sll_real64 :: tau_max = 1.0_f64
sll_real64 :: slope_min
sll_real64 :: slope_max

SLL_ALLOCATE(x(n),ierr)
SLL_ALLOCATE(y(n),ierr)

do i = 1, n
  x(i) = (i-1)*1.0_f64/(n-1)
end do

call initialize_bspline_1d(bsplines, n, k, tau_min, tau_max, SLL_HERMITE)

call build_system(bsplines)

SLL_ALLOCATE(gtau(n),ierr)
gtau = cos(2*sll_pi*bsplines%tau)
slope_min = -sin(2*sll_pi*tau_min)*2*sll_pi
slope_max = -sin(2*sll_pi*tau_max)*2*sll_pi
call compute_bspline_1d(bsplines, gtau, slope_min, slope_max)
do j = 1, 10000
  call interpolate_array_values( bsplines, n, x, y)
end do

print*, "values error = ", maxval(abs(y-cos(2*sll_pi*x)))
call interpolate_array_derivatives( bsplines, n, x, y)
print*, "derivatives error = ", maxval(abs(y+2*sll_pi*sin(2*sll_pi*x)))

SLL_ALLOCATE(htau(n),ierr)
htau = sin(2*sll_pi*bsplines%tau)
slope_min = cos(2*sll_pi*tau_min)*2*sll_pi
slope_max = cos(2*sll_pi*tau_max)*2*sll_pi
call compute_bspline_1d(bsplines, htau, slope_min, slope_max)
do j = 1, 10000
  call interpolate_array_values( bsplines, n, x, y)
end do

print*, "values error = ", maxval(abs(y-sin(2*sll_pi*x)))
call interpolate_array_derivatives( bsplines, n, x, y)
print*, "derivatives error = ", maxval(abs(y-2*sll_pi*cos(2*sll_pi*x)))

print*, 'PASSED'

end program test_deboor_hermite
