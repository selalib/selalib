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

type(sll_bspline_1d) :: spline

sll_real64, dimension(:), allocatable :: x
sll_real64, dimension(:), allocatable :: y
sll_int32,  parameter                 :: n = 1024
sll_int32                             :: ierr
sll_real64, dimension(:), allocatable :: gtau
sll_real64, dimension(:), allocatable :: htau
sll_real64                            :: err1
sll_real64                            :: err2

sll_int32                             :: i
sll_int32                             :: j
sll_int32,  parameter                 :: k = 4

sll_int32,  parameter :: m = 2
sll_real64 :: tau_min = 0.0_f64
sll_real64 :: tau_max = 1.0_f64
sll_real64 :: slope_min
sll_real64 :: slope_max
sll_real64 :: t0, t1, t2

SLL_ALLOCATE(x(n),ierr)
SLL_ALLOCATE(y(n),ierr)

do i = 1, n
  x(i) = (i-1)*1.0_f64/(n-1)
end do

call initialize_bspline_1d(spline, n, k, tau_min, tau_max, SLL_HERMITE)


SLL_ALLOCATE(gtau(n),ierr)
gtau = cos(2*sll_pi*spline%tau)
slope_min = -sin(2*sll_pi*tau_min)*2*sll_pi
slope_max = -sin(2*sll_pi*tau_max)*2*sll_pi
call cpu_time(t0)
do j = 1, 1000
  call compute_bspline_1d(spline, gtau, slope_min, slope_max)
end do
call cpu_time(t1)
do j = 1, 1000
  call interpolate_array_values( spline, n, x, y)
end do
call cpu_time(t2)
print*, ' time spent to compute interpolants : ', t1-t0
print*, ' time spent to interpolate values : ', t2-t1

print*, "values error = ", maxval(abs(y-cos(2*sll_pi*x)))
call interpolate_array_derivatives( spline, n, x, y)
print*, "derivatives error = ", maxval(abs(y+2*sll_pi*sin(2*sll_pi*x)))

SLL_ALLOCATE(htau(n),ierr)
htau = sin(2*sll_pi*spline%tau)
slope_min = cos(2*sll_pi*tau_min)*2*sll_pi
slope_max = cos(2*sll_pi*tau_max)*2*sll_pi
call update_bspline_1d(spline, htau, slope_min, slope_max)
call interpolate_array_values( spline, n, x, y)
print*, "values error = ", maxval(abs(y-sin(2*sll_pi*x)))
call interpolate_array_derivatives( spline, n, x, y)
print*, "derivatives error = ", maxval(abs(y-2*sll_pi*cos(2*sll_pi*x)))

err1 = 0.0_f64
call cpu_time(t0)
do i = 1, n
  err1 = err1 + abs(interpolate_value(spline,x(i))-sin(2*sll_pi*x(i))) 
end do
call cpu_time(t1)
err2 = 0.0_f64
do i = 1, n
  err2 = err2 + abs(interpolate_derivative(spline,x(i))-2*sll_pi*cos(2*sll_pi*x(i))) 
end do
call cpu_time(t2)

print*, "-------------------------------------------------"
print*, " values error = ", err1 / n
print*, " derivatives error = ", err2 / n
print*, ' time spent to interpolate_value      : ', t1-t0
print*, ' time spent to interpolate_derivative : ', t2-t1
print*, "-------------------------------------------------"



print*, 'PASSED'

end program test_deboor_hermite
