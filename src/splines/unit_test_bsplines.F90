! N,           the number of data points for the interpolation.
! K,           the order of the spline.
! TAU(N),      the data point abscissas. TAU should be strictly increasing.
! GTAU(N),     the data ordinates.
! T(N+K+2),    the knot sequence.
!
! Q((2*K-1)*(N+2)) is the triangular factorization
! of the coefficient matrix of the linear system for the B-coefficients 
! BCOEF(N+M), the B-spline coefficients of the interpolant.

program test_bsplines
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_assert.h"
#include "sll_boundary_condition_descriptors.h"
use sll_bsplines

implicit none

type(sll_bspline_1d),      pointer :: bspline


sll_real64, dimension(:), allocatable :: x
sll_real64, dimension(:), allocatable :: y
sll_int32,  parameter                 :: n = 256
sll_int32                             :: ierr
sll_real64, dimension(:), allocatable :: gtau
sll_real64, dimension(:), allocatable :: htau
sll_real64                            :: err1
sll_real64                            :: err2

sll_int32                             :: i
sll_int32                             :: j
sll_int32,  parameter                 :: d = 3
sll_real64                            :: h
sll_int32,  parameter                 :: nstep = 10000

sll_int32,  parameter :: m = 2
sll_real64 :: tau_min = 0.0_f64
sll_real64 :: tau_max = 1.0_f64
sll_real64 :: slope_min
sll_real64 :: slope_max
sll_real64 :: t0, t1, t2, t3

SLL_ALLOCATE(x(n),ierr)
SLL_ALLOCATE(y(n),ierr)
SLL_ALLOCATE(gtau(n),ierr)
SLL_ALLOCATE(htau(n),ierr)

h = 1.0_f64/(n-1)
do i = 1, n
  x(i) = (i-1)*h
end do

print*,'*** PERIODIC ***'
call test_process(SLL_PERIODIC)
print*,'*** HERMITE ***'
call test_process(SLL_HERMITE)

print*, 'PASSED'

contains

subroutine test_process(bc_type)

  sll_int32, intent(in) :: bc_type
  
  if (bc_type == SLL_PERIODIC) print*, "Periodic Bspline"
  bspline => new_bspline_1d( n, d, tau_min, tau_max, bc_type)
  print*, 'bspline allocated'
  
  gtau = cos(2*sll_pi*bspline%tau)
  slope_min = -sin(2*sll_pi*tau_min)*2*sll_pi
  slope_max = -sin(2*sll_pi*tau_max)*2*sll_pi
  call cpu_time(t0)
  call compute_bspline_1d(bspline, gtau, slope_min, slope_max)
  call cpu_time(t1)
  do j = 1,nstep
    call interpolate_array_values( bspline, n, x, y)
  end do
  print*, "average values error      = ", sum(abs(y-cos(2*sll_pi*x)))/n
  print*, "maximum values error      = ", maxval(abs(y-cos(2*sll_pi*x)))
  call cpu_time(t2)
  do j = 1,nstep
    call interpolate_array_derivatives( bspline, n, x, y)
  end do
  print*, "average derivatives error = ", sum(abs(y+2*sll_pi*sin(2*sll_pi*x)))/n
  print*, "maximum derivatives error = ", maxval(abs(y+2*sll_pi*sin(2*sll_pi*x)))
  call cpu_time(t3)

  print*, ' time spent to compute interpolants          : ', t1-t0
  print*, ' time spent to interpolate array values      : ', t2-t1
  print*, ' time spent to interpolate array derivatives : ', t3-t2
  
  htau = sin(2*sll_pi*bspline%tau)
  call update_bspline_1d(bspline, htau)
  call random_number(x)
  x = x * (tau_max-tau_min)
  call interpolate_array_values( bspline, n, x, y)
  print*, "L2 norm error = ", sqrt(sum((y-sin(2*sll_pi*x))**2*h))
  call interpolate_array_derivatives( bspline, n, x, y)
  print*, "H1 norm error = ", sqrt(sum((y-2*sll_pi*cos(2*sll_pi*x))**2*h))
  
  call cpu_time(t0)
  do j = 1,nstep
    err1 = 0.0_f64
    do i = 1, n
      err1 = err1 + abs(interpolate_value(bspline,x(i))-sin(2*sll_pi*x(i))) 
    end do
  end do
  call cpu_time(t1)
  do j = 1,nstep
    err2 = 0.0_f64
    do i = 1, n
      err2 = err2 + abs(interpolate_derivative(bspline,x(i))-2*sll_pi*cos(2*sll_pi*x(i))) 
    end do
  end do
  call cpu_time(t2)
  
  print*, "-------------------------------------------------"
  print*, " values error = ", err1 / n
  print*, " derivatives error = ", err2 / n
  print*, ' time spent in interpolate_value      : ', t1-t0
  print*, ' time spent in interpolate_derivative : ', t2-t1
  print*, "-------------------------------------------------"

end subroutine test_process


end program test_bsplines
