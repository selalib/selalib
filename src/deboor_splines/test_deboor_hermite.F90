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

implicit none


type(sll_bspline_1d) :: bsplines

sll_real64, dimension(:), allocatable :: x
sll_real64, dimension(:), allocatable :: y
sll_int32,  parameter                 :: n = 10
sll_int32,  parameter                 :: k = 4
sll_int32                             :: mflag
sll_int32                             :: ierr
sll_real64, dimension(:), allocatable :: gtau
sll_real64, dimension(:), allocatable :: htau

sll_int32                             :: i
sll_int32                             :: j
sll_int32, parameter                  :: nx = 1000

sll_int32, parameter :: m = 2
sll_real64 :: tau_min = 0.0_f64
sll_real64 :: tau_max = 1.0_f64
sll_real64 :: slope_min
sll_real64 :: slope_max

SLL_ALLOCATE(x(nx),ierr)
SLL_ALLOCATE(y(nx),ierr)

do i = 1, nx
  x(i) = (i-1)*1.0_f64/(nx-1)
end do

call initialize_bspline_1d(bsplines, n, k, tau_min, tau_max)

call build_system(bsplines)

SLL_ALLOCATE(gtau(n),ierr)
gtau = cos(2*sll_pi*bsplines%tau)
slope_min = -sin(2*sll_pi*tau_min)*2*sll_pi
slope_max = -sin(2*sll_pi*tau_max)*2*sll_pi
call compute_coefficients(bsplines, gtau, slope_min, slope_max)
do j = 1, 10000
  call interpolate_array_values( bsplines, nx, x, y)
end do

print*, "error = ", maxval(abs(y-cos(2*sll_pi*x)))

SLL_ALLOCATE(htau(n),ierr)
htau = sin(2*sll_pi*bsplines%tau)
slope_min = cos(2*sll_pi*tau_min)*2*sll_pi
slope_max = cos(2*sll_pi*tau_max)*2*sll_pi
call compute_coefficients(bsplines, htau, slope_min, slope_max)
do j = 1, 10000
  call interpolate_array_values( bsplines, nx, x, y)
end do

print*, "error = ", maxval(abs(y-sin(2*sll_pi*x)))

print*, 'PASSED'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains


subroutine build_system(this)

  type(sll_bspline_1d)    :: this 
  sll_real64              :: dbiatx(k,2)
  sll_real64              :: a(k,k)
  sll_real64              :: taui
  sll_int32               :: kpkm2
  sll_int32               :: left
  sll_int32               :: n
  sll_int32               :: k
  sll_int32               :: iflag
  sll_int32               :: j
  sll_int32               :: jj
  sll_int32               :: l
  
  n = this%n
  k = this%k
  a = 0.0_f64

  kpkm2     = 2*(k-1)
  left      = k
  this%q    = 0.0_f64
  dbiatx    = 0.0_f64
  
  SLL_ASSERT(m < n) 

  l = 0 ! index for the derivative

  do i = 1, n
      
    taui = this%tau(i)
    call interv( this%t, n+m+k, taui, left, mflag )


    if (i < n) then

      call bsplvb ( this%t, k, 1, taui, left, this%bcoef )
      jj = i-left+1+(left-k)*(k+k-1)+l
      do j = 1, k
        jj = jj + kpkm2
        this%q(jj) = this%bcoef(j)
      end do
   
      if ( i == 1 ) then   
        call bsplvd( this%t, k, taui, left, a, dbiatx, 2)
        l = l + 1
        jj = i-left+1+(left-k)*(k+k-1)+l
        do j = 1, k
          jj = jj + kpkm2
          this%q(jj) = dbiatx(j,2)
        end do
      end if

    else

      call bsplvd( this%t, k, taui, left, a, dbiatx, 2)
      jj = i-left+1+(left-k)*(k+k-1)+l
      do j = 1, k
        jj = jj + kpkm2
        this%q(jj) = dbiatx(j,2)
      end do
      l = l + 1
      
      call bsplvb ( this%t, k, 1, taui, left, this%bcoef )
      jj = i-left+1+(left-k)*(k+k-1)+l
      do j = 1, k
        jj = jj + kpkm2 
        this%q(jj) = this%bcoef(j)
      end do

    end if
 
  end do
  
  !Obtain factorization of A, stored again in Q.

  call banfac ( this%q, k+k-1, n+m, k-1, k-1, iflag )

end subroutine build_system

subroutine compute_coefficients( this, gtau, slope_min, slope_max)

  type(sll_bspline_1d)    :: this 
  sll_real64, intent(in)  :: gtau(:)
  sll_real64, intent(in)  :: slope_min
  sll_real64, intent(in)  :: slope_max

  sll_int32               :: n
  sll_int32               :: k

  n = this%n
  k = this%k


  SLL_ASSERT(size(gtau) == this%n)

  this%bcoef_spline(1)   = gtau(1)
  this%bcoef_spline(2)   = slope_min
  this%bcoef_spline(3:n) = gtau(2:n-1)
  this%bcoef_spline(n+1) = slope_max
  this%bcoef_spline(n+2) = gtau(n)

  call banslv ( this%q, k+k-1, n+m, k-1, k-1, this%bcoef_spline )
  
end subroutine compute_coefficients

subroutine interpolate_array_values( this,      &
                                     nx,        &
                                     x,         &
                                     y)
  type(sll_bspline_1d)    :: this 
  sll_int32,  intent(in)  :: nx
  sll_real64, intent(in)  :: x(nx)
  sll_real64, intent(out) :: y(nx)

  sll_int32               :: i
  sll_int32               :: n
  sll_int32               :: k

  n = this%n
  k = this%k

  do i = 1, nx
    y(i) = bvalue( this%t, this%bcoef_spline, n+m, k, x(i), 0)
  end do

end subroutine interpolate_array_values

end program test_deboor_hermite
