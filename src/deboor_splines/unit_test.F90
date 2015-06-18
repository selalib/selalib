program unit_test
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_assert.h"
use sll_module_deboor_splines_1d

implicit none

type :: sll_bsplines

  sll_int32               :: n
  sll_int32               :: k
  sll_real64, allocatable :: tau(:)
  sll_real64, pointer     :: t(:)
  sll_real64, allocatable :: q(:)
  sll_real64, pointer     :: bcoef(:)

end type sll_bsplines

type(sll_bsplines) :: bsplines

sll_real64, dimension(:), allocatable :: x
sll_real64, dimension(:), allocatable :: y
sll_int32                             :: n
sll_int32                             :: m
sll_int32                             :: k
sll_int32                             :: mflag
sll_int32                             :: ierr
sll_real64, dimension(:), allocatable :: gtau
sll_real64, dimension(:), allocatable :: htau

sll_int32                             :: i
sll_int32                             :: nx

sll_real64 :: tau_min = 0.0_f64
sll_real64 :: tau_max = 1.0_f64
sll_real64 :: slope_min
sll_real64 :: slope_max

! TAU(N),      the data point abscissas. TAU should be strictly increasing.
! TAU_der(M),  the node index to evaluate the derivative.
! GTAU(N),     the data ordinates.
! GTAU_der(M), the data ordinates.
! T(N+K+M),    the knot sequence.
! N,           the number of data points for the interpolation.
! M,           the number of data points for the derivative.
! K,           the order of the spline.
!
! Q((2*K-1)*(N+M)) is the triangular factorization
! of the coefficient matrix of the linear system for the B-coefficients 
! of the spline interpolant.  The B-coefficients for the interpolant 
! of an additional data set can be obtained without going through all 
! the calculations in this routine, simply by loading HTAU into BCOEF 
! and then executing the call:
!   call banslv ( q, 2*k-1, n+m, k-1, k-1, bcoef )
! Output, real ( kind = 8 ) BCOEF(N+M), the B-spline coefficients of 
! the interpolant.

nx = 100
SLL_ALLOCATE(x(nx),ierr)
SLL_ALLOCATE(y(nx),ierr)

do i = 1, nx
  x(i) = (i-1)*1.0_f64/(nx-1)
end do

k = 4
n = 10
m = 2

call initialize_bsplines(bsplines, n, k, tau_min, tau_max)

call compute_interpolants(bsplines)

SLL_ALLOCATE(gtau(n),ierr)
gtau = cos(2*sll_pi*bsplines%tau)
call interpolate_array_values( bsplines, gtau, nx, x, y)
do i = 1, nx
  write(12,*) x(i), y(i), cos(2*sll_pi*x(i))
end do

SLL_ALLOCATE(htau(n),ierr)
htau = sin(2*sll_pi*bsplines%tau)
call interpolate_array_values( bsplines, htau, nx, x, y)
do i = 1, nx
  write(13,*) x(i), y(i), sin(2*sll_pi*x(i))
end do

slope_min = -sin(2*sll_pi*tau_min)*2*sll_pi
slope_max = -sin(2*sll_pi*tau_max)*2*sll_pi

call compute_interpolants_with_slopes(bsplines)
call interpolate_array_values_with_slopes( bsplines, htau, nx, x, slope_min, slope_max, y)
    
print*, 'PASSED'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

subroutine initialize_bsplines(this, n, k, tau_min, tau_max)

  type(sll_bsplines), intent(out) :: this
  sll_int32         , intent(in)  :: n
  sll_int32         , intent(in)  :: k
  sll_real64        , intent(in)  :: tau_min
  sll_real64        , intent(in)  :: tau_max

  sll_int32                       :: ierr

  this%n = n
  this%k = k

  SLL_ALLOCATE(this%tau(n), ierr)
  do i = 1, n
    this%tau(i) = tau_min + (i-1) * (tau_max-tau_min) / (n-1)
  end do

  SLL_ALLOCATE(this%t(n+k), ierr)

  this%t(1:k) = this%tau(1)
  if ( mod(k,2) == 0 ) then
    do i = k+1,n
      this%t(i) = this%tau(i-k/2)
    end do
  else
    do i = k+1, n
      this%t(i) = 0.5*(this%tau(i-(k-1)/2)+this%tau(i-1-(k-1)/2))
    end do
  end if
  do i = n+1,n+k
    this%t(i) = this%tau(n)
  end do

  SLL_ALLOCATE(this%q(1:(2*k-1)*n),ierr)
  SLL_ALLOCATE(this%bcoef(n),ierr)

end subroutine initialize_bsplines

subroutine compute_interpolants(this)
  
  type(sll_bsplines)     :: this 

  sll_int32              :: n
  sll_int32              :: k
  sll_int32              :: j
  sll_int32              :: jj
  sll_int32              :: kpkm2
  sll_int32              :: left
  sll_int32              :: iflag
  sll_real64             :: taui

  n = this%n
  k = this%k

  kpkm2   = 2 * ( k - 1 )
  left    = k
  this%q  = 0.0D+00
  
  do i = 1, n !  Loop over I to construct the N interpolation equations.
    
    taui = this%tau(i)
    call interv( this%t, n+k, taui, left, mflag )
    call bsplvb ( this%t, k, 1, taui, left, this%bcoef )
    jj = i - left + 1 + ( left - k ) * ( k + k - 1 )
    do j = 1, k
      jj = jj + kpkm2
      this%q(jj) = this%bcoef(j)
    end do
      
  end do
  !
  !  Obtain factorization of A, stored again in Q.
  !
  call banfac ( this%q, k+k-1, n, k-1, k-1, iflag )
    
  if ( iflag == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINT - Fatal Error!'
    write ( *, '(a)' ) '  The linear system is not invertible!'
    stop
  end if

end subroutine compute_interpolants


subroutine interpolate_array_values( this, gtau, nx, x, y)

  type(sll_bsplines)      :: this 
  sll_real64, intent(in)  :: gtau(:)
  sll_int32,  intent(in)  :: nx
  sll_real64, intent(in)  :: x(nx)
  sll_real64, intent(out) :: y(nx)

  sll_int32              :: i
  
  SLL_ASSERT(size(gtau) == this%n)

  !  Solve A * BCOEF = GTAU by back substitution.
  this%bcoef(1:n) = gtau(1:n)
  
  call banslv( this%q, this%k+this%k-1, n, this%k-1, this%k-1, this%bcoef )

  do i = 1, nx
    y(i) = bvalue( this%t, this%bcoef, this%n, this%k, x(i), 0)
  end do

end subroutine interpolate_array_values

subroutine compute_interpolants_with_slopes(this)

  type(sll_bsplines)      :: this 
  sll_real64              :: dbiatx(k,2)
  sll_real64              :: a(k,k)
  sll_real64              :: taui
  sll_int32               :: kpkm2
  sll_int32               :: left
  sll_int32               :: n
  sll_int32, parameter    :: m = 2
  sll_int32               :: k
  sll_int32               :: iflag
  sll_int32               :: ierr
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

  l = 1 ! index for the derivative

  do i = 1, n-1
      
    taui = this%tau(i)
    call interv( this%t, n+m+k, taui, left, mflag )
    call bsplvb ( this%t, k, 1, taui, left, this%bcoef )
    jj = i - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
    
    do j = 1, k
      jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
      this%q(jj) = this%bcoef(j)
    end do
   
    if ( i == 1 ) then   
      call bsplvd( this%t, k, taui, left, a, dbiatx, 2)
      l = l + 1
      jj = i - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
      do j = 1, k
        jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
        this%q(jj) = dbiatx(j,2)
      end do
    end if
   
  end do
  
  taui = this%tau(n)
  call interv( this%t, n+m+k, taui, left, mflag )
  call bsplvd( this%t, k, taui, left, a, dbiatx, 2)
  jj = n - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
  do j = 1, k
     jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
     this%q(jj) = dbiatx(j,2)
  end do
  l = l + 1
  
  call bsplvb ( this%t, k, 1, taui, left, this%bcoef )
  jj = n - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
     
  do j = 1, k
    jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
    this%q(jj) = this%bcoef(j)
  end do
  
  !Obtain factorization of A, stored again in Q.

  call banfac ( this%q, k+k-1, n+m, k-1, k-1, iflag )

end subroutine compute_interpolants_with_slopes

subroutine interpolate_array_values_with_slopes( this,      &
                                                 gtau,      &
                                                 nx,        &
                                                 x,         &
                                                 slope_min, &
                                                 slope_max, &
                                                 y)

  type(sll_bsplines)      :: this 
  sll_real64, intent(in)  :: gtau(:)
  sll_int32,  intent(in)  :: nx
  sll_real64, intent(in)  :: x(nx)
  sll_real64, intent(out) :: y(nx)
  sll_real64, intent(in)  :: slope_min
  sll_real64, intent(in)  :: slope_max

  sll_int32               :: i
  sll_real64, allocatable :: bcoef_spline(:)

  n = this%n
  k = this%k

  SLL_ALLOCATE(bcoef_spline(n+k+m), ierr)

  SLL_ASSERT(size(gtau) == this%n)

  bcoef_spline(2)     = slope_min
  bcoef_spline(n+k-1) = slope_max

  call banslv ( this%q, k+k-1, n+m, k-1, k-1, bcoef_spline )
  
  do i = 1, nx
    y(i) = bvalue( this%t, this%bcoef, this%n, this%k, x(i), 0)
  end do

end subroutine interpolate_array_values_with_slopes

end program unit_test
