program unit_test
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_constants.h"
use sll_module_deboor_splines_1d

implicit none

sll_real64, dimension(:), pointer     :: t
sll_real64, dimension(:), allocatable :: f
sll_real64, dimension(:), allocatable :: df
sll_real64, dimension(:), allocatable :: tau
sll_int32,  dimension(:), allocatable :: tau_der
sll_int32                             :: n
sll_int32                             :: m
sll_int32                             :: k
sll_real64, dimension(:), allocatable :: q
sll_real64, dimension(:),   pointer   :: bcoef
sll_real64, dimension(:),   pointer   :: bcoef_spline 
sll_real64, dimension(:,:), pointer   :: bcoef_der
sll_real64, dimension(:,:), pointer   :: a
sll_int32                             :: iflag
sll_int32                             :: mflag
sll_int32                             :: ierr
sll_real64                            :: y
sll_int32                             :: i
sll_int32                             :: j
sll_int32                             :: l
sll_int32                             :: left
sll_int32                             :: kpkm2
sll_int32                             :: jj
sll_real64, dimension(:), allocatable ::  gtau ! (n)
sll_real64, dimension(:), allocatable ::  gtau_der ! (n)

sll_real64:: taui
sll_real64:: taui_der

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


k = 4
n = 10
m = 2

SLL_ALLOCATE(t(n+k+m),ierr)
SLL_ALLOCATE(f(n),ierr)
SLL_ALLOCATE(df(m),ierr)
SLL_ALLOCATE(tau(n),ierr)
SLL_ALLOCATE(gtau(n),ierr)
SLL_ALLOCATE(gtau_der(n),ierr)
SLL_ALLOCATE(tau_der(m),ierr)
SLL_ALLOCATE(q((2*k-1)*(n+m)),ierr)
SLL_ALLOCATE(bcoef(n+m),ierr)
SLL_ALLOCATE(bcoef_spline(n+m),ierr)
SLL_ALLOCATE(bcoef_der(k,2),ierr)
SLL_ALLOCATE(a(k,k),ierr); a = 0.0_f64

!PN : Set point values and abscissae
do i = 1 , n
 tau(i) = 0.0_f64 + (i-1)*1./(n-1)
 f(i) = cos(2*sll_pi*tau(i))
end do

!PN : Initialize Knots
t(1:k) = tau(1)
do i = k+1, n+k-2
  t(i) = tau(i-k+1)
end do
do i = n+k-1,n+k+m
  t(i) = tau(n)
end do

tau_der(1)   = 1
df(1) = -sin(2*sll_pi*tau(1))*2*sll_pi
tau_der(2)   = n
df(2) = -sin(2*sll_pi*tau(n))*2*sll_pi

call splint_der(tau,f,tau_der,df,t,n,m,k,q,bcoef,iflag ) 

!PN : interpolate values
print*, " Values "
do i = 1 , n
  y=bvalue( t, bcoef, n+m, k, tau(i), 0)
  print"(3(a8,f17.12))", 'approx',y,'exact', f(i), 'error', y-f(i)
end do

!PN : interpolate values of first derivative
!print*, " Values of first derivative "
do i = 1,n
  y=bvalue( t, bcoef, n+m, k, tau(i), 1)
  print"(3(a8,f17.12))",'approx',y                              &
                       ,'exact', -sin(2*sll_pi*tau(i))*2*sll_pi &
                       ,'error',y+sin(2*sll_pi*tau(i))*2*sll_pi
end do
    
kpkm2 = 2*(k-1)
left  = k
q(1:(2*k-1)*(n+m)) = 0.0_f64
bcoef_der(1:k,1:2) = 0.0_f64

! we must suppose that m is <= than n 
if (m > n) then
  print*, 'problem m must be < = at n'
  print*, 'value m =', m, 'value n =', n
  stop
end if
l = 1 ! index for the derivative
!
!  Loop over I to construct the N interpolation equations.
!
do i = 1, n-1
   
  taui = tau(i)
  
  !
  !  Find LEFT in the closed interval (I,I+K-1) such that
  !
  !    T(LEFT) <= TAU(I) < T(LEFT+1)
  !
  !  The matrix is singular if this is not possible.
  !  With help of the Schoenberg-Whitney theorem 
  !  we can prove that if the diagonal of the 
  !  matrix B_j(x_i) is null, we have a non-inversible matrix.  

  call interv( t, n+m+k, taui, left, mflag )

  !  The I-th equation enforces interpolation at TAUI, hence for all J,
  !
  !    A(I,J) = B(J,K,T)(TAUI).
  !
  !Only the K entries with J = LEFT-K+1,...,LEFT actually might be nonzero.
  !
  !These K numbers are returned, in BCOEF (used for temporary storage here),
  !  by the following.
  !
  
  call bsplvb ( t, k, 1, taui, left, bcoef )
  
  !  We therefore want BCOEF(J) = B(LEFT-K+J)(TAUI) to go into
  !  A(I,LEFT-K+J), that is, into Q(I-(LEFT+J)+2*K,(LEFT+J)-K) since
  !  A(I+J,J) is to go into Q(I+K,J), for all I, J, if we consider Q
  !  as a two-dimensional array, with  2*K-1 rows.  See comments in
  !  BANFAC.
  !
  !  In the present program, we treat Q as an equivalent
  !  one-dimensional array, because of fortran restrictions on
  !  dimension statements.
  !
  !  We therefore want  BCOEF(J) to go into the entry of Q with index:
  !
  !    I -(LEFT+J)+2*K + ((LEFT+J)-K-1)*(2*K-1)
  !    =  begin_ligne +  (begin_col -1) * number_coef_different_0
  !   = I-LEFT+1+(LEFT -K)*(2*K-1) + (2*K-2)*J
  !
  jj = i - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
  
  do j = 1, k
     jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
     q(jj) = bcoef(j)
  end do

  bcoef_spline(i+l-1) = gtau(i)
  if ( tau_der(l) == i ) then   
     taui_der = taui
     call bsplvd( t, k, taui_der, left, a, bcoef_der, 2)
     l = l + 1
     jj = i - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
     do j = 1, k
        jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
        q(jj) = bcoef_der(j,2)
     end do
  bcoef_spline(i+ l-1) = gtau_der(l-1)
  end if

end do

taui = tau(n)
call interv( t, n+m+k, taui, left, mflag )
if ( tau_der(l)== n ) then   
  taui_der = taui
  call bsplvd( t, k, taui_der, left, a, bcoef_der, 2)
  jj = n - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
  do j = 1, k
     jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
     q(jj) = bcoef_der(j,2)
  end do
  bcoef_spline(n+l-1) = gtau_der(l)
  l = l + 1
end if

!  calculates the value of all possibly nonzero b-splines at taui of order
!  k  =  max( jhigh , (j+1)*(index-1) )
 
call bsplvb ( t, k, 1, taui, left, bcoef )
jj = n - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
   
do j = 1, k
  jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
  q(jj) = bcoef(j)
end do
bcoef_spline(n+l-1) = gtau(n)
!
!  Obtain factorization of A, stored again in Q.
!
call banfac ( q, k+k-1, n+m, k-1, k-1, iflag )

if ( iflag == 2 ) then
   write ( *, '(a)' ) ' '
   write ( *, '(a)' ) 'SPLINT - Fatal Error!'
   write ( *, '(a)' ) '  The linear system is not invertible!'
   stop
end if
!
!  Solve A * BCOEF = GTAU by back substitution.
!
call banslv ( q, k+k-1, n+m, k-1, k-1, bcoef_spline )

print*, 'PASSED'
end program unit_test
