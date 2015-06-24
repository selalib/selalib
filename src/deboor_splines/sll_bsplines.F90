module sll_bsplines

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_boundary_condition_descriptors.h"

use sll_module_deboor_splines_1d, only: bsplvb, interv

implicit none 

private
  

!> @brief 
!> basic type for one-dimensional B-spline data. 
!> @details This should be
!> treated as an opaque type. No access to its internals is directly allowed.
type, public :: sll_bspline_1d

  sll_int32           :: n
  sll_int32           :: k
  sll_real64, pointer :: tau(:)
  sll_real64, pointer :: t(:)
  sll_real64, pointer :: q(:)
  sll_real64, pointer :: bcoef(:)
  sll_real64, pointer :: bcoef_spline(:)
  sll_int32           :: bc_type
  sll_real64, pointer :: dbiatx(:,:)
  sll_real64, pointer :: a(:,:)
  sll_real64, pointer :: aj(:)
  sll_real64, pointer :: dl(:)
  sll_real64, pointer :: dr(:)

end type sll_bspline_1d

public :: initialize_bspline_1d
public :: compute_bspline_1d
public :: update_bspline_1d
public :: interpolate_array_values
public :: interpolate_array_derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize_bspline_1d(this, n, k, tau_min, tau_max, bc_type)

  type(sll_bspline_1d), intent(out) :: this
  sll_int32         , intent(in)    :: n
  sll_int32         , intent(in)    :: k
  sll_real64        , intent(in)    :: tau_min
  sll_real64        , intent(in)    :: tau_max
  sll_int32         , intent(in)    :: bc_type

  sll_int32                         :: i
  sll_int32                         :: ierr
  sll_int32, parameter              :: m=2

  this%n = n
  this%k = k

  SLL_ALLOCATE(this%tau(n), ierr)
  do i = 1, n
    this%tau(i) = tau_min + (i-1) * (tau_max-tau_min) / (n-1)
  end do

  SLL_ALLOCATE(this%t(n+k+m), ierr)

  this%t(1:k)         = tau_min
  this%t(k+1:n+m)     = this%tau(2:n-1)
  this%t(n+m+1:n+m+k) = tau_max

  SLL_ALLOCATE(this%q(1:(2*k-1)*(n+m)), ierr)
  SLL_ALLOCATE(this%bcoef(n),           ierr)
  SLL_ALLOCATE(this%bcoef_spline(n+m),  ierr)

  allocate(this%a(k,k))
  allocate(this%dbiatx(k,m))
  allocate(this%aj(k))
  allocate(this%aj(k))
  allocate(this%dl(k))
  allocate(this%dr(k))

end subroutine initialize_bspline_1d

subroutine build_system(this)

  type(sll_bspline_1d)    :: this 

  sll_real64              :: taui
  sll_int32               :: kpkm2
  sll_int32               :: left
  sll_int32               :: n
  sll_int32               :: k
  sll_int32               :: iflag
  sll_int32               :: mflag
  sll_int32               :: i
  sll_int32               :: j
  sll_int32               :: jj
  sll_int32               :: l
  sll_int32, parameter    :: m=2
  
  n = this%n
  k = this%k
  this%a = 0.0_f64

  kpkm2       = 2*(k-1)
  left        = k
  this%q      = 0.0_f64
  this%dbiatx = 0.0_f64
  
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
        call bsplvd( this%t, k, taui, left, this%a, this%dbiatx, 2)
        l = l + 1
        jj = i-left+1+(left-k)*(k+k-1)+l
        do j = 1, k
          jj = jj + kpkm2
          this%q(jj) = this%dbiatx(j,2)
        end do
      end if

    else

      call bsplvd( this%t, k, taui, left, this%a, this%dbiatx, 2)
      jj = i-left+1+(left-k)*(k+k-1)+l
      do j = 1, k
        jj = jj + kpkm2
        this%q(jj) = this%dbiatx(j,2)
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

subroutine compute_bspline_1d_aux( this, gtau, slope_min, slope_max)

  type(sll_bspline_1d)    :: this 
  sll_real64, intent(in)  :: gtau(:)
  sll_real64, optional    :: slope_min
  sll_real64, optional    :: slope_max

  sll_int32               :: n
  sll_int32               :: k
  sll_int32, parameter    :: m = 2

  n = this%n
  k = this%k

  SLL_ASSERT(size(gtau) == this%n)

  this%bcoef_spline(1)   = gtau(1)
  if (present(slope_min)) then
    this%bcoef_spline(2) = slope_min
  else
    this%bcoef_spline(2) = 0.0_f64
  end if
  this%bcoef_spline(3:n) = gtau(2:n-1)
  if (present(slope_max)) then
    this%bcoef_spline(n+1) = slope_max
  else
    this%bcoef_spline(n+1) = 0.0_f64
  end if
  this%bcoef_spline(n+2) = gtau(n)

  call banslv ( this%q, k+k-1, n+m, k-1, k-1, this%bcoef_spline )
  
end subroutine compute_bspline_1d_aux


!> @brief
!>  produces the B-spline coefficients of an interpolating spline.
!> @details
!>   The spline is of order K with knots T(1:N+K), and takes on the 
!>   value GTAU(I) at TAU(I), for I = 1 to N.
!>   The I-th equation of the linear system A * BCOEF = B 
!>   for the B-spline coefficients of the interpolant enforces interpolation
!>   at TAU(1:N).
!>   The matrix A is generated row by row and stored, diagonal by diagonal,
!>   in the rows of the array Q, with the main diagonal going
!>   into row K.  
!>   The banded system is then solved by a call to BANFAC, which 
!>   constructs the triangular factorization for A and stores it again in
!>   Q, followed by a call to BANSLV, which then obtains the solution
!>   BCOEF by substitution.
!>   The B-coefficients for the interpolant 
!>   of an additional data set can be obtained without going through all 
!>   the calculations in this routine, simply by calling the subroutine
!>   update_bspline_1d
!>
!>   @param[inout], the bspline object with knots and data point abscissas.
!>   The triangular factorization of the coefficient matrix of the linear 
!>   system for the B-coefficients is computed.
!>   @param[in]  gtau(N), the data ordinates.
!>   @param[in]  slope_min, the derivative at left boundary.
!>   @param[in]  slope_max, the derivative at right boundary.
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
subroutine compute_bspline_1d(this, gtau, slope_min, slope_max)

  type(sll_bspline_1d)    :: this 
  sll_real64, intent(in)  :: gtau(:)
  sll_real64, optional    :: slope_min
  sll_real64, optional    :: slope_max

  call build_system(this)

  call compute_bspline_1d_aux( this, gtau, slope_min, slope_max)

end subroutine compute_bspline_1d

!>@brief
!> produces the B-spline coefficients of an interpolating spline.
!> @details
!> If interpolants are already computed and if you change only 
!> the data set and knots and points positions did not change
!> use this routine.
subroutine update_bspline_1d(this, htau, slope_min, slope_max)

  type(sll_bspline_1d)    :: this 
  sll_real64, intent(in)  :: htau(:)
  sll_real64, optional    :: slope_min
  sll_real64, optional    :: slope_max

  call compute_bspline_1d_aux( this, htau, slope_min, slope_max)

end subroutine update_bspline_1d

!> @brief returns the values of the images of a collection of abscissae,
!> represented by a 1D array in another output array. The spline coefficients
!> used are stored in the spline object pointer.
!> @param[in] x input double-precison element array containing the 
!> abscissae to be interpolated.
!> @param[out] y output double-precision element array containing the 
!> results of the interpolation.
!> @param[in] n the number of elements of the input array which are to be
!> interpolated.
!> @param[inout] spline the spline object pointer, duly initialized and 
!> already operated on by the compute_bspline_1d() subroutine.
subroutine interpolate_array_values( this, n, x, y)

type(sll_bspline_1d)    :: this 
sll_int32,  intent(in)  :: n
sll_real64, intent(in)  :: x(n)
sll_real64, intent(out) :: y(n)

sll_int32               :: i
sll_int32               :: j
sll_int32               :: ilo
sll_int32               :: left
sll_int32               :: jc
sll_int32               :: jcmax
sll_int32               :: jcmin
sll_int32               :: jj
sll_int32               :: mflag
sll_int32               :: ierr
sll_int32               :: k
sll_int32               :: l
sll_int32               :: nmk
sll_int32,  parameter   :: m = 2
sll_real64, allocatable :: aj(:)
sll_real64, allocatable :: dl(:)
sll_real64, allocatable :: dr(:)
sll_real64              :: xi

k = this%k
nmk = n+m+k

SLL_CLEAR_ALLOCATE(aj(1:k), ierr)
SLL_CLEAR_ALLOCATE(dl(1:k), ierr)
SLL_CLEAR_ALLOCATE(dr(1:k), ierr)

do l = 1, n

  xi = x(l)

  call interv ( this%t, nmk, xi, i, mflag )
  
  if ( mflag /= 0 ) return

  if ( k <= 1 ) then
    y(l) = this%bcoef_spline(i)
    cycle
  end if
  
  if ( k <= i ) then
    do j = 1, k-1
      dl(j) = xi - this%t(i+1-j)
    end do
    jcmin = 1
  else
    jcmin = 1-(i-k)
    do j = 1, i
      dl(j) = xi - this%t(i+1-j)
    end do
    do j = i, k-1
      aj(k-j) = 0.0_f64
      dl(j) = dl(i)
    end do
  end if
  
  if ( n+m < i ) then
    jcmax = nmk-i
    do j = 1, jcmax
      dr(j) = this%t(i+j) - xi
    end do
    do j = jcmax, k-1
      aj(j+1) = 0.0_f64
      dr(j) = dr(jcmax)
    end do
  else
    jcmax = k
    do j = 1, k-1
      dr(j) = this%t(i+j) - xi
    end do
  end if
  
  do jc = jcmin, jcmax
    aj(jc) = this%bcoef_spline(i-k+jc)
  end do
  
  do j = 1, k-1
    ilo = k-j
    do jj = 1, k-j
      aj(jj) = (aj(jj+1)*dl(ilo)+aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
      ilo = ilo - 1
    end do
  end do
  
  y(l) = aj(1)

end do

DEALLOCATE(aj,dl,dr)

end subroutine interpolate_array_values

!> @brief returns the values of the derivatives evaluated at a collection of 
!> abscissae stored by a 1D array in another output array. The spline coefficients
!> used are stored in the spline object pointer.
!> @param[in] x input double-precison element array containing the 
!> abscissae to be interpolated.
!> @param[out] y output double-precision element array containing the 
!> results of the interpolation.
!> @param[in] n the number of elements of the input array which are to be
!> interpolated.
!> @param[inout] spline the spline object pointer, duly initialized and 
!> already operated on by the compute_bspline_1d() subroutine.
subroutine interpolate_array_derivatives( this, n, x, y)

type(sll_bspline_1d)    :: this 
sll_int32,  intent(in)  :: n
sll_real64, intent(in)  :: x(n)
sll_real64, intent(out) :: y(n)

sll_int32               :: i
sll_int32               :: j
sll_int32               :: ilo
sll_int32               :: left
sll_int32               :: jc
sll_int32               :: jcmax
sll_int32               :: jcmin
sll_int32               :: jj
sll_int32               :: mflag
sll_int32               :: ierr
sll_int32               :: k
sll_int32               :: l
sll_int32               :: nmk
sll_int32,  parameter   :: m = 2
sll_real64, allocatable :: aj(:)
sll_real64, allocatable :: dl(:)
sll_real64, allocatable :: dr(:)
sll_real64              :: xi
sll_int32               :: jderiv = 1

k = this%k
nmk = n+m+k

SLL_ALLOCATE(aj(k), ierr)
SLL_ALLOCATE(dl(k), ierr)
SLL_ALLOCATE(dr(k), ierr)

do l = 1, n

  xi = x(l)
  if ( k <= jderiv ) cycle
  !
  !  Find I so that 1 <= I < N+K and T(I) < T(I+1) and T(I) <= X < T(I+1).
  !
  !  If no such I can be found, X lies outside the support of the
  !  spline F and  BVALUE = 0.  The asymmetry in this choice of I makes F
  !  right continuous.
  !
  call interv ( this%t, n+m+k, xi, i, mflag )
  
  if ( mflag /= 0 ) return
  !
  !  If K = 1 (and JDERIV = 0), BVALUE = BCOEF(I).
  !
  if ( k <= 1 ) then
    y(l) = this%bcoef_spline(i)
    cycle
  end if
  !
  !  Store the K B-spline coefficients relevant for the knot interval
  !  ( T(I),T(I+1) ) in AJ(1),...,AJ(K) and compute DL(J) = X - T(I+1-J),
  !  DR(J) = T(I+J)-X, J=1,...,K-1.  Set any of the AJ not obtainable
  !  from input to zero.
  !
  !  Set any T's not obtainable equal to T(1) or to T(N+K) appropriately.
  !
  jcmin = 1
  
  if ( k <= i ) then
    do j = 1, k-1
      dl(j) = x(l) - this%t(i+1-j)
    end do
  else
    jcmin = 1 - ( i - k )
    do j = 1, i
      dl(j) = xi - this%t(i+1-j)
    end do
    do j = i, k-1
      aj(k-j) = 0.0_f64
      dl(j) = dl(i)
    end do
  end if
  
  jcmax = k
  if ( n + m < i ) then
    jcmax = k + n + m - i
    do j = 1, jcmax
      dr(j) = this%t(i+j) - xi
    end do
    do j = jcmax, k-1
      aj(j+1) = 0.0_8
      dr(j) = dr(jcmax)
    end do
  else
    do j = 1, k-1
      dr(j) = this%t(i+j) - xi
    end do
  end if
  
  do jc = jcmin, jcmax
    aj(jc) = this%bcoef_spline(i-k+jc)
  end do
  !  Difference the coefficients JDERIV times.
  do j = 1, jderiv
    ilo = k - j
    do jj = 1, k - j
      aj(jj) = ((aj(jj+1)-aj(jj))/(dl(ilo)+dr(jj)))*(k-j)
      ilo = ilo - 1
    end do
  end do
  !
  !  Compute value at X in (T(I),T(I+1)) of JDERIV-th derivative,
  !  given its relevant B-spline coefficients in AJ(1),...,AJ(K-JDERIV).
  !
  do j = jderiv+1, k-1
    ilo = k-j
    do jj = 1, k-j
      aj(jj) = (aj(jj+1)*dl(ilo)+aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
      ilo = ilo - 1
    end do
  end do
  
  y(l) = aj(1)

end do

DEALLOCATE(aj,dl,dr)

end subroutine interpolate_array_derivatives

!    Values at X of all the relevant B-splines of order K:K+1-NDERIV
!    are generated via BSPLVB and stored temporarily in DBIATX.
!
!    Then the B-spline coefficients of the required derivatives
!    of the B-splines of interest are generated by differencing,
!    each from the preceding one of lower order, and combined with
!    the values of B-splines of corresponding order in DBIATX
!    to produce the desired values.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!Input, real ( kind = 8 ) T(LEFT+K), the knot sequence.  It is assumed that
!    T(LEFT) < T(LEFT+1).  Also, the output is correct only if
!    T(LEFT) <= X <= T(LEFT+1).
!
!    Input, integer K, the order of the B-splines to be evaluated.
!
!    Input, real ( kind = 8 ) X, the point at which these values are sought.
!
!    Input, integer LEFT, indicates the left endpoint of the interval of
!    interest.  The K B-splines whose support contains the interval
!    ( T(LEFT), T(LEFT+1) ) are to be considered.
!
!    Workspace, real ( kind = 8 ) A(K,K).
!
!    Output, real ( kind = 8 ) DBIATX(K,NDERIV).  DBIATX(I,M) contains
!    the value of the (M-1)st derivative of the (LEFT-K+I)-th B-spline
!    of order K for knot sequence T, I=M,...,K, M=1,...,NDERIV.
!
!    Input, integer NDERIV, indicates that values of B-splines and their
!    derivatives up to but not including the NDERIV-th are asked for.
!
!> Calculates the nonvanishing B-splines and derivatives at X.
subroutine bsplvd ( t, k, x, left, a, dbiatx, nderiv )
    
sll_int32,  intent(in)    :: k
sll_real64, intent(in)    :: x
sll_int32,  intent(in)    :: left
sll_real64, intent(in)    :: t(left+k)
sll_real64, intent(inout) :: a(k,k)
sll_int32,  intent(in)    :: nderiv
sll_real64, intent(out)   :: dbiatx(k,nderiv)

sll_real64 :: factor
sll_real64 :: fkp1mm
sll_int32  :: i
sll_int32  :: ideriv
sll_int32  :: il
sll_int32  :: j
sll_int32  :: jlow
sll_int32  :: jp1mid
sll_int32  :: ldummy
sll_int32  :: m
sll_int32  :: mhigh

mhigh = max ( min ( nderiv, k ), 1 )
!
!  MHIGH is usually equal to NDERIV.
!
call bsplvb ( t, k+1-mhigh, 1, x, left, dbiatx )

if ( mhigh == 1 ) return
!
!  The first column of DBIATX always contains the B-spline values
!  for the current order.  These are stored in column K+1-current
!  order before BSPLVB is called to put values for the next
!  higher order on top of it.
!
ideriv = mhigh
do m = 2, mhigh
  jp1mid = 1
  do j = ideriv, k
     dbiatx(j,ideriv) = dbiatx(jp1mid,1)
     jp1mid = jp1mid + 1
  end do
  ideriv = ideriv - 1
  call bsplvb ( t, k+1-ideriv, 2, x, left, dbiatx )
end do
!
!  At this point, B(LEFT-K+I, K+1-J)(X) is in DBIATX(I,J) for
!  I=J,...,K and J=1,...,MHIGH ('=' NDERIV).
!
!  In particular, the first column of DBIATX is already in final form.
!
!  To obtain corresponding derivatives of B-splines in subsequent columns,
!  generate their B-representation by differencing, then evaluate at X.
!
jlow = 1
do i = 1, k
   a(jlow:k,i) = 0.0D+00
   jlow = i
   a(i,i) = 1.0D+00
end do
!
!  At this point, A(.,J) contains the B-coefficients for the J-th of the
!  K B-splines of interest here.
!
do m = 2, mhigh
  fkp1mm = real ( k + 1 - m, kind = 8 )
  il = left
  i = k
  !
  !  For J = 1,...,K, construct B-coefficients of (M-1)st derivative of
  !  B-splines from those for preceding derivative by differencing
  !  and store again in  A(.,J).  The fact that  A(I,J) = 0 for
  !  I < J is used.
  !
  do ldummy = 1, k+1-m
     factor = fkp1mm / ( t(il+k+1-m) - t(il) )
     !  The assumption that T(LEFT) < T(LEFT+1) makes denominator
     !  in FACTOR nonzero.
     a(i,1:i) = ( a(i,1:i) - a(i-1,1:i) ) * factor
     il = il - 1
     i = i - 1
  end do
  !  For I = 1,...,K, combine B-coefficients A(.,I) with B-spline values
  !  stored in DBIATX(.,M) to get value of (M-1)st derivative of
  !  I-th B-spline (of interest here) at X, and store in DBIATX(I,M).
  !
  !  Storage of this value over the value of a B-spline
  !  of order M there is safe since the remaining B-spline derivatives
  !  of the same order do not use this value due to the fact
  !  that  A(J,I) = 0  for J < I.
  do i = 1, k
     jlow = max ( i, m )
     dbiatx(i,m) = dot_product ( a(jlow:k,i), dbiatx(jlow:k,m) )
  end do
end do
end subroutine bsplvd
  

end module sll_bsplines
