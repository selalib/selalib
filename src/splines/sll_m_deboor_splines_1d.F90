module sll_m_deboor_splines_1d

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

use sll_m_boundary_condition_descriptors

implicit none 
  
!this derived type is created to avoid the save attribute you found
!in module deboor splines 1d
!Some deboor functions are clone copied because we don't want to have
!side effects with general coordinates elliptic solver that uses
!module deboor splines
type, public :: deboor_type
   sll_int32 :: ilo = 1
   sll_int32 :: j   = 1
   sll_real64, dimension(20) :: deltal
   sll_real64, dimension(20) :: deltar
end type deboor_type


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bsplvd ( deboor, t, k, x, left, a, dbiatx, nderiv )

type(deboor_type) :: deboor
sll_int32 :: k
sll_int32 :: left
sll_int32 :: nderiv

sll_real64, dimension(k,k):: a!(k,k)
sll_real64, dimension(k,nderiv), intent(out) :: dbiatx!(k,nderiv)
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
!  sll_real64 sum1  ! this one is not used...
sll_real64,dimension(left+k):: t ! (left+k)
sll_real64:: x

mhigh = max ( min ( nderiv, k ), 1 )
!
!  MHIGH is usually equal to NDERIV.
!
call bsplvb ( deboor, t, k+1-mhigh, 1, x, left, dbiatx )

if ( mhigh == 1 ) then
   return
end if
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
   
   call bsplvb ( deboor, t, k+1-ideriv, 2, x, left, dbiatx )
   
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
      !
      !  The assumption that T(LEFT) < T(LEFT+1) makes denominator
      !  in FACTOR nonzero.
      !
      a(i,1:i) = ( a(i,1:i) - a(i-1,1:i) ) * factor
      
      il = il - 1
      i = i - 1
   
   end do
   !
   !  For I = 1,...,K, combine B-coefficients A(.,I) with B-spline values
   !  stored in DBIATX(.,M) to get value of (M-1)st derivative of
   !  I-th B-spline (of interest here) at X, and store in DBIATX(I,M).
   !
   !  Storage of this value over the value of a B-spline
   !  of order M there is safe since the remaining B-spline derivatives
   !  of the same order do not use this value due to the fact
   !  that  A(J,I) = 0  for J < I.
   !
   do i = 1, k
      
      jlow = max ( i, m )
      
      dbiatx(i,m) = dot_product ( a(jlow:k,i), dbiatx(jlow:k,m) )
      
   end do

end do
return
end subroutine bsplvd

subroutine interv( deboor, xt, lxt, x, left, mflag )
    
type(deboor_type)                     :: deboor
sll_real64, dimension(:), intent(in)  :: xt
sll_int32,                intent(in)  :: lxt
sll_real64,               intent(in)  :: x
sll_int32,                intent(out) :: left
sll_int32,                intent(out) :: mflag
sll_int32                             :: ihi

sll_int32 :: istep
sll_int32 :: middle

SLL_ASSERT(size(xt) == lxt)

ihi = deboor%ilo + 1
if ( lxt <= ihi ) then
   if ( xt(lxt) <= x ) go to 110
   if ( lxt <= 1 ) then
      mflag = -1
      left = 1
      return
   end if
   deboor%ilo = lxt - 1
   ihi = lxt
end if
if ( xt(ihi) <= x ) go to 20
if ( xt(deboor%ilo) <= x ) then
   mflag = 0
   left = deboor%ilo
   return
end if
!
!  Now X < XT(ILO).  Decrease ILO to capture X.
!
istep = 1
    
10  continue
    
ihi = deboor%ilo
deboor%ilo = ihi - istep

if ( 1 < deboor%ilo ) then
   if ( xt(deboor%ilo) <= x ) go to 50
   istep = istep * 2
   go to 10
end if

deboor%ilo = 1

if ( x < xt(1) ) then
   mflag = -1
   left = 1
   return
end if

go to 50
!
!  Now XT(IHI) <= X.  Increase IHI to capture X.
!
20  continue
    
istep = 1
    
30  continue
    
deboor%ilo = ihi
ihi = deboor%ilo + istep

if ( ihi < lxt ) then
   if ( x < xt(ihi) ) go to 50
   istep = istep * 2
   go to 30
end if

if ( xt(lxt) <= x ) go to 110
!
!  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
!
ihi = lxt

50 continue
    
do
   
  middle = ( deboor%ilo + ihi ) / 2
  if ( middle == deboor%ilo ) then
    mflag = 0
    left = deboor%ilo
    return
  end if
  !
  !  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
  !
  if ( xt(middle) <= x ) then
    deboor%ilo = middle
  else
    ihi = middle
  end if
   
end do
!
!  Set output and return.
!

110 continue
    
mflag = 1

if ( x == xt(lxt) ) mflag = 0

do left = lxt, 1, -1
  if ( xt(left) < xt(lxt) ) return
end do

end subroutine interv

subroutine bsplvb ( deboor, t, jhigh, index, x, left, biatx )

type(deboor_type) :: deboor
sll_int32:: jhigh

sll_real64,dimension(jhigh):: biatx !(jhigh)
sll_int32:: i
sll_int32:: index
sll_int32:: left
sll_real64:: saved
sll_real64,dimension(left+jhigh), intent(in) :: t!() left+jhigh
sll_real64:: term
sll_real64:: x

if ( index == 1 ) then
   deboor%j = 1
   biatx(1) = 1.0_8
   if ( jhigh <= deboor%j ) then
      return
   end if
end if

if ( t(left+1) <= t(left) ) then
   print*,'x=',x
   write ( *, '(a)' ) ' '
   write ( *, '(a)' ) 'BSPLVB - Fatal error!'
   write ( *, '(a)' ) '  It is required that T(LEFT) < T(LEFT+1).'
   write ( *, '(a,i8)' ) '  But LEFT = ', left
   write ( *, '(a,g14.6)' ) '  T(LEFT) =   ', t(left)
   write ( *, '(a,g14.6)' ) '  T(LEFT+1) = ', t(left+1)
   stop
end if

do
   
  deboor%deltar(deboor%j) = t(left+deboor%j) - x
  deboor%deltal(deboor%j) = x - t(left+1-deboor%j)
   
  saved = 0.0_f64
  do i = 1, deboor%j
     term = biatx(i) / ( deboor%deltar(i) + deboor%deltal(deboor%j+1-i) )
     biatx(i) = saved + deboor%deltar(i) * term
     saved = deboor%deltal(deboor%j+1-i) * term
  end do

  biatx(deboor%j+1) = saved
  deboor%j = deboor%j + 1
  
  if ( jhigh <= deboor%j ) exit

end do

return
end subroutine bsplvb
  
function bvalue( deboor, t, coeff_splines, n, k, x, jderiv )
    
implicit none

type(deboor_type)        :: deboor
sll_real64               :: bvalue

sll_real64, dimension(:) :: t     
sll_real64, dimension(:) :: coeff_splines 
sll_int32                :: n
sll_int32                :: k
sll_real64               :: x
sll_int32                :: jderiv

sll_real64, dimension(k) :: aj
sll_real64, dimension(k) :: dl
sll_real64, dimension(k) :: dr

sll_int32 :: i
sll_int32 :: ilo
sll_int32 :: j
sll_int32 :: jc
sll_int32 :: jcmax
sll_int32 :: jcmin
sll_int32 :: jj
sll_int32 :: mflag

bvalue = 0.0_8

aj(:)=0.0_8
dl(:)=0.0_8
dr(:)=0.0_8

if ( k <= jderiv ) return

!  Find I so that 1 <= I < N+K and T(I) < T(I+1) and T(I) <= X < T(I+1).
!
!  If no such I can be found, X lies outside the support of the
!  spline F and  BVALUE = 0.  The asymmetry in this choice of I makes F
!  right continuous.

call interv ( deboor, t, n+k, x, i, mflag )

if ( mflag /= 0 ) return
!
!  If K = 1 (and JDERIV = 0), BVALUE = BCOEF(I).
!
if ( k <= 1 ) then
  bvalue = coeff_splines(i)
  return
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
    dl(j) = x - t(i+1-j)
  end do
   
else
   
  jcmin = 1 - ( i - k )
  do j = 1, i
    dl(j) = x - t(i+1-j)
  end do
  do j = i, k-1
    aj(k-j) = 0.0_8
    dl(j) = dl(i)
  end do
   
end if

jcmax = k

if ( n < i ) then
   
  jcmax = k + n - i
  do j = 1, k + n - i
    dr(j) = t(i+j) - x
  end do
   
  do j = k+n-i, k-1
    aj(j+1) = 0.0_8
    dr(j) = dr(k+n-i)
  end do
   
else
   
  do j = 1, k-1
    dr(j) = t(i+j) - x
  end do
   
end if

do jc = jcmin, jcmax
  aj(jc) = coeff_splines(i-k+jc)
end do
!
!  Difference the coefficients JDERIV times.
!
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

bvalue = aj(1)

return

end function bvalue
  
subroutine splint ( deboor, tau, gtau, t, n, k, q, coeff_splines, iflag )
implicit none

type(deboor_type)                     :: deboor
sll_real64, dimension(:), intent(in)  :: tau
sll_real64, dimension(:), intent(in)  :: gtau 
sll_real64, dimension(:), intent(in)  :: t
sll_int32,                intent(in)  :: n
sll_int32,                intent(in)  :: k
sll_real64, dimension(:), intent(out) :: q
sll_real64, dimension(:), intent(out) :: coeff_splines
sll_int32,                intent(out) :: iflag

sll_int32  :: j
sll_int32  :: jj
sll_int32  :: kpkm2
sll_int32  :: left
sll_int32  :: ilp1mx
sll_real64 :: taui
sll_int32  :: i

SLL_ASSERT(size(tau)  == n)
SLL_ASSERT(size(gtau) == n)
SLL_ASSERT(size(t)    == n+k)
SLL_ASSERT(size(q)    == (2*k-1)*n)

kpkm2 = 2 * ( k - 1 )
left  = k
q     = 0.0_f64

!  Loop over I to construct the N interpolation equations.
do i = 1, n
   
  taui = tau(i)
  ilp1mx = min ( i + k, n + 1 )
  !
  !  Find LEFT in the closed interval (I,I+K-1) such that
  !
  !    T(LEFT) <= TAU(I) < T(LEFT+1)
  !
  !  The matrix is singular if this is not possible.
  !
  left = max ( left, i )
  
  if ( taui < t(left) ) then
    iflag = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINT - Fatal Error!'
    write ( *, '(a)' ) '  The linear system is not invertible!'
    return
  end if

  do while ( t(left+1) <= taui )
     
    left = left + 1
     
    if ( left < ilp1mx ) then
      cycle
    end if
     
    left = left - 1
     
    if ( t(left+1) < taui ) then
      iflag = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINT - Fatal Error!'
      write ( *, '(a)' ) '  The linear system is not invertible!'
      return
    end if
     
    exit
     
  end do
  !
  !  The I-th equation enforces interpolation at TAUI, hence for all J,
  !
  !    A(I,J) = B(J,K,T)(TAUI).
  !
  !Only the K entries with J = LEFT-K+1,...,LEFT actually might be nonzero.
  !
  !These K numbers are returned, in BCOEF 
  ! (used for temporary storage here),
  !  by the following.
  !
  call bsplvb ( deboor, t, k, 1, taui, left, coeff_splines )
  !
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
  !   = I-LEFT+1+(LEFT -K)*(2*K-1) + (2*K-2)*J
  !
  jj = i - left + 1 + ( left - k ) * ( k + k - 1 )
  
  do j = 1, k
    jj = jj + kpkm2
    q(jj) = coeff_splines(j)
  end do
  
end do
!
!  Obtain factorization of A, stored again in Q.
!
call banfac ( q, k+k-1, n, k-1, k-1, iflag )

if ( iflag == 2 ) then
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPLINT - Fatal Error!'
  write ( *, '(a)' ) '  The linear system is not invertible!'
  return
end if
!
!  Solve 
!
!    A * BCOEF = GTAU
!
!  by back substitution.
!
coeff_splines(1:n) = gtau(1:n)

call banslv ( q, k+k-1, n, k-1, k-1, coeff_splines )

return
end subroutine splint

subroutine splint_der( deboor,       &
                       tau,          &
                       gtau,         &
                       tau_der,      &
                       gtau_der,     &
                       t,            &
                       n,            &
                       m,            &
                       k,            &
                       q,            &
                       bcoef_spline, &
                       iflag )
    
type(deboor_type)                     :: deboor
sll_real64, dimension(:), intent(in)  :: tau
sll_real64, dimension(:), intent(in)  :: gtau
sll_int32,  dimension(:), intent(in)  :: tau_der
sll_real64, dimension(:), intent(in)  :: gtau_der 
sll_real64, dimension(:), intent(in)  :: t
sll_int32,                intent(in)  :: n
sll_int32,                intent(in)  :: m
sll_real64, dimension(:), intent(out) :: q
sll_real64, dimension(:), intent(out) :: bcoef_spline 
sll_int32,                intent(out) :: iflag


sll_real64,dimension(n+m)  :: bcoef
sll_int32                  :: i
sll_int32                  :: j
sll_int32                  :: jj
sll_int32                  :: k
sll_int32                  :: l
sll_int32                  :: mflag
sll_int32                  :: kpkm2
sll_int32                  :: left
sll_real64                 :: taui
sll_real64                 :: taui_der
sll_real64, dimension(k,k) :: a
sll_real64, dimension(k,2) :: bcoef_der

kpkm2     = 2*(k-1)
left      = k
q         = 0.0_f64
a         = 0.0_f64
bcoef_der = 0.0_f64

! we must suppose that m is <= than n 
SLL_ASSERT(m <= n)
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

  call interv( deboor, t, n+m+k, taui, left, mflag )

  !
  !  The I-th equation enforces interpolation at TAUI, hence for all J,
  !
  !    A(I,J) = B(J,K,T)(TAUI).
  !
  !Only the K entries with J = LEFT-K+1,...,LEFT actually might be nonzero.
  !
  !These K numbers are returned, in BCOEF 
  ! (used for temporary storage here),
  !  by the following.
  
  call bsplvb ( deboor, t, k, 1, taui, left, bcoef )
   
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

  bcoef_spline(i+ l-1) = gtau(i)
  if ( tau_der(l) == i ) then   
    taui_der = taui
      
    call bsplvd( deboor, t, k, taui_der, left, a, bcoef_der, 2)

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
call interv( deboor, t, n+m+k, taui, left, mflag )
if ( tau_der(l)== n ) then   
  taui_der = taui
      
  call bsplvd( deboor, t, k, taui_der, left, a, bcoef_der, 2)
      
  jj = n - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
   
  do j = 1, k
    jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
    q(jj) = bcoef_der(j,2)
  end do
  bcoef_spline(n+ l-1) = gtau_der(l)
  l = l + 1
      
end if

call bsplvb ( deboor, t, k, 1, taui, left, bcoef )
jj = n - left + 1 + ( left - k ) * ( k + k - 1 ) + l - 1
   
do j = 1, k
  jj = jj + kpkm2  ! kpkm2 = 2*(k-1)
  q(jj) = bcoef(j)
end do
bcoef_spline(n+l-1) = gtau(n)

!  Obtain factorization of A, stored again in Q.
!
call banfac ( q, k+k-1, n+m, k-1, k-1, iflag )

if ( iflag == 2 ) then
   write ( *, '(a)' ) ' '
   write ( *, '(a)' ) 'SPLINT - Fatal Error!'
   write ( *, '(a)' ) '  The linear system is not invertible!'
   return
end if
!
!  Solve 
!
!    A * BCOEF = GTAU
!
!  by back substitution.

call banslv ( q, k+k-1, n+m, k-1, k-1, bcoef_spline )

return
end subroutine splint_der


end module sll_m_deboor_splines_1d
