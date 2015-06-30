program test_bsplines_2d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_utilities.h"
use sll_bsplines

  sll_int32, parameter :: nx=7
  sll_int32, parameter :: kx=3
  sll_int32, parameter :: ny=6
  sll_int32, parameter :: ky=4

  sll_int32            :: i
  sll_int32            :: j
  sll_int32            :: iflag
  sll_int32            :: jj
  sll_int32            :: lefty
  sll_int32            :: mflag
  sll_int32            :: ilo

  sll_real64            :: gtau(nx,ny)
  sll_real64            :: bcoef(nx,ny)
  sll_real64            :: taux(nx)
  sll_real64            :: tauy(ny)
  sll_real64            :: tx(nx+kx)
  sll_real64            :: ty(ny+ky)
  sll_real64            :: work1(ny,nx)
  sll_real64            :: work2(nx)
  sll_real64            :: work3(nx*ny)
  
  
  ! set up data points and knots
  ! in x, interpolate between knots by parabolic splines, using
  ! not-a-knot end condition
  
  do i=1,nx
    taux(i) = float(i)
  end do
  do i=1,kx
    tx(i)    = taux(1)
    tx(nx+i) = taux(nx)
  end do
  kp1 = kx+1  
  do i=kp1,nx
    tx(i) = (taux(i-kx+1) + taux(i-kx+2))/2.
  end do
  
  ! in y, interpolate at knots by cubic splines, using not-a-knot
  ! end condition
  
  do i=1,ny
    tauy(i) = float(i)
  end do
  do i=1,ky
    ty(i) = tauy(1)
    ty(ny+i) = tauy(ny)
  end do
  kp1 = ky+1  
  do i=kp1,ny
    ty(i) = tauy(i-ky+2)
  end do
  
  ! generate and print out function values
  print 620,(tauy(i),i=1,ny)
  do i=1,nx
    do j=1,ny
      bcoef(i,j) = g(taux(i),tauy(j))
    end do
    print 632,taux(i),(bcoef(i,j),j=1,ny)
  end do
  
  ! construct b-coefficients of interpolant
  call spli2d(taux,bcoef,tx,nx,kx,ny,work2,work3,work1,iflag)
  call spli2d(tauy,work1,ty,ny,ky,nx,work2,work3,bcoef,iflag)
  
  ! evaluate interpolation error at mesh points and print out
  print 640,(tauy(j),j=1,ny)
  ilo = ky
  do j=1,ny
    call interv(ty,ny+1,tauy(j),lefty,mflag)
    do i=1,nx
      do jj=1,ky
        work2(jj)=bvalue(tx,bcoef(:,lefty-ky+jj),nx,kx,taux(i),0)
      end do
      gtau(i,j) = g(taux(i),tauy(j)) - &
                   bvalue(ty(lefty-ky+1:),work2,ky,ky,tauy(j),0)
    end do
  end do
  do i=1,nx
    print 632,taux(i),(gtau(i,j),j=1,ny)
  end do

  print*, 'PASSED'

  620 format(' given data'//6f13.1)
  632 format(f5.1,6e13.5)
  640 format(//' interpolation error'//6f13.1)

contains

!     
!> @param[in]  tau array of length  n , containing data point abscissae.  
!>   (strictly  increasing)
!> @param[in]  gtau array of length  n , containing data point ordinates, j=1,...,m
!> @param[in]  t knot sequence, of length  n+k
!> @param[in]  n number of data points and dimension of spline space s(k,t)
!> @param[in]  k order of spline
!> @param[in]  m number of data sets
!> @param[in]  work a vector of length n
!> @param[out] q array of size  (2*k-1)*n , containing the triangular factoriz- 
!>               ation of the coefficient matrix of the linear system for the b-
!>               coefficients of the spline interpolant.
!>               the b-coeffs for the interpolant of an additional data set  
!>               (tau(i),htau(i)), i=1,...,n  with the same data abscissae can  
!>               be obtained without going through all the calculations in this 
!>               routine, simply by loading  htau  into  bcoef  and then execut-
!>       ing the    call banslv ( q, 2*k-1, n, k-1, k-1, bcoef )  
!> @param[out] bcoef the b-coefficients of the interpolant, of length  n  
!> @param[out] iflag an sll_int32 indicating success (= 1)  or failure (= 2)
!>       the linear system to be solved is (theoretically) invertible if
!>       and only if    
!>             t(i) .lt. tau(i) .lt. tau(i+k),    all i.    
!>       violation of this condition is certain to lead to  iflag = 2 . 
!>    
!>@brief Produces the b-spline coeff.s  bcoef(j,.)  of the spline of
!>  order  k  with knots  t (i), i=1,..., n + k , which takes on the
!>  value  gtau (i,j)  at  tau (i), i=1,..., n , j=1,..., m .
!> @details
!>    the i-th equation of the linear system  a*bcoef = b  for the b-co-
!> effs of the interpolant enforces interpolation at  tau(i), i=1,...,n.
!> hence,  b(i) = gtau(i), all i, and  a  is a band matrix with  2k-1   
!>  bands (if it is invertible).    
!>    the matrix  a  is generated row by row and stored, diagonal by di-
!> agonal, in the  rows  of the array  q , with the main diagonal go-
!> ing into row  k .  see comments in the program below.    
!>    the banded system is then solved by a call to  banfac (which con- 
!> structs the triangular factorization for  a  and stores it again in  
!>  q ), followed by a call to  banslv (which then obtains the solution 
!>  bcoef  by substitution).  
!>    banfac  does no pivoting, since the total positivity of the matrix
!> a  makes this unnecessary. 
!>    
subroutine spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )

sll_real64, intent(in)  :: tau(n)
sll_real64, intent(in)  :: gtau(n,m)
sll_real64, intent(in)  :: t(n+k)
sll_int32,  intent(in)  :: n
sll_int32,  intent(in)  :: k 
sll_int32,  intent(in)  :: m
sll_real64, intent(out) :: work(n)
sll_real64, intent(out) :: q((2*k-1)*n)
sll_real64, intent(out) :: bcoef(m,n)

sll_int32 :: km1
sll_int32 :: kpkm2
sll_int32 :: left
sll_int32 :: np1 
sll_int32 :: jj
sll_int32 :: j
sll_int32 :: i
sll_int32 :: ilp1mx
sll_int32 :: iflag

sll_real64 :: taui 

np1   = n + 1 
km1   = k - 1 
kpkm2 = 2*km1     
left  = k    
q     = 0.0_f64  ! zero out all entries of q
     
! loop over i to construct the  n  interpolation equations 
do i=1,n 

  taui = tau(i)  
  ilp1mx = min0(i+k,np1)     
  ! find  left  in the closed interval (i,i+k-1) such that     
  ! t(left) .le. tau(i) .lt. t(left+1)   
  ! matrix is singular if this is not possible   
  left = max0(left,i)  
  if (taui .lt. t(left))      stop ' linear system not invertible '
  15 if (taui .lt. t(left+1)) go to 16
  left = left + 1   
  if (left .lt. ilp1mx)       go to 15
  left = left - 1
  if (taui .gt. t(left+1))    stop ' linear system not invertible '

  ! *** the i-th equation enforces interpolation at taui, hence    
  ! a(i,j) = b(j,k,t)(taui), all j. only the  k  entries with  j = 
  ! left-k+1,...,left actually might be nonzero. these  k  numbers 
  ! are returned, in  work  (used for temp.storage here), by the   
  ! following
  16 call bsplvb ( t, k, 1, taui, left, work )                     
  ! we therefore want  work(j) = b(left -k+j)(taui) to go into     
  ! a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since  
  ! a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q 
  ! as a two-dim. array , with  2*k-1  rows (see comments in 
  ! banfac). in the present program, we treat  q  as an equivalent 
  ! one-dimensional array (because of fortran restrictions on
  ! dimension statements) . we therefore want  work(j) to go into  
  ! entry    
  !     i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)   
  !            =  i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j   
  ! of  q .  
  jj = i-left+1 + (left-k)*(k+km1) 
  do j=1,k    
    jj = jj+kpkm2     
    q(jj) = work(j)
  end do
end do
     
! ***obtain factorization of  a  , stored again in  q.  

call banfac ( q, k+km1, n, km1, km1, iflag )    

if (iflag == 1) then  ! solve  a*bcoef = gtau  by backsubstitution  
  do j=1,m
    do i=1,n
      work(i) = gtau(i,j)
    end do
    call banslv ( q, k+km1, n, km1, km1, work ) 
    do i=1,n                                   
      bcoef(j,i) = work(i)
    end do
  end do
end if   

end subroutine spli2d 

function g (x , y)

  sll_real64 :: x
  sll_real64 :: y
  sll_real64 :: g
  
  g = max(x-3.5_f64,0.0_f64)**2 + max(y-3.0_f64,0.0_f64)**3

end function g

subroutine interv( xt, lxt, x, left, mflag )
 
sll_real64, intent(in)  :: xt(:)
sll_int32,  intent(in)  :: lxt
sll_real64, intent(in)  :: x
sll_int32,  intent(out) :: left
sll_int32,  intent(out) :: mflag

sll_int32               :: ihi
sll_int32, save         :: ilo = 1
sll_int32               :: istep
sll_int32               :: middle

ihi = ilo + 1
if ( lxt <= ihi ) then
   if ( xt(lxt) <= x ) then
      go to 110
   end if
   if ( lxt <= 1 ) then
      mflag = -1
      left = 1
      return
   end if
   ilo = lxt - 1
   ihi = lxt
end if
if ( xt(ihi) <= x ) then
   go to 20
end if
if ( xt(ilo) <= x ) then
   mflag = 0
   left = ilo
   return
end if
!
!  Now X < XT(ILO).  Decrease ILO to capture X.
!
istep = 1
10  continue
ihi = ilo
ilo = ihi - istep
if ( 1 < ilo ) then
   if ( xt(ilo) <= x ) then
      go to 50
   end if
   istep = istep * 2
   go to 10
end if
ilo = 1
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
ilo = ihi
ihi = ilo + istep
if ( ihi < lxt ) then
   if ( x < xt(ihi) ) then
      go to 50
   end if
   istep = istep * 2
   go to 30
end if
if ( xt(lxt) <= x ) then
   go to 110
end if
!
!  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
!
ihi = lxt
50  continue
do
  middle = ( ilo + ihi ) / 2
  if ( middle == ilo ) then
     mflag = 0
     left = ilo
     return
  end if
  !
  !  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
  !
  if ( xt(middle) <= x ) then
     ilo = middle
  else
     ihi = middle
  end if
   
end do
!
!  Set output and return.
!
110 continue
mflag = 1
if ( x == xt(lxt) ) then
   mflag = 0
end if
do left = lxt, 1, -1
   if ( xt(left) < xt(lxt) ) then
      return
   end if
end do

end subroutine interv

!> Evaluates B-splines at a point X with a given knot sequence.
subroutine bsplvb ( t, jhigh, index, x, left, biatx )
  
sll_int32,  intent(in)  :: jhigh
sll_int32,  intent(in)  :: index
sll_real64, intent(in)  :: x
sll_int32,  intent(in)  :: left
sll_real64, intent(in)  :: t(left+jhigh)
sll_real64, intent(out) :: biatx (jhigh)

sll_int32,  parameter             :: jmax = 20
sll_real64, save, dimension(jmax) :: deltal
sll_real64, save, dimension(jmax) :: deltar
sll_int32                         :: i
sll_int32, save                   :: j = 1
sll_real64                        :: saved
sll_real64                        :: term


if ( index == 1 ) then
  j = 1
  biatx(1) = 1.0_f64
  if ( jhigh <= j ) then
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
   
  deltar(j) = t(left+j) - x
  deltal(j) = x - t(left+1-j)
  
  saved = 0.0_f64
  do i = 1, j
     term = biatx(i) / ( deltar(i) + deltal(j+1-i) )
     biatx(i) = saved + deltar(i) * term
     saved = deltal(j+1-i) * term
  end do

  biatx(j+1) = saved
  j = j + 1
  
  if ( jhigh <= j ) exit

end do

end subroutine bsplvb
  
!> Evaluates a derivative of a spline from its B-spline representation.
function bvalue( t, bcoef, n, k, x, jderiv ) result(res)
    
sll_real64, intent(in)  :: t(:)
sll_real64, intent(in)  :: bcoef(:)
sll_int32,  intent(in)  :: n
sll_int32,  intent(in)  :: k
sll_real64, intent(in)  :: x
sll_int32,  intent(in)  :: jderiv

sll_real64              :: res

sll_real64, allocatable :: aj(:)
sll_real64, allocatable :: dl(:)
sll_real64, allocatable :: dr(:)

sll_int32 :: i
sll_int32 :: ilo
sll_int32 :: j
sll_int32 :: jc
sll_int32 :: jcmax
sll_int32 :: jcmin
sll_int32 :: jj
sll_int32 :: mflag
sll_int32 :: ierr

res = 0.0_f64

SLL_CLEAR_ALLOCATE(aj(1:k), ierr)
SLL_CLEAR_ALLOCATE(dl(1:k), ierr)
SLL_CLEAR_ALLOCATE(dr(1:k), ierr)

if ( k <= jderiv ) return
!
!  Find I so that 1 <= I < N+K and T(I) < T(I+1) and T(I) <= X < T(I+1).
!
!  If no such I can be found, X lies outside the support of the
!  spline F and  BVALUE = 0.  The asymmetry in this choice of I makes F
!  right continuous.
!
call interv ( t, n+k, x, i, mflag )

if ( mflag /= 0 ) return
!
!  If K = 1 (and JDERIV = 0), BVALUE = BCOEF(I).
!
if ( k <= 1 ) then
  res = bcoef(i)
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
    aj(k-j) = 0.0_f64
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
    aj(j+1) = 0.0_f64
    dr(j) = dr(k+n-i)
  end do
else
  do j = 1, k-1
    dr(j) = t(i+j) - x
  end do
end if

do jc = jcmin, jcmax
  aj(jc) = bcoef(i-k+jc)
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

res = aj(1)

end function bvalue
  
end program test_bsplines_2d
