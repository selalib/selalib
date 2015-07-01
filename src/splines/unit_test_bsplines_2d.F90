program test_bsplines_2d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_utilities.h"
use sll_bsplines
use sll_boundary_condition_descriptors

  implicit none

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

  sll_real64           :: gtau(nx,ny)
  sll_real64           :: bcoef(nx,ny)
  sll_real64           :: taux(nx)
  sll_real64           :: tauy(ny)
  sll_real64           :: tx(nx+kx)
  sll_real64           :: ty(ny+ky)
  sll_real64           :: work1(ny,nx)
  sll_real64           :: work2(nx)
  sll_real64           :: work3(nx*ny)
  sll_real64           :: workx(nx)
  sll_real64           :: worky(ny)

  type(sll_bspline_2d), pointer :: bspline_2d
  
  
  bspline_2d => new_bspline_2d( nx, kx-1, 1.0_f64, nx*1.0_f64, SLL_PERIODIC, &
                                ny, ky-1, 1.0_f64, ny*1.0_f64, SLL_PERIODIC  )

  ! set up data points and knots
  ! in x, interpolate between knots by parabolic splines, using
  ! not-a-knot end condition
  taux = bspline_2d%bs1%tau
  tx   = bspline_2d%bs1%t
  tauy = bspline_2d%bs2%tau
  ty   = bspline_2d%bs2%t
  
  
  ! generate and print out function values
  print 620,(tauy(i),i=1,ny)
  do i=1,nx
    do j=1,ny
      bcoef(i,j) = g(taux(i),tauy(j))
    end do
    print 632,taux(i),(bcoef(i,j),j=1,ny)
  end do
  
  call build_system( bspline_2d%bs1 )  
  call build_system( bspline_2d%bs2 )
  call update_bspline_2d( bspline_2d, bcoef) 
  bcoef = bspline_2d%bcoef
 
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
