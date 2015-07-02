program test_bsplines_2d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_utilities.h"
use sll_bsplines
use sll_boundary_condition_descriptors

  implicit none

  sll_int32, parameter :: nx=17
  sll_int32, parameter :: kx=3
  sll_int32, parameter :: ny=11
  sll_int32, parameter :: ky=4

  sll_int32            :: i
  sll_int32            :: j
  sll_int32            :: jj
  sll_int32            :: left
  sll_int32            :: leftx
  sll_int32            :: lefty
  sll_int32            :: mflag
  sll_int32            :: ilo
  sll_int32            :: jlo
  sll_int32            :: klo
  sll_int32            :: ierr

  sll_real64           :: gtau(nx,ny)
  sll_real64           :: bcoef(nx,ny)
  sll_real64           :: taux(nx)
  sll_real64           :: tauy(ny)
  sll_real64           :: tx(nx+kx)
  sll_real64           :: ty(ny+ky)
  sll_real64           :: work(nx)

  sll_real64, allocatable :: aj(:)
  sll_real64, allocatable :: dl(:)
  sll_real64, allocatable :: dr(:)

  type(sll_bspline_2d), pointer :: bspline_2d

  sll_real64 :: t0, t1, t2
  
  
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
  
  call cpu_time(t0)
  call compute_bspline_2d( bspline_2d, bcoef) 
  bcoef = bspline_2d%bcoef
  call cpu_time(t1)
 
  ! evaluate interpolation error at mesh points and print out
  print 640,(tauy(j),j=1,ny)

  SLL_CLEAR_ALLOCATE(aj(1:max(kx,ky)), ierr)
  SLL_CLEAR_ALLOCATE(dl(1:max(kx,ky)), ierr)
  SLL_CLEAR_ALLOCATE(dr(1:max(kx,ky)), ierr)

  jlo = ky
  do j=1,ny
    call interv(ty,ny+1,tauy(j),lefty,jlo,mflag)
    ilo = kx
    do i=1,nx
      call interv (tx,nx+kx,taux(i),leftx, ilo, mflag )
      do jj=1,ky
        work(jj)=bvalue(tx,bcoef(:,lefty-ky+jj),nx,kx,taux(i),leftx)
      end do
      call interv(ty(lefty-ky+1:),ky+ky,tauy(j),left,klo,mflag)
      gtau(i,j) = bvalue(ty(lefty-ky+1:),work,ky,ky,tauy(j),left)
    end do
  end do
  call cpu_time(t2)

  do i=1,nx
    print 632,taux(i),(gtau(i,j)-g(taux(i),tauy(j)),j=1,ny)
  end do

  do j=1,ny
  do i=1,nx
    gtau(i,j)=gtau(i,j)-g(taux(i),tauy(j))
  end do
  end do

  print*, 'Max error                    ', maxval(abs(gtau))
  print*, 'Average error                ', sum(abs(gtau))/(nx*ny)
  print*, 'Time to compute interpolants ', t1-t0
  print*, 'Time to interpolate values   ', t2-t1
  print*, 'PASSED'

  620 format(' given data'//11f8.1)
  632 format(f5.1,11f8.4)
  640 format(//' interpolation error'//11f8.1)

contains

function g (x , y)

  sll_real64 :: x
  sll_real64 :: y
  sll_real64 :: g
  
  g = max(x-0.5_f64*nx,0.0_f64)**2 + max(y-0.5_f64*ny,0.0_f64)**3

end function g

!> Evaluates a derivative of a spline from its B-spline representation.
function bvalue( t, bcoef, n, k, x, i ) result(res)
    
sll_real64, intent(in)  :: t(:)
sll_real64, intent(in)  :: bcoef(:)
sll_int32,  intent(in)  :: n
sll_int32,  intent(in)  :: k
sll_real64, intent(in)  :: x
sll_int32,  intent(in)  :: i

sll_real64              :: res


sll_int32 :: ilo
sll_int32 :: j
sll_int32 :: jc
sll_int32 :: jcmax
sll_int32 :: jcmin
sll_int32 :: jj

res = bcoef(i)
jcmin = 1
if ( k <= i ) then
  do j = 1, k-1
    dl(j) = x - t(i+1-j)
  end do
else
  jcmin = 1-(i-k)
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
  jcmax = k+n-i
  do j = 1, k+n-i
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
do j = 1, k-1
  ilo = k-j
  do jj = 1, k-j
    aj(jj) = (aj(jj+1)*dl(ilo)+aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
    ilo = ilo - 1
  end do
end do

res = aj(1)

end function bvalue
  
end program test_bsplines_2d
