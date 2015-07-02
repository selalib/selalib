program test_bsplines_2d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_utilities.h"
use sll_bsplines
use sll_boundary_condition_descriptors

  implicit none

  sll_int32, parameter    :: nx=2800
  sll_int32, parameter    :: kx=3
  sll_int32, parameter    :: ny=2400
  sll_int32, parameter    :: ky=4

  sll_int32               :: i
  sll_int32               :: j
  sll_int32               :: jj
  sll_int32               :: left
  sll_int32               :: leftx
  sll_int32               :: lefty
  sll_int32               :: mflag
  sll_int32               :: ilo
  sll_int32               :: jlo
  sll_int32               :: klo
  sll_int32               :: llo
  sll_int32               :: ierr
  sll_int32               :: jjj
  sll_int32               :: kkk
  sll_int32               :: jc, jcmin, jcmax

  sll_real64              :: xi, xj
  sll_real64              :: gtau(nx,ny)
  sll_real64, pointer     :: bcoef(:,:)
  sll_real64, pointer     :: taux(:)
  sll_real64, pointer     :: tauy(:)
  sll_real64, pointer     :: tx(:)
  sll_real64, pointer     :: ty(:)
  sll_real64              :: work(nx)

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
  taux => bspline_2d%bs1%tau
  tx   => bspline_2d%bs1%t
  tauy => bspline_2d%bs2%tau
  ty   => bspline_2d%bs2%t
  
  ! generate and print out function values
  !print 620,(tauy(i),i=1,ny)
  do i=1,nx
    do j=1,ny
      gtau(i,j) = g(taux(i),tauy(j))
    end do
    !print 632,taux(i),(bcoef(i,j),j=1,ny)
  end do
  
  call cpu_time(t0)
  call compute_bspline_2d( bspline_2d, gtau) 
  bcoef => bspline_2d%bcoef
  call cpu_time(t1)
 
  ! evaluate interpolation error at mesh points and print out
  !print 640,(tauy(j),j=1,ny)

  SLL_CLEAR_ALLOCATE(aj(1:max(kx,ky)), ierr)
  SLL_CLEAR_ALLOCATE(dl(1:max(kx,ky)), ierr)
  SLL_CLEAR_ALLOCATE(dr(1:max(kx,ky)), ierr)

  jlo = ky
  do j=1,ny
    xj = tauy(j)
    call interv(ty,ny+ky,xj,lefty,jlo,mflag)
    ilo = kx
    klo = jlo
    do i=1,nx
      xi = taux(i)
      call interv (tx,nx+kx,xi,leftx, ilo, mflag )
      do jj=1,ky
        jcmin = 1
        if ( kx <= leftx ) then
          do jjj = 1, kx-1
            dl(jjj) = xi - tx(leftx+1-jjj)
          end do
        else
          jcmin = 1-(leftx-kx)
          do jjj = 1, leftx
            dl(jjj) = xi - tx(leftx+1-jjj)
          end do
          do jjj = leftx, kx-1
            aj(kx-jjj) = 0.0_f64
            dl(jjj) = dl(leftx)
          end do
        end if
        jcmax = kx
        if ( nx < leftx ) then
          jcmax = kx+nx-leftx
          do jjj = 1, kx+nx-leftx
            dr(jjj) = tx(leftx+jjj) - xi
          end do
          do jjj = kx+nx-leftx, kx-1
            aj(jjj+1) = 0.0_f64
            dr(jjj) = dr(kx+nx-leftx)
          end do
        else
          do jjj = 1, kx-1
            dr(jjj) = tx(leftx+jjj) - xi
          end do
        end if
        do jc = jcmin, jcmax
          aj(jc) = bcoef(leftx-kx+jc,lefty-ky+jj)
        end do
        do jjj = 1, kx-1
          llo = kx-jjj
          do kkk = 1, kx-jjj
            aj(kkk) = (aj(kkk+1)*dl(llo)+aj(kkk)*dr(kkk))/(dl(llo)+dr(kkk))
            llo = llo - 1
          end do
        end do
        work(jj) = aj(1)
      end do

      call interv(ty(lefty-ky+1:),ky+ky,xj,left,klo,mflag)
      jcmin = 1
      if ( ky <= left ) then
        do jjj = 1, ky-1
          dl(jjj) = xj - ty(lefty-ky+left+1-jjj)
        end do
      else
        jcmin = 1-(left-ky)
        do jjj = 1, left
          dl(jjj) = xj - ty(lefty-ky+left+1-jjj)
        end do
        do jjj = left, ky-1
          aj(ky-jjj) = 0.0_f64
          dl(jjj) = dl(left)
        end do
      end if
      jcmax = ky
      if ( ky < left ) then
        jcmax = ky+ky-left
        do jjj = 1, ky+ky-left
          dr(jjj) = ty(lefty-ky+left+jjj) - xj
        end do
        do jjj = ky+ky-left, ky-1
          aj(jjj+1) = 0.0_f64
          dr(jjj) = dr(ky+ky-left)
        end do
      else
        do jjj = 1, ky-1
          dr(jjj) = ty(lefty-ky+left+jjj) - xj
        end do
      end if
      do jc = jcmin, jcmax
        aj(jc) = work(left-ky+jc)
      end do
      do jjj = 1, ky-1
        llo = ky-jjj
        do kkk = 1, ky-jjj
          aj(kkk) = (aj(kkk+1)*dl(llo)+aj(kkk)*dr(kkk))/(dl(llo)+dr(kkk))
          llo = llo - 1
        end do
      end do
      gtau(i,j) = aj(1)
    end do
  end do
  call cpu_time(t2)

  !do i=1,nx
  !  print 632,taux(i),(gtau(i,j)-g(taux(i),xj),j=1,ny)
  !end do

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

end program test_bsplines_2d
