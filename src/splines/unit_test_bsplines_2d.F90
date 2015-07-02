program test_bsplines_2d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_utilities.h"
use sll_bsplines
use sll_boundary_condition_descriptors

  implicit none

  sll_int32, parameter    :: nx=14
  sll_int32, parameter    :: kx=3
  sll_int32, parameter    :: ny=11
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
  sll_real64              :: htau(nx,ny)
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
  tauy => bspline_2d%bs2%tau
  tx   => bspline_2d%bs1%t
  ty   => bspline_2d%bs2%t
  
  ! generate and print out function values
  print 620,(tauy(i),i=1,ny)
  do i=1,nx
    do j=1,ny
      gtau(i,j) = g(taux(i),tauy(j))
    end do
    print 632,taux(i),(gtau(i,j),j=1,ny)
  end do
  
  call cpu_time(t0)
  call compute_bspline_2d( bspline_2d, gtau) 
  bcoef => bspline_2d%bcoef
  call cpu_time(t1)
 
  ! evaluate interpolation error at mesh points and print out
  call interpolate_array_values_2d(bspline_2d, nx, ny, gtau, htau)

  call cpu_time(t2)

  print 640,(tauy(j),j=1,ny)
  do i=1,nx
    print 632,taux(i),(htau(i,j)-g(taux(i),tauy(j)),j=1,ny)
  end do

  do j=1,ny
    do i=1,nx
      gtau(i,j)=htau(i,j)-g(taux(i),tauy(j))
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
