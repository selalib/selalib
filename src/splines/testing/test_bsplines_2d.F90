!PN This test is only for bsplines 2d
!PN It is used to test different boundary conditions
!PN Data values are also displayed on the screen
!PN This test is dedicated for debug not for continuous integration
program test_bsplines_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_periodic

  use sll_m_bsplines, only: &
    compute_bspline_2d, &
    interpolate_array_values_2d, &
    new_bspline_2d, &
    sll_bspline_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_int32, parameter    :: nx=14
sll_int32, parameter    :: kx=3
sll_int32, parameter    :: ny=11
sll_int32, parameter    :: ky=4

sll_int32               :: i
sll_int32               :: j

sll_real64              :: gtau(nx,ny)
sll_real64              :: htau(nx,ny)
sll_real64, pointer     :: taux(:)
sll_real64, pointer     :: tauy(:)


type(sll_bspline_2d), pointer :: bspline_2d

sll_real64 :: t0, t1, t2

bspline_2d => new_bspline_2d( nx, kx-1, 1.0_f64, nx*1.0_f64, SLL_PERIODIC, &
                              ny, ky-1, 1.0_f64, ny*1.0_f64, SLL_PERIODIC  )

! set up data points and knots in x, interpolate between knots by parabolic 
! splines, using not-a-knot end condition

taux => bspline_2d%bs1%tau
tauy => bspline_2d%bs2%tau

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
call cpu_time(t1)

! evaluate interpolation error at mesh points and print out
call interpolate_array_values_2d(bspline_2d, nx, ny, gtau, htau, 0, 0)
!do j=1,ny
!  do i=1,nx
!    htau(i,j) = interpolate_value_2d(bspline_2d, taux(i), tauy(j), 0, 0)
!  end do
!end do

print 630,(tauy(j),j=1,ny)
do i=1,nx
  print 632,taux(i),(htau(i,j),j=1,ny)
end do

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
print*, 'Average error                ', sum(abs(gtau))/real(nx*ny,f64)
print*, 'Time to compute interpolants ', t1-t0
print*, 'Time to interpolate values   ', t2-t1
print*, 'PASSED'

620 format(' given data'//5x,11f8.1)
630 format(' interpolated data'//5x,11f8.1)
632 format(f5.1,11f8.4)
640 format(//' interpolation error'//5x,11f8.1)

contains

  function g(x,y)
  
    sll_real64 :: x
    sll_real64 :: y
    sll_real64 :: g
    
    g = max(x-0.5_f64*nx,0.0_f64)**2 + max(y-0.5_f64*ny,0.0_f64)*3.0
  
  end function g

end program test_bsplines_2d
