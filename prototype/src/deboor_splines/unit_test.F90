program test_deboor_splines
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_constants.h"
use sll_module_deboor_splines_1d

implicit none

sll_int32  :: dim_point
sll_int32  :: dim_point_der
sll_int32  :: deg_spline
sll_int32  :: ordre
sll_int32  :: iflag
sll_int32  :: ierr
sll_real64 :: d
sll_real64 :: res
sll_int32  :: i
sll_int32  :: j
sll_real64 :: tstart
sll_real64 :: tend

sll_real64, dimension(:), allocatable :: knots
sll_real64, dimension(:), allocatable :: value_func
sll_real64, dimension(:), allocatable :: value_der
sll_real64, dimension(:), allocatable :: points,points_test
sll_int32,  dimension(:), allocatable :: points_der
sll_real64, dimension(:), allocatable :: matrices
sll_real64, dimension(:), allocatable :: bcoef


deg_spline    = 3
dim_point     = 1024
dim_point_der = 2
ordre         = deg_spline+1

SLL_ALLOCATE(knots(dim_point+ordre+dim_point_der),ierr)
SLL_ALLOCATE(value_func(dim_point),ierr)
SLL_ALLOCATE(value_der(dim_point_der),ierr)
SLL_ALLOCATE(points(dim_point),ierr)
SLL_ALLOCATE(points_test(dim_point),ierr)
SLL_ALLOCATE(points_der(dim_point_der),ierr)
SLL_ALLOCATE(matrices((2*ordre-1)*(dim_point+dim_point_der)),ierr)
SLL_ALLOCATE(bcoef(dim_point+ dim_point_der),ierr)

call cpu_time(tstart)

do j = 1, 10000

  do i = 1 , dim_point
    points(i)     = (i-1.)/(dim_point-1)
    value_func(i) = cos(2*sll_pi*points(i))
  end do

  d = sum( abs( value_func(2:dim_point)-value_func(1:dim_point-1)))
  points_test(1) =  points(1)
  points_test(dim_point) =  points(dim_point)
  
  do i = 2,dim_point
    points_test(i) = points_test(i-1) +  abs( value_func(i)-value_func(i-1))/d
  end do

  knots(1:ordre) = points(1)

  do i = ordre+1, dim_point + deg_spline-1
    knots(i) = points(i-deg_spline)
  end do
  
  do i = dim_point +deg_spline,dim_point + ordre + dim_point_der
    knots(i) = points(dim_point)
  end do

  points_der(1) = 1
  value_der(1)  = -sin(2*sll_pi*points(1))*2*sll_pi
  points_der(2) = dim_point
  value_der(2)  = -sin(2*sll_pi*points(dim_point))*2*sll_pi

  call splint_der( points,          &
                   value_func,      &
                   points_der,      &
                   value_der,       &
                   knots,           &
                   dim_point,       &
                   dim_point_der,   &
                   ordre,           &
                   matrices,        &
                   bcoef,           &
                   iflag )

  do i = 1 , dim_point
    res=bvalue( knots, bcoef, dim_point+dim_point_der, ordre, points(i), 0)
  end do

  do i = 1,dim_point
    res=bvalue( knots, bcoef, dim_point+dim_point_der, ordre, points(i), 1)
  end do

end do

call cpu_time(tend)

print*, tend - tstart

DEALLOCATE(knots)
DEALLOCATE(value_func)
DEALLOCATE(value_der)
DEALLOCATE(points)
DEALLOCATE(points_der)
DEALLOCATE(matrices)
DEALLOCATE(bcoef)

print*, 'PASSED'
end program test_deboor_splines
