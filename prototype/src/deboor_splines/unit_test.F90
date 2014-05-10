program unit_test
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_constants.h"
use sll_module_deboor_splines_1d

implicit none



sll_real64, dimension(:), pointer :: knots
sll_real64, dimension(:), pointer :: value_func
sll_real64, dimension(:), pointer :: value_der
sll_real64, dimension(:), pointer :: points,points_test
sll_int32, dimension(:), pointer :: points_der
sll_int32  :: dim_point, dim_point_der
sll_int32  :: deg_spline,ordre
sll_real64, dimension(:), pointer :: matrices
sll_real64, dimension(:), pointer :: bcoef
sll_int32 :: iflag,ierr
sll_real64 :: d,res
sll_int32 :: i,j



deg_spline    = 3
dim_point     = 10
dim_point_der = 2
ordre         = deg_spline +1

SLL_ALLOCATE(knots(dim_point + ordre + dim_point_der),ierr)
SLL_ALLOCATE(value_func(dim_point),ierr)
SLL_ALLOCATE(value_der(dim_point_der),ierr)
SLL_ALLOCATE(points(dim_point),ierr)
SLL_ALLOCATE(points_test(dim_point),ierr)
SLL_ALLOCATE(points_der(dim_point_der),ierr)
SLL_ALLOCATE(matrices((2*ordre-1)*(dim_point+dim_point_der)),ierr)
SLL_ALLOCATE(bcoef(dim_point+ dim_point_der),ierr)



do i = 1 , dim_point
   points(i) = 0.0_f64 + (i-1)*1./(dim_point-1)
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

call splint_der(&
     points,value_func,&
     points_der,value_der,&
     knots,dim_point,dim_point_der,ordre,&
     matrices, bcoef, iflag )



do i = 1 , dim_point
   res=bvalue( knots, bcoef, dim_point+dim_point_der, ordre, points(i), 0)
   !print*, 'approx',res,'exact', value_func(i)
end do
do i = 1,dim_point
   res=bvalue( knots, bcoef, dim_point+dim_point_der, ordre, points(i), 1)
   !print*, 'approx',res,'exact', value_der(i)
end do
   

SLL_DEALLOCATE(knots,ierr)
SLL_DEALLOCATE(value_func,ierr)
SLL_DEALLOCATE(value_der,ierr)
SLL_DEALLOCATE(points,ierr)
SLL_DEALLOCATE(points_der,ierr)
SLL_DEALLOCATE(matrices,ierr)
SLL_DEALLOCATE(bcoef,ierr)

print*, 'PASSED'
end program unit_test
