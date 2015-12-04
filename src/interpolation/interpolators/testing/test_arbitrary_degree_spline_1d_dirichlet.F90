program arbitrary_degree_splines_1d_dirichlet
#include "sll_working_precision.h"

use sll_m_arbitrary_degree_spline_interpolator_1d
use sll_m_boundary_condition_descriptors
use sll_m_constants, only : &
     sll_pi

implicit none

#define NPTS 65
#define SPL_DEG 3
#define XMIN 0.0_f64
#define XMAX 1.0_f64

type(sll_arbitrary_degree_spline_interpolator_1d) :: interpolator

sll_real64, dimension(NPTS) :: x
sll_real64, dimension(NPTS) :: y
sll_real64, dimension(NPTS) :: y_int, dy_int
sll_real64, dimension(NPTS) :: y_ref, dy_ref

sll_int32  :: i 
sll_real64 :: h
sll_real64 :: normL2
sll_real64 :: normH1
  
h = (XMAX-XMIN)/real(NPTS-1,f64)
  
print *, '***********************************************************'
print *, '              Dirichlet'
print *, '***********************************************************'
  
do i=1,NPTS
  y(i)  = f(XMIN + (i-1)*h)
end do
call random_number(x)
x = x * (XMAX-XMIN)
  
call interpolator%initialize(NPTS,XMIN,XMAX,SLL_DIRICHLET,SLL_DIRICHLET,SPL_DEG)
call set_values_at_boundary1d(interpolator,value_left=1.0_f64,value_right=1.0_f64)

call interpolator%compute_interpolants(y)
  
normL2 = 0.0_f64
normH1 = 0.0_f64
do i=1,NPTS
  y_int(i)  = interpolator%interpolate_from_interpolant_value(x(i))
  y_ref(i)  = f(x(i))
  dy_int(i) = interpolator%interpolate_from_interpolant_derivative_eta1(x(i))
  dy_ref(i) = df(x(i))
  write(10,*) x(i), y_int(i), y_ref(i)
  write(11,*) x(i), dy_int(i), dy_ref(i)
end do

normL2 = sum((y_int-y_ref)**2*h)
normH1 = sum((dy_int-dy_ref)**2*h)
  
print*,'--------------------------------------------'
print*,' Average error in nodes', sum(abs(y_int-y_ref))/NPTS
print*,' Max     error in nodes', maxval(abs(y_int-y_ref))
print*,'--------------------------------------------'
print*,' Average error in nodes first derivative',sum(abs(dy_int-dy_ref))/NPTS
print*,' Max     error in nodes first derivative',maxval(abs(dy_int-dy_ref))
print*,'--------------------------------------------'
print*,' Norm L2 error ', sqrt(normL2), h**(SPL_DEG)
print*,'--------------------------------------------'
print*,' Norm H1 error ', sqrt(normH1), h**(SPL_DEG-2)
print*,'--------------------------------------------'

if(( sqrt(normL2) <= h**(SPL_DEG)) .AND. &
   ( sqrt(normH1) <= h**(SPL_DEG-2))) then
  print *, 'PASSED'
else
  print *, 'FAILED'
end if

call sll_delete(interpolator)

contains

function f(x)

  sll_real64 :: x
  sll_real64 :: f

  f = sin(2.0_f64*sll_pi*x)+1

end function f

function df(x)

  sll_real64 :: x
  sll_real64 :: df

  df = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*x)

end function df


end program arbitrary_degree_splines_1d_dirichlet
