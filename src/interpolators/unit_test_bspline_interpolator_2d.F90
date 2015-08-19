program bspline_2d_dirichlet

#include "sll_working_precision.h"
#include "sll_constants.h"

use sll_module_bspline_interpolator_2d
use sll_boundary_condition_descriptors

implicit none

#define NPTS1    65
#define NPTS2    65
#define SPL_DEG1 3
#define SPL_DEG2 3
#define X1MIN    0.0_f64
#define X1MAX    1.0_f64
#define X2MIN    0.0_f64
#define X2MAX    1.0_f64

type(sll_bspline_interpolator_2d), pointer :: interpolator

sll_real64, dimension(NPTS1,NPTS2) :: x
sll_real64, dimension(NPTS1,NPTS2) :: y
sll_real64, dimension(NPTS1,NPTS2) :: g_int, dg_dx_int, dg_dy_int
sll_real64, dimension(NPTS1,NPTS2) :: g_ref, dg_dx_ref, dg_dy_ref

sll_int32  :: i 
sll_int32  :: j 
sll_real64 :: h1
sll_real64 :: h2
sll_real64 :: normL2
sll_real64 :: normH1
  
h1 = (X1MAX-X1MIN)/real(NPTS1-1,f64)
h2 = (X2MAX-X2MIN)/real(NPTS2-1,f64)
  
print *, '***********************************************************'
print *, '              Dirichlet'
print *, '***********************************************************'
  
do j=1,NPTS2
do i=1,NPTS1
  y(i,j)  = f(X1MIN+(i-1)*h1,X2MIN+(j-1)*h2)
end do
end do
  
!call interpolator%initialize(NPTS1,X1MIN,X2MAX,SPL_DEG1,SLL_DIRICHLET)
!
!call interpolator%compute_interpolants(y)
!  
normL2 = 0.0_f64
normH1 = 0.0_f64
do j=1,NPTS2
do i=1,NPTS1
  g_int(i,j)  = interpolator%interpolate_value(x(i,j),y(i,j))
  g_ref(i,j)  = f(x(i,j),y(i,j))
  dg_dx_int(i,j) = interpolator%interpolate_derivative_eta1(x(i,j),y(i,j))
  dg_dx_ref(i,j) = df_dx(x(i,j),y(i,j))
  write(10,*) x(i,j), y(i,j), g_int(i,j), g_ref(i,j)
  write(11,*) x(i,j), y(i,j), dg_dx_int(i,j), dg_dx_ref(i,j)
  write(12,*) x(i,j), y(i,j), dg_dy_int(i,j), dg_dy_ref(i,j)
end do
end do
!
!normL2 = sum((y_int-y_ref)**2*h)
!normH1 = sum((dy_int-dy_ref)**2*h)
!  
!print*,'--------------------------------------------'
!print*,' Average error in nodes', sum(abs(y_int-y_ref))/NPTS
!print*,' Max     error in nodes', maxval(abs(y_int-y_ref))
!print*,'--------------------------------------------'
!print*,' Average error in nodes first derivative',sum(abs(dy_int-dy_ref))/NPTS
!print*,' Max     error in nodes first derivative',maxval(abs(dy_int-dy_ref))
!print*,'--------------------------------------------'
!print*,' Norm L2 error ', sqrt(normL2), h**(SPL_DEG)
!print*,'--------------------------------------------'
!print*,' Norm H1 error ', sqrt(normH1), h**(SPL_DEG-2)
!print*,'--------------------------------------------'
!
!if(( sqrt(normL2) <= h**(SPL_DEG)) .AND. &
!   ( sqrt(normH1) <= h**(SPL_DEG-2))) then
!  print *, 'PASSED'
!else
!  print *, 'FAILED'
!end if
!
!call sll_delete(interpolator)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

  function f(x,y)
  
    sll_real64 :: x
    sll_real64 :: y
    sll_real64 :: f
  
    f = sin(2.0_f64*sll_pi*x)*sin(2.0_f64*sll_pi*y)
  
  end function f
  
  function df_dx(x,y)
  
    sll_real64 :: x
    sll_real64 :: y
    sll_real64 :: df_dx
  
    df_dx = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*x)*sin(2.0_f64*sll_pi*y)
  
  end function df_dx

  function df_dy(x,y)
  
    sll_real64 :: x
    sll_real64 :: y
    sll_real64 :: df_dy
  
    df_dy = 2.0_f64*sll_pi*cos(2.0_f64*sll_pi*y)*sin(2.0_f64*sll_pi*x)
  
  end function df_dy

end program bspline_2d_dirichlet
