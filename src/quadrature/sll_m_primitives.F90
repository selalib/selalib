!> @ingroup utilities
!> Functions to compute primitive of 1d function.
module sll_m_primitives
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_s_function_to_primitive, &
    sll_s_primitive_to_function

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

!> Compute primitive of f 
!> @param[in]  f  1d array of double float
!> @param[in]  x  nodes positions
!> @param[in]  n  nodes number
!> @param[out] xm mean function
subroutine sll_s_function_to_primitive(f,x,n,xm)

sll_real64,dimension(:),intent(inout) :: f
sll_real64,dimension(:),intent(in)    :: x
sll_int32,              intent(in)    :: n
sll_real64,             intent(out)   :: xm

sll_int32  :: i
sll_real64 :: dx
sll_real64 :: tmp
sll_real64 :: tmp2

dx = 1._f64/real(n,f64)
    
!from f compute the mean
xm=0.0_f64
do i=1,n
  xm=xm+f(i)*(x(i+1)-x(i))
enddo

f(1) = (f(1)-xm)*(x(2)-x(1))
tmp  = f(1)
f(1) = 0.0_f64
do i=2,n
  f(i) = (f(i)-xm)*(x(i+1)-x(i))
  tmp2 = f(i)
  f(i) = f(i-1) + tmp
  tmp  = tmp2
enddo    
f(n+1) = f(n)+tmp

end subroutine sll_s_function_to_primitive

!> Compute function value from primitive 
!> @param[inout]  f  1d array of double float
!> @param[in]     x  nodes positions
!> @param[in]     n  nodes number
!> @param[in]     xm mean function
subroutine sll_s_primitive_to_function(f,x,n,xm)

sll_real64, dimension(:), intent(inout) :: f
sll_real64, dimension(:), intent(in)    :: x
sll_int32,                intent(in)    :: n
sll_real64,               intent(in)    :: xm

sll_int32                               :: i
sll_real64                              :: tmp
sll_real64                              :: dx 

dx = 1._f64/real(n,f64)

tmp=f(1)
do i=1,n-1
  f(i)=f(i+1)-f(i)+xm*(x(i+1)-x(i))
enddo
f(n)=tmp-f(n)+xm*(x(1)+1._f64-x(n))

!from mean compute f
do i=1,n
  f(i)=f(i)/(x(i+1)-x(i))
enddo

f(n+1) = f(1)

end subroutine sll_s_primitive_to_function

end module sll_m_primitives
