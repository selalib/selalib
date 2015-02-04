module sll_primitives
#include "sll_working_precision.h"

contains

subroutine function_to_primitive(f,node_positions,N,M)

sll_real64,dimension(:),intent(inout) :: f
sll_real64,dimension(:),intent(in)    :: node_positions
sll_int32,              intent(in)    :: N
sll_real64,             intent(out)   :: M

sll_int32  :: i
sll_real64 :: dx,tmp,tmp2

dx = 1._f64/real(N,f64)
    
!from f compute the mean
M=0._f64
do i=1,N
  M=M+f(i)*(node_positions(i+1)-node_positions(i))
enddo

f(1)=(f(1)-M)*(node_positions(2)-node_positions(1))
tmp=f(1)
f(1)=0._f64
do i=2,N!+1
  f(i)=(f(i)-M)*(node_positions(i+1)-node_positions(i))
  tmp2=f(i)
  f(i)=f(i-1)+tmp
  tmp=tmp2
enddo    
f(N+1)=f(N)+tmp

end subroutine function_to_primitive

subroutine primitive_to_function(f,node_positions,N,M)

sll_real64, dimension(:), intent(inout) :: f
sll_real64, dimension(:), intent(in)    :: node_positions
sll_int32,                intent(in)    :: N
sll_real64,               intent(in)    :: M

sll_int32                               :: i
sll_real64                              :: tmp
sll_real64                              :: dx 

dx = 1._f64/real(N,f64)

tmp=f(1)
do i=1,N-1
  f(i)=f(i+1)-f(i)+M*(node_positions(i+1)-node_positions(i))
enddo
f(N)=tmp-f(N)+M*(node_positions(1)+1._f64-node_positions(N))

!from mean compute f
do i=1,N
  f(i)=f(i)/(node_positions(i+1)-node_positions(i))
enddo

f(N+1) = f(1)

end subroutine primitive_to_function

end module sll_primitives
