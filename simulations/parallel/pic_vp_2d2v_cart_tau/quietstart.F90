module m_quietstart
#include "sll_working_precision.h"
use m_zone

implicit none

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sll_real64 function bit_reversing( n )

sll_int32  :: n
sll_int32  :: deci
sll_int32  :: k
sll_int32  :: div
sll_real64 :: miroir

k      = 25
miroir = 0.d0
deci   = n

if (deci > 33554432) then 
  stop 'entier trop grand (superieur a 33554432=2**25)'
else 
  do while (k >= 0)
    div = 2**k
    if (deci/div == 1) then
      miroir = miroir + 2.**(-k-1)
      deci = deci - div
    endif
    k = k-1
 enddo
endif

bit_reversing = miroir

end function bit_reversing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sll_real64 function trinary_reversing( n )

sll_int32  :: deci
sll_int32  :: k
sll_int32  :: div
sll_int32  :: n
sll_real64 :: miroir

k = 16
miroir = 0.d0
deci = n

if (deci > 43046721) then 
  stop 'entier trop grand (superieur a 43046721=3**16)'
else 
  do while (k >= 0)
    div = 3**k
    if (deci/div == 1) then
      miroir = miroir + 3.**(-k-1)
      deci = deci - div
    else if (deci/div == 2) then
      miroir = miroir + 2 * 3.**(-k-1)
      deci = deci - 2 * div
    endif
    k = k-1
  enddo
endif

trinary_reversing = miroir

end function trinary_reversing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sll_real64 function penta_reversing( n )

sll_int32  :: deci
sll_int32  :: k
sll_int32  :: div
sll_int32  ::  n
sll_real64 :: miroir

k = 11
miroir = 0.d0
deci = n

if (deci > 48828125) then 
  print*,'entier trop grand (superieur a 48828125=5**11)'
  stop
else 
  do while (k >= 0)
    div = 5**k
    if (deci/div == 1) then
      miroir = miroir + 5.**(-k-1)
      deci = deci - div
    else if (deci/div == 2) then
      miroir = miroir + 2 * 5.**(-k-1)
      deci = deci - 2 * div
    else if (deci/div == 3) then
      miroir = miroir + 3 * 5.**(-k-1)
      deci = deci - 3 * div
    else if (deci/div == 4) then
      miroir = miroir + 4 * 5.**(-k-1)
      deci = deci - 4 * div
    endif
    k = k-1
  enddo
endif

penta_reversing = miroir

end function penta_reversing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dichotomie_x( a, b, R, eps ) 

! il faut D(a)<R<D(b), on cherche x tq R=D(x), resu dans a 

sll_real64 :: a
sll_real64 :: b
sll_real64 :: R
sll_real64 :: eps
sll_real64 :: x
sll_real64 :: D

D = ( kx*a + alpha * sin(kx*a) ) / (2*pi)
do while ( D<R-eps .or. D>R+eps )
  x = (a+b)/2
  D = ( kx*x + alpha * sin(kx*x) ) / (2*pi)
  if ( D<R-eps ) then
    a = x
  else if ( D>R+eps ) then 
    b = x
  else
    a = x
  endif
end do

end subroutine dichotomie_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module m_quietstart
