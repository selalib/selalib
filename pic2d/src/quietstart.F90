module quietstart
#include "sll_working_precision.h"
use zone

implicit none

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sll_real64 function bit_reversing( n )

sll_int32 :: n, deci, k, div
sll_real64 :: miroir

!a = 0
k = 25
miroir = 0.d0
deci = n

if (deci > 33554432) then 
   print*,'entier trop grand (superieur a 33554432=2**25)'
   stop
else 
   do while (k >= 0)
      div = 2**k
      if (deci/div == 1) then
        ! a = a + 10**k
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

sll_int32 :: deci, k, div, n
sll_real64 :: miroir

!a = 0
k = 16
miroir = 0.d0
deci = n

if (deci > 43046721) then 
   print*,'entier trop grand (superieur a 43046721=3**16)'
   stop
else 
   do while (k >= 0)
      div = 3**k
      if (deci/div == 1) then
        ! a = a + 10**k
         miroir = miroir + 3.**(-k-1)
         deci = deci - div
      else if (deci/div == 2) then
        ! a = a + 2 * 10**k
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

sll_int32 :: deci, k, div, n
sll_real64 :: miroir

!a = 0
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
        ! a = a + 10**k
         miroir = miroir + 5.**(-k-1)
         deci = deci - div
      else if (deci/div == 2) then
        ! a = a + 2 * 10**k
         miroir = miroir + 2 * 5.**(-k-1)
         deci = deci - 2 * div
      else if (deci/div == 3) then
        ! a = a + 3 * 10**k
         miroir = miroir + 3 * 5.**(-k-1)
         deci = deci - 3 * div
      else if (deci/div == 4) then
        ! a = a + 4 * 10**k
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

real(kind = prec) :: a, b, R, eps, x, D

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

end module quietstart
