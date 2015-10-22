!> hamm pseudo-random generator
!>
!>\author 
!>\date created: 
module sll_m_hammersley
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none

contains 

!! !> @brief Returns
!! !> @param[in]
  function suite_hamm (n,b)
    sll_real64 :: suite_hamm
    sll_int32  :: n,m
    sll_int32  :: b
    sll_int32  :: u
    sll_real64 :: h,k
    
    k=0
    h=0
    m=n
    do while (m>0)
       k = k+1
       u = m / b
       h = h +  b**(-k) * (m - u*b)
       m = u
    enddo
    suite_hamm = h 
  end function suite_hamm
  
  function vandercorput(n, p1, p2)!  p1>p2  !
    sll_int32, intent(in) :: n
    sll_int32, intent(in) :: p1, p2
    sll_real64 :: vandercorput, s
    sll_int32 :: m

    m=n
    vandercorput = 0.d0
    s=1
    do while (m>0)
       s=s/p1
       vandercorput=vandercorput+s*mod(p2*mod(m,p1),p1)
       m=int(m/p1)
    enddo
  end function vandercorput

end module sll_m_hammersley
