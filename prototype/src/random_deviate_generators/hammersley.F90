!> hamm pseudo-random generator
!>
!>\author 
!>\date created: 
module hammersley
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

end module hammersley
