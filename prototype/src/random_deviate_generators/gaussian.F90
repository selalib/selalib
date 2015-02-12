!> gaussian random generator
!>
!>\author 
!>\date created: 
module gaussian
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none

contains 

  !>Returns a value of random variable following the Gaussian probability 
  !>density with zero mean and unit variance, using standart function "random" 
  !>as the source of uniform deviates.
  !>It's a copy-and-paste from Numerical Recipes
  !>\author Sever Hirstoaga
  !>\date created: 2012-05-31
  function gaussian_deviate()
    sll_real64 :: gaussian_deviate

    sll_real64 :: rsq, v1, v2
    sll_real64, save :: g
    logical, save :: g_stored=.false.
    
    if (g_stored) then
       gaussian_deviate=g
       g_stored=.false.
    else
       do
          call random_number(v1)
          call random_number(v2)
          v1=2._f64*v1 - 1._f64
          v2=2._f64*v2 - 1._f64
          rsq = v1**2+v2**2
          if (rsq>0._f64.and.rsq<1._f64) exit
       enddo
       rsq = sqrt(-2._f64*log(rsq)/rsq)
       gaussian_deviate=v1*rsq
       g=v2*rsq
       g_stored=.true.
    endif
  end function gaussian_deviate

end module gaussian
