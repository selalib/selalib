program test_hermite_aligned_interpolation_2d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use sll_m_hermite_aligned_interpolation_2d
  use sll_m_constants
  implicit none
  
  type(sll_hermite_aligned_interpolation_2d), pointer :: interp

  interp => new_hermite_aligned_interpolation_2d()
  
  print *,'#PASSED'



end program test_hermite_aligned_interpolation_2d
