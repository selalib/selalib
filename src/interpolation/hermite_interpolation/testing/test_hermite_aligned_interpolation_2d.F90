program test_hermite_aligned_interpolation_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"

  use sll_m_hermite_aligned_interpolation_2d, only: &
    new_hermite_aligned_interpolation_2d, &
    sll_hermite_aligned_interpolation_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  type(sll_hermite_aligned_interpolation_2d), pointer :: interp

  interp => new_hermite_aligned_interpolation_2d()
  
  print *,'#PASSED'



end program test_hermite_aligned_interpolation_2d
