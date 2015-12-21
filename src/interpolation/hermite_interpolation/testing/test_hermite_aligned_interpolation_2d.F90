program test_hermite_aligned_interpolation_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"

  use sll_m_hermite_aligned_interpolation_2d, only: &
    sll_f_new_hermite_aligned_interpolation_2d, &
    sll_t_hermite_aligned_interpolation_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  type(sll_t_hermite_aligned_interpolation_2d), pointer :: interp

  interp => sll_f_new_hermite_aligned_interpolation_2d()
  
  print *,'#PASSED'



end program test_hermite_aligned_interpolation_2d
