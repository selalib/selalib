program accumulators_unit_test
  use sll_electric_field_2d_accumulator
#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none

  type(electric_field_2d_accumulator), dimension(:,:), allocatable :: efield
  sll_int32 :: ierr
  SLL_ALLOCATE(efield(16,16),ierr)

  print *, 'PASSED'
end program accumulators_unit_test
