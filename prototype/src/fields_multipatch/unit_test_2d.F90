program unit_test_fields_multipatch
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_module_scalar_field_2d_multipatch
  use sll_coordinate_transformation_multipatch_module
  implicit none

  
  type(sll_coordinate_transformation_multipatch_2d), pointer :: T
  type(sll_scalar_field_multipatch_2d), pointer              :: F

  T => new_coordinate_transformation_multipatch_2d("circle_5mp_info.nml")
  print *, 'initialized multipatch transformation'
  
  F => new_scalar_field_multipatch_2d("test_field_multipatch", T)
  

     
  print *, 'PASSED'
  
end program unit_test_fields_multipatch

