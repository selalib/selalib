program unit_test_fields_multipatch
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_module_scalar_field_2d_multpatch
  use sll_coordinate_transformation_multipatch_module
  implicit none

  
  class(sll_coordinate_transformation_multipatch_2d), pointer :: T
 
  T => new_coordinate_transformation_multipatch_2d("circle_5mp_info.nml")

  
  print *, 'initialized mesh 2D'
  
  ! coordinate transformation
  T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22, &
       params_identity )
  print *, 'initialized transformation'
  
     
  print *, 'PASSED'
  
end program unit_test_fields_multipatch

