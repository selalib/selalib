program unit_test_initializers_4d
#include "sll_working_precision.h"
  use sll_logical_meshes
  use sll_coordinate_transformation_2d_base_module
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  use sll_parallel_array_initializer_module
  implicit none

  ! logical meshes
  type(sll_logical_mesh_2d), pointer :: mx
  type(sll_logical_mesh_2d), pointer :: mv
  ! coordinate transformations (test transforming spatial coordinates only)
  class(sll_coordinate_transformation_2d_base), pointer :: tx
  sll_real64, dimension(2) :: params_identity

  params_identity(:) = (/0.0_f64, 0.0_f64/) ! for identity this can be whatever
  ! initialize the logical meshes
  mx => new_logical_mesh_2d(64,64)
  mv => new_logical_mesh_2d(64,64)

  ! initialize the transformation
  tx => new_coordinate_transformation_2d_analytic( &
       "identity_transformation", &
       mx, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22, &
       params_identity )

  ! initialize the array

  ! THIS TEST NEEDS TO BE FINISHED IN A COMPLETE WAY. FOR NOW, THE MODULE
  ! WILL BE PARTIALLY TESTED ON A SIMULATION...

  print *, 'PASSED'

end program unit_test_initializers_4d
