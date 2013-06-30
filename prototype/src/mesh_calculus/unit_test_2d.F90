program test_mesh_calculus
#include "sll_working_precision.h" 
#include "sll_assert.h"
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_logical_meshes
  use sll_mesh_calculus_2d_module
  implicit none

#define NCELLS1 1
#define NCELLS2 1

  type(sll_logical_mesh_2d), pointer                    :: m
  class(sll_coordinate_transformation_2d_base), pointer :: Ta
  sll_real64 :: volume, length_east

  m => new_logical_mesh_2d( NCELLS1, NCELLS2)

  Ta => new_coordinate_transformation_2d_analytic( &
       "map_a", &
       m, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22 )

  volume = cell_volume(Ta,1,1,3)
  print *, 'volume = ', volume

  length_east = edge_length_eta1_plus(Ta,1,1,3)
  print *, 'length east edge = ', length_east


  print *, 'PASSED'

end program test_mesh_calculus
