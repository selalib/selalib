program unit_test_2d
#include "sll_working_precision.h"
  use sll_coordinate_transformation_multipatch_module
  implicit none

  type(sll_coordinate_transformation_multipatch_2d) :: mp

  call read_from_file_multipatch_2d(mp, "circle_5mp_info.nml")

  print *, 'connectivity patch 1, face 3', get_connectivity(mp, 1, 3)
  print *, 'connectivity patch 3, face 1', get_connectivity(mp, 3, 1)
  print *, "PASSED"

end program unit_test_2d
