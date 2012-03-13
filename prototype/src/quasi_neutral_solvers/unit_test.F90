! gfortran qnefspl.f90 bsplvd.f90 bsplvb.f90 test_qnefspl.f90 -llapack -ldfftpack
program test_quasi_neutral
  use sll_qns_2d_with_finite_diff
  use sll_collective
#include "sll_remap.h"
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
  implicit none

end program test_quasi_neutral
