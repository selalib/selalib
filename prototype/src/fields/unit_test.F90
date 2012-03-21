program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
  use numeric_constants
  use geometry_functions
  !use sll_mapped_mesh_base
  use sll_mapped_meshes
  implicit none
  
  type(sll_mapped_mesh_2d_analytic), target :: mesh
  class(sLL_mapped_mesh_2d_base), pointer :: m
  sll_int32 :: nc1, nc2
  procedure(polar_x1), pointer :: px1, px2, pjac11, pjac12, pjac21, pjac22

  nc1 = 10
  nc2 = 10
  px1 => sinprod_x1
  px2 => sinprod_x2
  pjac11 => sinprod_jac11
  pjac12 => sinprod_jac12
  pjac21 => sinprod_jac21
  pjac22 => sinprod_jac22
  call mesh%initialize( nc1, nc2, px1, px2, pjac11, pjac12, pjac21, pjac22)
  m => mesh

  print*, m%x1_at_node(5,3), m%x1(.3_f64, .4_f64)

  print *, 'Successful, exiting program.'
  

  
end program unit_test
