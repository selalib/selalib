program unit_test
#include "sll_working_precision.h"
#include "sll_memory.h"
  use numeric_constants
  use geometry_functions
  use sll_mesh_2d
  implicit none
  
  type(mesh_2d_analytic), target :: mesh
  class(mesh_2d), pointer :: m
  sll_int32 :: nc1, nc2
  procedure(polar_x1), pointer :: px1, px2, pjac

  nc1 = 10
  nc2 = 10
  px1 => sinprod_x1
  px2 => sinprod_x2
  pjac => sinprod_jac
  call new_mesh_2d_analytic ( mesh, nc1, nc2, px1, px2, pjac)
  m => mesh

  print*, m%x1_at_node(5,3), m%x1(.3_f64, .4_f64), m%x1_array(5,3)

  print *, 'Successful, exiting program.'
  

  
end program unit_test
