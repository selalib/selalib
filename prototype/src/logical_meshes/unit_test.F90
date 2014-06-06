program unit_test_logical_meshes
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_logical_meshes
  implicit none

  type(sll_logical_mesh_2d), pointer :: m2d
  type(sll_logical_mesh_3d), pointer :: m3d
  type(sll_logical_mesh_4d), pointer :: m4d

  type(sll_logical_mesh_1d), pointer :: m1d_x
  type(sll_logical_mesh_1d), pointer :: m1d_y
  type(sll_logical_mesh_2d), pointer :: m2d_xy

  type(sll_logical_mesh_2d), pointer :: mx
  type(sll_logical_mesh_2d), pointer :: mv
  type(sll_logical_mesh_4d), pointer :: mxv

  m1d_x => new_logical_mesh_1d(64)
  m1d_y => new_logical_mesh_1d(32,eta_min=-1.0_f64,eta_max=+1.0_f64)

  call sll_display(m1d_x)
  call sll_display(m1d_y)

#ifndef STDF95
  m2d_xy => m1d_x * m1d_y
  call sll_display(m2d_xy)
  call delete(m2d_xy)
#endif

  call delete(m1d_x)
  call delete(m1d_y)

  m2d => new_logical_mesh_2d(100,100)
  m3d => new_logical_mesh_3d(100,100,100)
  m4d => new_logical_mesh_4d(32,32,32,32, eta1_min=-1.0_f64, eta1_max = 2.0_f64)

  call sll_display(m2d)
  call sll_display(m3d)
  call sll_display(m4d)

  call delete(m2d)
  call delete(m3d)
  call delete(m4d)

  mx => new_logical_mesh_2d(100,100, 0.0_f64, 12.56_f64, 0.0_f64, 12.56_f64)
  mv => new_logical_mesh_2d(64,64,-6.0_f64,6.0_f64,-6.0_f64,6.0_f64)

#ifndef STDF95
  mxv => mx * mv
  call sll_display(mxv)
  call delete(mxv)
#endif

  call delete(mx)
  call delete(mv)

  call unit_test_logical_meshes_test_1d_functionality()

  print *, 'PASSED'



contains

subroutine unit_test_logical_meshes_test_1d_functionality()
  type(sll_logical_mesh_1d), pointer :: mesh
  sll_real64, dimension(9) :: nodes
  sll_real64, dimension(9) :: celln
  mesh=>new_logical_mesh_1d(8, 0.0_f64, 1.0_f64)
  nodes=sll_mesh_nodes(mesh)
  celln=sll_cell(mesh,nodes)
print *, nodes
!  call sll_display(nodes, "(F8.4)")
!  call sll_display(celln, "(F8.4)")



endsubroutine

end program unit_test_logical_meshes


