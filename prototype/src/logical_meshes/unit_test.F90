program unit_test_logical_meshes
  use sll_logical_meshes
  implicit none

  type(sll_logical_mesh_2d), pointer :: m2d
  type(sll_logical_mesh_4d), pointer :: m4d

  m2d => new_logical_mesh_2d(100,100)
  m4d => new_logical_mesh_4d(32,32,32,32, eta1_min=-1.0_f64, eta1_max = 2.0_f64)

  call delete(m2d)
  call delete(m4d)

  print *, 'PASSED'

end program unit_test_logical_meshes
