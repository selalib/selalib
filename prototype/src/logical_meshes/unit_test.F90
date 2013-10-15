program unit_test_logical_meshes
  use sll_logical_meshes
  implicit none


  type(sll_logical_mesh_2d), pointer :: m2d
  type(sll_logical_mesh_3d), pointer :: m3d
  type(sll_logical_mesh_4d), pointer :: m4d

  type(sll_logical_mesh_1d), pointer :: m1d_x
  type(sll_logical_mesh_1d), pointer :: m1d_y
  type(sll_logical_mesh_2d), pointer :: m2d_xy

  
  m1d_x => new_logical_mesh_1d(64)
  m1d_y => new_logical_mesh_1d(32,eta_min=-1.0_f64,eta_max=+1.0_f64)
  m2d_xy => m1d_x * m1d_y

  call display(m1d_x)
  call display(m1d_y)
  call display(m2d_xy)
  
  m2d => new_logical_mesh_2d(100,100)
  m3d => new_logical_mesh_3d(100,100,100)
  m4d => new_logical_mesh_4d(32,32,32,32, eta1_min=-1.0_f64, eta1_max = 2.0_f64)

  call display(m2d)

  call delete(m2d)
  call delete(m3d)
  call delete(m4d)

  print *, 'PASSED'

end program unit_test_logical_meshes
