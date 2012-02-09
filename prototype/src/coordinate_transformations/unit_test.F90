program unit_test_experimental
  use sll_mapped_mesh
  use geometry_functions
  implicit none

  type(map_2D), pointer     :: map
  type(mapped_mesh_2D_scalar), pointer :: m2d

  map => new_map_2D( ANALYTIC_MAP )
  call initialize_map_2D( &
       map, &
       256, &
       256, &
       x1_func=polar_x1, &
       x2_func=polar_x2, &
       j11_func=polar_jac11, &
       j12_func=polar_jac12, &
       j21_func=polar_jac21, &
       j22_func=polar_jac22 )

  m2d => new_mesh_2D_scalar( &
       NODE_CENTERED_MESH, &
       PERIODIC_MESH_BC, &
       HERMITE_MESH_BC, &
       map )

  print *, get_m2ds_node(m2d,16,16)
  call set_m2ds_node(m2d,16,16,5.0_f64)
  print *, get_m2ds_node(m2d,16,16)

  call delete(m2d)
  call delete(map)

end program unit_test_experimental
