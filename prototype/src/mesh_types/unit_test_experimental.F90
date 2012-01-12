program unit_test_experimental
  use sll_mesh_types_experimental
  use geometry_functions
  implicit none

  type(mapping_2D), pointer     :: map
  type(mesh_2D_scalar), pointer :: m2d

  map => new_mapping_2D( 256, 256, ANALYTIC_MAP )
  call initialize_mapping_2D( &
       map, &
       polar_jac11, &
       polar_jac12, &
       polar_jac21, &
       polar_jac22, &
       polar_x1, &
       polar_x2 )

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
