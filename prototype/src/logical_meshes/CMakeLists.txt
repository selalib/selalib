add_library(sll_logical_meshes STATIC
  sll_logical_meshes.F90
  )

add_dependencies( sll_logical_meshes
  sll_working_precision
  sll_memory
  )

add_executable( test_logical_meshes unit_test.F90)
target_link_libraries( test_logical_meshes
  sll_logical_meshes
  sll_working_precision
  sll_memory
  )