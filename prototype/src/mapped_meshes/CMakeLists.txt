add_library(sll_mapped_meshes STATIC
	sll_mapped_mesh_base.F90
	sll_mapped_meshes.F90
  )

add_executable(test_mapped_meshes unit_test.F90)
target_link_libraries(test_mapped_meshes
	sll_mapped_meshes
	interpolators
	sll_splines
	geometry_funcs
	sll_constants
	sll_low_level_file_io
	sll_memory
	sll_assertion
  )
