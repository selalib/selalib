add_library(sll_scalar_field_2d STATIC sll_scalar_field_2d.F90)

add_executable(test_scalar_field unit_test.F90)
target_link_libraries(test_scalar_field
	sll_mapped_mesh
	sll_low_level_io
	sll_splines
	geometry_functions
	sll_constants
	sll_utilities
	sll_assertion
	sll_precision
	sll_memory)

