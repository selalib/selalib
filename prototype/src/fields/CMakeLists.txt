add_library(sll_scalar_field_2d STATIC sll_scalar_field_2d.F90)
add_library(sll_scalar_field_1d STATIC sll_scalar_field_1d.F90)

add_executable(test_scalar_field unit_test.F90)
target_link_libraries(test_scalar_field
	sll_scalar_field_2d
	sll_mapped_meshes
	sll_splines
	sll_low_level_file_io
	sll_constants
	sll_utilities
	geometry_functions
	sll_assertion
	sll_memory
	sll_working_precision
)