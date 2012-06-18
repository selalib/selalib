add_library(sll_advection_field STATIC
	sll_advection_field.F90)
add_executable(test_advection_field unit_test.F90)
target_link_libraries(test_advection_field
	sll_advection_field
	sll_scalar_field_2d
        sll_scalar_field_1d
	sll_low_level_file_io
	sll_splines
	geometry_functions
	sll_mapped_meshes
	sll_utilities
	sll_memory
	sll_assertion
	sll_constants)
