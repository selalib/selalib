add_library(sll_simulation_1d STATIC
	simulation_1d.F90
  )

add_executable(test_simulation_1d user_geometry_functions.F90 unit_test.F90)

target_link_libraries(test_simulation_1d
	sll_mapped_meshes
	interpolators
	sll_splines
	geometry_functions
	sll_distribution_function
	sll_scalar_field_2d	
	sll_constants
	sll_low_level_file_io
        sll_utilities
	sll_memory
	sll_assertion
  )

