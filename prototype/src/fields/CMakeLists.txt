ADD_LIBRARY(geometry_functions STATIC geometry_functions.F90)

ADD_DEPENDENCIES( geometry_functions sll_assertion
	                             sll_memory
	                             sll_working_precision)

IF(NOT STDF95)
add_library(sll_scalar_field_2d STATIC 
  sll_scalar_field_2d_base.F90
  sll_scalar_field_2d.F90 
  sll_scalar_field_2d_alternative.F90
  scalar_field_initializer_base.F90
  landau_2d_initializer.F90
  gaussian_2d_initializer.F90
  tsi_2d_initializer.F90	
	)		

add_dependencies(sll_scalar_field_2d
	sll_mapped_meshes
	sll_interpolators
	sll_splines
	sll_file_io
	sll_constants
	sll_utilities
	geometry_functions
	sll_assertion
	sll_memory
	sll_working_precision
)

add_library(sll_scalar_field_1d STATIC sll_scalar_field_1d.F90
		                       scalar_field_initializer_base.F90)
add_dependencies(sll_scalar_field_1d
	sll_mapped_meshes
	sll_interpolators
	sll_splines
	sll_file_io
	sll_constants
	sll_utilities
	geometry_functions
	sll_assertion
	sll_memory
	sll_working_precision
)

add_library(sll_scalar_field STATIC scalar_field_initializer_base.F90 
  distribution_function_initializer_4d.F90
  distribution_function_initializer_6d.F90
)

add_dependencies( sll_scalar_field
                  sll_remap
                  sll_constants
		  sll_utilities
		  sll_assertion
		  sll_memory
		  sll_working_precision
)

add_executable(test_scalar_field unit_test.F90)
target_link_libraries(test_scalar_field
	sll_scalar_field_2d
	sll_mapped_meshes
	sll_interpolators
	sll_splines
	sll_file_io
	sll_constants
	sll_utilities
	geometry_functions
	sll_assertion
	sll_memory
	sll_working_precision
)
ENDIF(NOT STDF95)
