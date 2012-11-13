IF(Fortran_COMPILER_NAME STREQUAL "xlf")
  set(CMAKE_Fortran_FLAGS "-qxlf2003=polymorphic ${CMAKE_Fortran_FLAGS}")
ENDIF()

add_library(sll_distribution_function STATIC
	sll_distribution_function.F90
	initial_functions.F90)
add_dependencies(sll_distribution_function
	sll_working_precision
	sll_constants
	sll_memory
	sll_assertion
	sll_scalar_field_2d
)

IF (NOT STDF95)
add_executable(test_distribution_function unit_test.F90)
target_link_libraries(test_distribution_function
	sll_distribution_function
	geometry_functions
	sll_scalar_field_2d
	sll_low_level_file_io
	sll_interpolators
	sll_splines
	sll_mapped_meshes
	sll_utilities
	sll_memory
	sll_assertion
	sll_constants
	sll_working_precision
)
ENDIF (NOT STDF95)
