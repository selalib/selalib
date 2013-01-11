IF(Fortran_COMPILER_NAME STREQUAL "xlf")
  set(CMAKE_Fortran_FLAGS "-qxlf2003=polymorphic ${CMAKE_Fortran_FLAGS}")
ENDIF()

add_library(sll_linrood STATIC
	sll_linrood.F90)
IF (NOT STDF95)
add_executable(test_linrood unit_test.F90)
target_link_libraries(test_linrood
        sll_linrood
        WENO_interp
        WENO_recon
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
	sll_constants)
ENDIF (NOT STDF95)
