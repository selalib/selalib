SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
IF(STDF95)
 add_library(sll_coordinate_transformations STATIC
        sll_coordinate_transformations_2d.F90
	sll_common_coordinate_transformations.F90
  )
ELSE()
 add_library(sll_coordinate_transformations STATIC
        sll_coordinate_transformation_2d_base.F90
        sll_coordinate_transformations_2d.F90
	sll_common_coordinate_transformations.F90
  )
ENDIF()
add_dependencies(sll_coordinate_transformations
	sll_interpolators
	sll_splines
        sll_utilities	
	sll_constants
	sll_file_io
	sll_logical_meshes
	sll_utilities
	sll_memory
	sll_assertion
	sll_working_precision
)

add_executable(test_coordinate_transformations_2d unit_test_2d.F90)
target_link_libraries(test_coordinate_transformations_2d
	sll_coordinate_transformations
	sll_interpolators
	sll_splines
        sll_utilities	
	sll_constants
	sll_file_io
	sll_logical_meshes
	sll_utilities
	sll_memory
	sll_assertion
	sll_working_precision
  )

