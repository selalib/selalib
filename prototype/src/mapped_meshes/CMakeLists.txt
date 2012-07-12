IF(Fortran_COMPILER_NAME STREQUAL "xlf")
  set(CMAKE_Fortran_FLAGS "-qxlf2003=polymorphic ${CMAKE_Fortran_FLAGS}")
ENDIF()

IF(STDF95)
 add_library(sll_mapped_meshes STATIC
	sll_mapped_meshes_1d.F90
	sll_mapped_meshes_2d.F90
  )
ELSE()
 add_library(sll_mapped_meshes STATIC
	sll_mapped_meshes_1d_base.F90
	sll_mapped_meshes_2d_base.F90
	sll_mapped_meshes_1d.F90
	sll_mapped_meshes_2d.F90
	sll_mapped_meshes_2d_cartesian.F90
  )
ENDIF()
add_dependencies(sll_mapped_meshes
	sll_interpolators
	sll_splines
	geometry_functions
	sll_constants
	sll_low_level_file_io
	sll_memory
	sll_assertion
	sll_working_precision
)

add_executable(test_mapped_meshes_1d unit_test_1D.F90)
target_link_libraries(test_mapped_meshes_1d
	sll_mapped_meshes
	sll_interpolators
	sll_splines
	geometry_functions
	sll_constants
	sll_low_level_file_io
	sll_memory
	sll_assertion
	sll_working_precision
  )

add_executable(test_mapped_meshes_2d unit_test.F90)
target_link_libraries(test_mapped_meshes_2d
	sll_mapped_meshes
	sll_interpolators
	sll_splines
	geometry_functions
	sll_constants
	sll_low_level_file_io
	sll_memory
	sll_assertion
	sll_working_precision
  )

IF(NOT STDF95)
 add_executable(test_mapped_meshes_cartesian unit_test_cartesian.F90)
 target_link_libraries(test_mapped_meshes_cartesian
	sll_mapped_meshes
	sll_interpolators
	sll_splines
	geometry_functions
	sll_constants
	sll_low_level_file_io
	sll_memory
	sll_assertion
 )
ENDIF()
