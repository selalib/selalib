IF(Fortran_COMPILER_NAME STREQUAL "xlf")
  set(CMAKE_Fortran_FLAGS "-qxlf2003=polymorphic ${CMAKE_Fortran_FLAGS}")
ENDIF()

 add_library(sll_time_splitting STATIC
	sll_time_splitting_base.F90
	sll_const_coef_adv_2d.F90
  )

add_dependencies(sll_time_splitting
	sll_interpolators
	sll_splines
	sll_low_level_file_io
	sll_memory
	sll_assertion
	sll_working_precision
)

add_executable(test_time_splitting unit_test.F90)
target_link_libraries(test_time_splitting
        sll_time_splitting
	sll_interpolators
	sll_low_level_file_io
	sll_splines
	sll_memory
	sll_assertion
	sll_working_precision
  )
