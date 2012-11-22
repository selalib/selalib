IF(Fortran_COMPILER_NAME STREQUAL "xlf")
  set(CMAKE_Fortran_FLAGS "-qxlf2003=polymorphic ${CMAKE_Fortran_FLAGS}")
ENDIF()

IF(STDF95)
add_library(sll_time_splitting STATIC
	sll_vlasov_poisson_adv_cart2d.F90
	sll_const_coef_adv_2d.F90
  )

ELSE()
add_library(sll_time_splitting STATIC
	sll_time_splitting_base.F90
	sll_const_coef_adv_2d.F90
	sll_vlasov_poisson_adv_cart2d.F90
)
ENDIF()

add_dependencies(sll_time_splitting
        sll_distribution_function
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

IF(FFTPACK_ENABLED)
  add_executable(vp1d_bsl_time_split vp1d_bsl_time_split.F90)
  target_link_libraries(vp1d_bsl_time_split
    sll_mapped_meshes
    sll_interpolators
    sll_splines
    geometry_functions
    sll_distribution_function
    sll_scalar_field_2d
    sll_time_splitting
    sll_constants
    sll_low_level_file_io
    sll_memory
    sll_assertion
    sll_working_precision
    sll_utilities
    sll_poisson_1d_periodic
    fftpack_module
  )
ENDIF()
