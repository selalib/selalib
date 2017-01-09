# Executables not added to ctest
IF(HDF5_IS_PARALLEL)
SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(HDF5_IS_PARALLEL)

IF( NOT PGI)
  ADD_EXECUTABLE(test_dg_fields test_dg_fields.F90)
  TARGET_LINK_LIBRARIES(test_dg_fields sll_fields)
  ADD_EXECUTABLE( test_scalar_field_1d test_scalar_field_1d.F90  )
  TARGET_LINK_LIBRARIES( test_scalar_field_1d sll_fields )
  ADD_TEST( NAME scalar_field_1d COMMAND test_scalar_field_1d )
  SET_TESTS_PROPERTIES( scalar_field_1d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
ENDIF( NOT PGI)

ADD_EXECUTABLE( test_scalar_field_2d test_scalar_field_2d.F90 helper_functions.F90 )
TARGET_LINK_LIBRARIES( test_scalar_field_2d sll_fields  )
ADD_TEST( NAME scalar_field_2d COMMAND test_scalar_field_2d )
SET_TESTS_PROPERTIES( scalar_field_2d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ADD_EXECUTABLE( test_plot_array_4d test_plot_array_4d.F90 )
TARGET_LINK_LIBRARIES( test_plot_array_4d sll_fields  )
ADD_TEST( NAME plot_array_4d COMMAND test_plot_array_4d )
SET_TESTS_PROPERTIES( plot_array_4d PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
