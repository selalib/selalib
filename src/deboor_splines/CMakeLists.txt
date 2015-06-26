
ADD_LIBRARY( deboor_splines STATIC
             sll_bsplines.F90
             sll_module_deboor_splines_2d.F90
             sll_module_deboor_splines_1d.F90 )

ADD_LIBRARY(banded_linear_solver STATIC banfac.f90 banslv.f90)

TARGET_LINK_LIBRARIES( deboor_splines banded_linear_solver 
                       sll_utilities sll_memory sll_working_precision sll_assert)
  
IF(BUILD_TESTING)

  ADD_EXECUTABLE(test_deboor_splines_1d unit_test.F90)
  TARGET_LINK_LIBRARIES( test_deboor_splines_1d deboor_splines )

  ADD_TEST( NAME deboor_splines COMMAND test_deboor_splines_1d   )
  SET_TESTS_PROPERTIES( deboor_splines PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ENDIF (BUILD_TESTING)
