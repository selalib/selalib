ADD_LIBRARY(deboor_splines STATIC
  sll_module_deboor_splines_2d.F90
  sll_module_deboor_splines_1d.F90				   
  banded_linear_system_solver.F90
  knots.f90)	

TARGET_LINK_LIBRARIES( deboor_splines sll_memory sll_utilities )
  
ADD_EXECUTABLE(test_deboor unit_test.F90)
TARGET_LINK_LIBRARIES( test_deboor deboor_splines  sll_utilities sll_file_io)

# Ctest
ADD_TEST( NAME deboor_spline COMMAND test_deboor )
SET_TESTS_PROPERTIES( deboor_spline PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
