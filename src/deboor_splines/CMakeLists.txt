
ADD_LIBRARY(deboor_splines STATIC
  sll_module_deboor_splines_2d.F90
  sll_module_deboor_splines_1d.F90				   )

ADD_LIBRARY(banded_linear_solver STATIC
  banfac.f90
  banslv.f90)

TARGET_LINK_LIBRARIES( deboor_splines banded_linear_solver sll_memory sll_working_precision)
  
ADD_EXECUTABLE(test_deboor_splines_1d unit_test.F90)
TARGET_LINK_LIBRARIES( test_deboor_splines_1d deboor_splines )

ADD_EXECUTABLE(test_deboor_dirichlet test_deboor_dirichlet.F90)
TARGET_LINK_LIBRARIES( test_deboor_dirichlet deboor_splines sll_memory sll_assert)
ADD_EXECUTABLE(test_deboor_hermite test_deboor_hermite.F90)
TARGET_LINK_LIBRARIES( test_deboor_hermite deboor_splines sll_memory sll_assert)
