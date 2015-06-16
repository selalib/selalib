
ADD_LIBRARY(deboor_splines STATIC
  sll_module_deboor_splines_2d.F90
  sll_module_deboor_splines_1d.F90				   
  splopt.f90
  knots.f90)	

ADD_LIBRARY(banded_linear_solver STATIC
  banfac.f90
  banslv.f90)

ADD_LIBRARY(original_deboor_splines STATIC
  interv.f
  bsplvd.f
  bsplvb.f)

TARGET_LINK_LIBRARIES( deboor_splines banded_linear_solver sll_memory sll_working_precision)
  
ADD_EXECUTABLE(test_deboor unit_test.F90)
TARGET_LINK_LIBRARIES( test_deboor deboor_splines sll_utilities sll_file_io)

