
ADD_LIBRARY(deboor_splines STATIC
  sll_module_deboor_splines_2d.F90
  sll_module_deboor_splines_1d.F90				   
  banfac.f90
  banslv.f90	
  splopt.f90
  knots.f90)	

TARGET_LINK_LIBRARIES( deboor_splines sll_memory)
  
ADD_EXECUTABLE(test_deboor unit_test.F90)
TARGET_LINK_LIBRARIES( test_deboor deboor_splines  sll_utilities sll_file_io)
