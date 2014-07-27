#ADD_LIBRARY(deboor_splines STATIC
  #dvalue2d.f90
  #dvalue1d.f90
  #bvalue2d.f90
  #spli2d_perper.f90
  #spli2d_per.f90
  #spli1d_per.f90	
  #spli1d_dir.f90
  #spli2d_custom.f90
  #splint.f90
 # splopt.f90
 # knots.f90
  #spli2d.f90
 # banfac.f90
 # banslv.f90
  #bsplvd.f90
  #bsplvb.f90
  #bvalue.f90
 # interv.f90
#)


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
