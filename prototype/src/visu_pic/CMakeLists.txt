ADD_LIBRARY(sll_visu_pic STATIC sll_visu_pic.F90)

TARGET_LINK_LIBRARIES(sll_visu_pic
  sll_utilities
  sll_file_io_parallel
  sll_assertion
  sll_working_precision )

ADD_EXECUTABLE(test_visu_pic 
  unit_test.F90 
  biot_savart.F90)

TARGET_LINK_LIBRARIES(test_visu_pic 
  sll_visu_pic)
