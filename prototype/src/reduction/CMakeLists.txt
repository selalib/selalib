IF(NOT STDF95)

   ADD_LIBRARY(sll_reduction STATIC sll_reduction.F90)
   TARGET_LINK_LIBRARIES(sll_reduction sll_logical_meshes) 

   ADD_EXECUTABLE(test_reduction unit_test.F90)
   TARGET_LINK_LIBRARIES(test_reduction sll_reduction)

ENDIF(NOT STDF95)
