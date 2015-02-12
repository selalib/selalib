ADD_LIBRARY(sll_assertion STATIC sll_assert.F90)
ADD_EXECUTABLE(test_assert unit_test.F90)
TARGET_LINK_LIBRARIES(test_assert sll_assertion) 
