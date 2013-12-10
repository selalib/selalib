ADD_LIBRARY(sll_timer STATIC sll_timer.F90)
ADD_EXECUTABLE(test_timer unit_test.F90)
TARGET_LINK_LIBRARIES(test_timer sll_timer)
