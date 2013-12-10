ADD_LIBRARY(sll_collective STATIC sll_collective.F90)
TARGET_LINK_LIBRARIES(sll_collective sll_utilities)

ADD_EXECUTABLE(test_collective unit_test.F90)
TARGET_LINK_LIBRARIES(test_collective sll_collective sll_timer)

