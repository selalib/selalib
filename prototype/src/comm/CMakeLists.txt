ADD_LIBRARY(sll_comm STATIC sll_comm.F90)
TARGET_LINK_LIBRARIES(sll_comm sll_collective)

ADD_EXECUTABLE(test_comm unit_test.F90)
TARGET_LINK_LIBRARIES(test_comm
  sll_comm 
  sll_collective 
  sll_timer)

