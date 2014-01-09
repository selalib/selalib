IF(NOT STDF95)

   ADD_LIBRARY(sll_point_to_point_communications 
     STATIC sll_point_to_point_comms.F90)

   TARGET_LINK_LIBRARIES(sll_point_to_point_communications
     sll_collective)

   ADD_EXECUTABLE(test_p2p_comms_1d unit_test_1d.F90)
   TARGET_LINK_LIBRARIES(test_p2p_comms_1d
     sll_point_to_point_communications 
     sll_collective 
     sll_timer)

   ADD_EXECUTABLE(test_p2p_comms_2d unit_test_2d.F90)
   TARGET_LINK_LIBRARIES(test_p2p_comms_2d
     sll_point_to_point_communications 
     sll_collective 
     sll_timer)

ENDIF(NOT STDF95)
