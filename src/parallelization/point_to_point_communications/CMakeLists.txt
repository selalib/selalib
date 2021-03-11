add_library(sll_point_to_point_communications STATIC
            sll_m_point_to_point_comms.F90)

target_link_libraries(sll_point_to_point_communications sll_collective)

# Ctests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
