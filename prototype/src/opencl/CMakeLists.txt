IF(FORTRANCL_FOUND)
   add_executable(test_opencl sum.F90)
   target_link_libraries(test_opencl ${FORTRANCL_LIBRARIES} ${OPENCL_LIBRARIES})
ENDIF()
