IF(FORTRANCL_ENABLED)

   FIND_PACKAGE(FORTRANCL QUIET)

   IF(FORTRANCL_FOUND)

      MESSAGE(STATUS "FORTRANCL FOUND")
      INCLUDE_DIRECTORIES(${FORTRANCL_INCLUDE_DIRS})

      add_executable(test_opencl sum.F90)
      target_link_libraries(test_opencl ${FORTRANCL_LIBRARIES} ${OPENCL_LIBRARIES})
      add_executable(vlasov_gpu tridiag.F90)
      target_link_libraries(vlasov_gpu ${FORTRANCL_LIBRARIES} ${OPENCL_LIBRARIES})

   ENDIF(FORTRANCL_FOUND)

ENDIF(FORTRANCL_ENABLED)
