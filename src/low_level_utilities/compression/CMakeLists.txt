add_library(sll_compression STATIC zfp_wrapper.c sll_m_compression.F90)
set(COMPRESSION_DEPS sll_working_precision)

if(ZFP_LIBRARIES)
  list(APPEND COMPRESSION_DEPS ${ZFP_LIBRARIES})
endif()

target_link_libraries(sll_compression ${COMPRESSION_DEPS})

if(ZFP_LIBRARIES)
  if(BUILD_TESTING)
    add_subdirectory(testing)
  endif(BUILD_TESTING)
endif()
