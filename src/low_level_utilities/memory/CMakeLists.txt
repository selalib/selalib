set(MEMORY_SRC sll_m_memory.F90)
add_library(sll_memory STATIC ${MEMORY_SRC})
target_link_libraries(sll_memory sll_working_precision sll_constants)

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
