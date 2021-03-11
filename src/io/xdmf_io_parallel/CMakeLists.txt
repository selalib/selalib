set(LIB sll_xdmf_io_parallel)
set(SRC sll_xdmf_io_parallel.F90 sll_m_xdmf_light_parallel.F90)
set(DEP sll_xdmf_io sll_collective)

add_library(${LIB} STATIC ${SRC})
target_link_libraries(${LIB} ${DEP})

# Tests
if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
