set(REMAP_SRC sll_remap.F90)
add_library(sll_remap STATIC ${REMAP_SRC})
add_executable(test_remap unit_test.F90)
target_link_libraries(test_remap sll_remap)
