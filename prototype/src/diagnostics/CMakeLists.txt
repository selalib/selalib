# use all sll_*F90 files
file(GLOB DIAGNOS_SRC "sll_*.F90")
add_library(sll_diagnostics STATIC ${DIAGNOS_SRC})
add_executable(test_diagnostics unit_test.F90)
target_link_libraries(test_diagnostics sll_diagnostics sll_memory sll_assert)
