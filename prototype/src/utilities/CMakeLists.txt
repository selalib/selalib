set(UTILITIES_SRC misc_utils.F90)
add_library(sll_utilities STATIC ${UTILITIES_SRC})
add_executable(test_utilities unit_test.F90)
target_link_libraries(test_utilities sll_utilities)
