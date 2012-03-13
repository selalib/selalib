set(CONSTANTS_SRC physical_constants.F90 numeric_constants.F90)
add_library(sll_constants STATIC ${CONSTANTS_SRC})
add_executable(test_constants unit_test.F90)
target_link_libraries(test_constants sll_constants)
