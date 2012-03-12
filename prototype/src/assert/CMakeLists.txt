set(ASSERT_SRC sll_assert.F90)
add_library(sll_assert STATIC ${ASSERT_SRC})
add_executable(test_assert unit_test.F90)
target_link_libraries(test_assert sll_assert sll_memory)
