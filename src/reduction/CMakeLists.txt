ADD_LIBRARY(sll_reduction STATIC sll_reduction.F90)
TARGET_LINK_LIBRARIES(sll_reduction sll_meshes) 

ADD_EXECUTABLE(test_reduction unit_test.F90)
TARGET_LINK_LIBRARIES(test_reduction sll_reduction)

# Ctest
ADD_TEST( NAME reduction COMMAND test_reduction )
SET_TESTS_PROPERTIES( reduction PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )
