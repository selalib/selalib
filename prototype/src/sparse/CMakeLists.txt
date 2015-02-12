ADD_LIBRARY(sll_sparse STATIC sll_sparse.F90)

TARGET_LINK_LIBRARIES(sll_sparse sll_utilities sll_memory)

ADD_EXECUTABLE(test_sparse unit_test_sparse.F90)

TARGET_LINK_LIBRARIES(test_sparse sll_sparse)
