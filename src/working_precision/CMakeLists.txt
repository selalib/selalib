SET(PRECISION_SRC sll_working_precision.F90)
ADD_LIBRARY(sll_working_precision STATIC ${PRECISION_SRC})
ADD_EXECUTABLE(test_precision unit_test.F90)
TARGET_LINK_LIBRARIES(test_precision sll_working_precision)
