ADD_LIBRARY(sll_integration STATIC gauss_legendre.F90
                                   gauss_lobatto.F90 )

TARGET_LINK_LIBRARIES(sll_integration sll_splines)

ADD_EXECUTABLE(test_integration unit_test.F90 test_func_module.F90)

TARGET_LINK_LIBRARIES(test_integration sll_integration)
