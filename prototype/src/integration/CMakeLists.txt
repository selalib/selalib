ADD_LIBRARY(sll_integration STATIC gauss_legendre.F90
                                   gauss-lobatto.F90 lob.f gauss.f)

ADD_DEPENDENCIES(sll_integration sll_splines 
                                 sll_assertion 
                                 sll_working_precision)

ADD_EXECUTABLE(test_integration unit_test.F90 test_func_module.F90)

TARGET_LINK_LIBRARIES(test_integration sll_integration sll_working_precision sll_memory sll_assertion sll_constants )
