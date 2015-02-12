ADD_EXECUTABLE(test_newton_raphson newton_raphson.F90
                                   test_jacobian_and_f_module.F90 
                                   unit_test.F90 )

TARGET_LINK_LIBRARIES(test_newton_raphson sll_working_precision)
