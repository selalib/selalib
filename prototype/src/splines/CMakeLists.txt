ADD_LIBRARY( sll_splines STATIC sll_cubic_splines.F90
                                sll_quintic_splines.F90
                                sll_odd_degree_splines.F90
                                sll_arbitrary_degree_splines.F90
                                cubic_nonuniform_splines.F90)

TARGET_LINK_LIBRARIES(sll_splines sll_utilities 
                                  sll_boundary_condition_descriptors
                                  ${LAPACK_LIBRARIES})


ADD_EXECUTABLE(test_splines unit_test.F90  
                            test_func_module.F90 
                            util_constants.F90
                            test_processes_module.F90 )

TARGET_LINK_LIBRARIES(test_splines sll_splines)

ADD_EXECUTABLE(test_non_unif_splines unit_test_nonuniform.F90 )
TARGET_LINK_LIBRARIES(test_non_unif_splines sll_splines)

ADD_EXECUTABLE(test_arbitrary_degree_splines unit_test_arbitrary_degree.F90)
TARGET_LINK_LIBRARIES(test_arbitrary_degree_splines sll_splines)

ADD_EXECUTABLE(test_quintic_splines test_sll_quintic_splines.F90)
TARGET_LINK_LIBRARIES(test_quintic_splines sll_splines)

ADD_EXECUTABLE(test_odd_degree_splines test_sll_odd_degree_splines.F90)
TARGET_LINK_LIBRARIES(test_odd_degree_splines sll_splines)

