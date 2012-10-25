ADD_LIBRARY(sll_splines STATIC sll_splines.F90
                               sll_quintic_splines.F90
                               sll_odd_degree_splines.F90
                               sll_arbitrary_degree_splines.F90)

ADD_DEPENDENCIES(sll_splines sll_memory 
                             sll_utilities
                             sll_assertion 
                             sll_working_precision
                             sll_numerical_utilities)


SET(SPL_LIBS sll_splines 
             sll_numerical_utilities 
             sll_memory 
             sll_assertion 
             sll_constants 
             ${LAPACK_LIBRARIES})


#uniform
ADD_EXECUTABLE(test_splines unit_test.F90  
                            test_func_module.F90 
                            test_processes_module.F90 )

TARGET_LINK_LIBRARIES(test_splines ${SPL_LIBS})

#non uniform

ADD_EXECUTABLE(test_non_unif_splines unit_test_nonuniform.F90 
                                     cubic_nonuniform_splines.F90)

TARGET_LINK_LIBRARIES(test_non_unif_splines ${SPL_LIBS})

#arbitrary degree

ADD_EXECUTABLE(test_arbitrary_degree_splines unit_test_arbitrary_degree.F90)

TARGET_LINK_LIBRARIES(test_arbitrary_degree_splines ${SPL_LIBS})

#quintic
ADD_EXECUTABLE(test_quintic_splines test_sll_quintic_splines.F90)
TARGET_LINK_LIBRARIES(test_quintic_splines ${SPL_LIBS})

#odd degree
ADD_EXECUTABLE(test_odd_degree_splines test_sll_odd_degree_splines.F90 
                                       sll_odd_degree_splines.F90)

TARGET_LINK_LIBRARIES(test_odd_degree_splines ${SPL_LIBS})

