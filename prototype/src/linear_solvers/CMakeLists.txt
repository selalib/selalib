ADD_LIBRARY(sll_linear_solvers STATIC sll_linear_solvers.F90)

TARGET_LINK_LIBRARIES(sll_linear_solvers sll_memory)



IF(PASTIX_ENABLED)
   ADD_DEFINITIONS(-DPASTIX)
   ADD_EXECUTABLE(test_linear_solvers unit_test_linear_solvers.F90)
   TARGET_LINK_LIBRARIES(test_linear_solvers sll_linear_solvers 
                                             sll_pastix 
                                             sll_collective
                                             sll_sparse)
ELSE()
   ADD_EXECUTABLE(test_linear_solvers unit_test_linear_solvers.F90)
   TARGET_LINK_LIBRARIES(test_linear_solvers sll_linear_solvers 
                                             sll_sparse)
ENDIF()
