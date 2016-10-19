SET(CLAPP_LIBRARIES fema glt disco spl plaf clappio)

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include/spl)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include/fema)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include/disco)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include/plaf)
ADD_EXECUTABLE(test_fema_poisson_2d_periodic test_poisson_2d_periodic.F90)
ADD_DEPENDENCIES(test_fema_poisson_2d_periodic fema_lib)
TARGET_LINK_LIBRARIES(test_fema_poisson_2d_periodic ${CLAPP_LIBRARIES} )

SET(PARAMETERS_FILE "parameters_sim_periodic_1d1d.nml")
SET(TENSOR_FILE     "parameters_tensor_2d.nml")
SET(BSPLINE_FILE    "parameters_bspline_periodic.nml")
SET(SOLVER_FILE1    "solver_driver.nml" )
SET(SOLVER_FILE2    "solver_mgmres_nullspace.nml" )
FILE(COPY ${PARAMETERS_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY ${BSPLINE_FILE}    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY ${TENSOR_FILE}     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY ${SOLVER_FILE1}    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY ${SOLVER_FILE2}    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

ADD_TEST(NAME test_fema_poisson_2d_periodic 
	COMMAND test_fema_poisson_2d_periodic ${PARAMETERS_FILE} ${SOLVER_FILE1})
