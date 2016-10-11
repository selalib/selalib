
SET(SPIGA_LIBRARIES
spi_mapping
spi_poisson 
spi_bspline 
spi_pde 
spi_fem 
spi_numbering 
spi_field 
spi_space 
spi_basis 
spi_mesh 
spi_utilities 
spi_parameters 
spi_global
disco
spl
plaf
clappio)

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include/disco)
ADD_EXECUTABLE(test_spiga_poisson_2d_periodic test_poisson_2d_periodic.F90)
ADD_DEPENDENCIES(test_spiga_poisson_2d_periodic spiga_lib)
TARGET_LINK_LIBRARIES(test_spiga_poisson_2d_periodic ${SPIGA_LIBRARIES} )

SET(PARAMETERS_FILE "parameters_2d.nml")
SET(SOLVER_FILE1    "solver_driver.nml" )
SET(SOLVER_FILE2    "solver_mgmres_nullspace.nml" )
FILE(COPY ${PARAMETERS_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY ${SOLVER_FILE1}    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY ${SOLVER_FILE2}    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

ADD_TEST(NAME test_spiga_poisson_2d_periodic 
	COMMAND test_spiga_poisson_2d_periodic ${PARAMETERS_FILE} ${SOLVER_FILE1})
