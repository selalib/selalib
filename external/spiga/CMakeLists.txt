
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
spi_global)

SET(DISCO_LIBRARIES
dsc_field 
dsc_space 
dsc_fem 
dsc_basis 
dsc_mesh 
dsc_quadratures 
dsc_global )

SET(PLAF_LIBRARIES
plf_vector 
plf_linear_solver 
plf_preconditioner 
plf_linear_operator 
plf_matrix 
plf_matrix_market 
plf_ddm 
plf_numbering 
plf_utilities 
plf_global 
plf_graph 
plf_solvers 
plf_sparsekit 
plf_lapack )

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
ADD_EXECUTABLE(test_spiga_poisson_2d_periodic test_poisson_2d_periodic.F90)
ADD_DEPENDENCIES(test_spiga_poisson_2d_periodic spiga)
TARGET_LINK_LIBRARIES(test_spiga_poisson_2d_periodic
	${SPIGA_LIBRARIES} ${DISCO_LIBRARIES} ${PLAF_LIBRARIES})

SET(PARAMETERS_FILE "parameters_2d.nml")
SET(SOLVER_FILE     "solver_mgmres_nullspace.nml" )
FILE(COPY ${PARAMETERS_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY ${SOLVER_FILE}     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

ADD_TEST(NAME test_spiga_poisson_2d_periodic 
	COMMAND test_spiga_poisson_2d_periodic ${PARAMETERS_FILE} ${SOLVER_FILE})
