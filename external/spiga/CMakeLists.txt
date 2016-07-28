FILE(GLOB SOLVER_FILES ${CMAKE_CURRENT_SOURCE_DIR} "*.nml")
FILE(COPY ${SOLVER_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
ADD_EXECUTABLE(test_spiga_poisson_2d_periodic test_poisson_2d_periodic.f90)
TARGET_LINK_LIBRARIES(test_spiga_poisson_2d_periodic
spi_poisson 
spi_pde 
spi_fem 
spi_field 
spi_space 
spi_basis 
spi_mapping
spi_mesh 
spi_parameters 
spi_global
dsc_field 
dsc_space 
dsc_fem 
dsc_basis 
dsc_mesh 
dsc_quadratures 
dsc_global 
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
plf_lapack 
spi_utilities 
spi_space 
spi_numbering 
spi_bspline 
spi_mesh 
spi_global 
dsc_field 
dsc_space 
dsc_fem 
dsc_basis 
dsc_mesh 
dsc_quadratures 
dsc_global 
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
plf_lapack
)

SET(PARAMETERS_FILE "${CMAKE_CURRENT_BINARY_DIR}/data/parameters_2d.nml")
SET(SOLVER_FILE    "solver_mgmres_nullspace.nml" )

ADD_TEST(NAME test_spiga_poisson_2d_periodic 
	COMMAND test_spiga_poisson_2d_periodic ${PARAMETERS_FILE} ${SOLVER_FILE})
