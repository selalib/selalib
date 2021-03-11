# Save .mod files locally (this does not pollute the main "modules" directory)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(
  test_qn_solver_2d_polar
  test_qn_solver_2d_polar.F90
  m_test_qn_solver_2d_polar_base.F90
  m_test_qn_solver_2d_polar_annulus_dirichlet.F90
  m_test_qn_solver_2d_polar_annulus_neumann_mode0.F90
  m_test_qn_solver_2d_polar_disk_dirichlet.F90)

target_link_libraries(test_qn_solver_2d_polar sll_quasi_neutral_solvers)

add_test(NAME qn_solver_2d_polar COMMAND test_qn_solver_2d_polar)
set_tests_properties(qn_solver_2d_polar PROPERTIES PASS_REGULAR_EXPRESSION
                                                   "PASSED" TIMEOUT 10)

# -------------------------------------------------------------------------------
# Test FEM quasi-neutrality solver based on smooth polar splines (stencil
# version)
# -------------------------------------------------------------------------------

# Copy reference output files
file(GLOB OUTPUT_DATA "*.h5")
file(COPY ${OUTPUT_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(TEST_NAME test_qn_solver_2d_fem_sps_stencil)

# Copy reference input file
file(GLOB INPUT_FILE ${TEST_NAME}.nml)
file(COPY ${INPUT_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(${TEST_NAME} ${TEST_NAME}.F90)
target_link_libraries(${TEST_NAME} sll_quasi_neutral_solvers)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
set_tests_properties(${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
                                             TIMEOUT 60)

# compile executable to produce grid for manufactured solution
add_executable(${TEST_NAME}_grid ${TEST_NAME}_grid.F90)
target_link_libraries(${TEST_NAME}_grid sll_quasi_neutral_solvers)
