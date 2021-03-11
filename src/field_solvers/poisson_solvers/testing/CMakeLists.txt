# Save .mod files locally (this does not pollute the main "modules" directory)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(test_poisson_1d_polar test_poisson_1d_polar.F90)
target_link_libraries(test_poisson_1d_polar sll_poisson_solvers)
add_test(NAME poisson_1d_polar COMMAND test_poisson_1d_polar)
set_tests_properties(poisson_1d_polar PROPERTIES PASS_REGULAR_EXPRESSION
                                                 "PASSED")

add_executable(
  test_poisson_2d_polar
  test_poisson_2d_polar.F90
  m_test_poisson_2d_polar_base.F90
  m_test_poisson_2d_polar_annulus_dirichlet.F90
  m_test_poisson_2d_polar_annulus_neumann_mode0.F90
  m_test_poisson_2d_polar_disk_dirichlet.F90)

target_link_libraries(test_poisson_2d_polar sll_poisson_solvers)
add_test(NAME poisson_2d_polar COMMAND test_poisson_2d_polar)
set_tests_properties(poisson_2d_polar PROPERTIES PASS_REGULAR_EXPRESSION
                                                 "PASSED")

add_executable(test_poisson_1d_periodic test_poisson_1d_periodic.F90)
target_link_libraries(test_poisson_1d_periodic sll_poisson_solvers)
add_test(NAME poisson_1d_periodic COMMAND test_poisson_1d_periodic)
set_tests_properties(poisson_1d_periodic PROPERTIES PASS_REGULAR_EXPRESSION
                                                    "PASSED")

add_executable(test_poisson_1d_hmf test_poisson_1d_hmf.F90)
target_link_libraries(test_poisson_1d_hmf sll_poisson_solvers)
add_test(NAME poisson_1d_hmf COMMAND test_poisson_1d_hmf)
set_tests_properties(poisson_1d_hmf PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

# Cannot be compiled ADD_EXECUTABLE( test_poisson_1d_fem test_poisson_1d_fem.F90
# ) TARGET_LINK_LIBRARIES( test_poisson_1d_fem sll_poisson_solvers )

add_executable(test_poisson_2d_periodic test_poisson_2d_periodic.F90)
target_link_libraries(test_poisson_2d_periodic sll_poisson_solvers)
add_test(NAME poisson_2d_periodic COMMAND test_poisson_2d_periodic)
set_tests_properties(poisson_2d_periodic PROPERTIES PASS_REGULAR_EXPRESSION
                                                    "PASSED")

if(FFT_LIB MATCHES FFTW)
  add_executable(test_poisson_3d_periodic test_poisson_3d_periodic.F90)
  target_link_libraries(test_poisson_3d_periodic sll_poisson_solvers sll_fft)
  add_test(NAME poisson_3d_periodic COMMAND test_poisson_3d_periodic)
  set_tests_properties(poisson_3d_periodic PROPERTIES PASS_REGULAR_EXPRESSION
                                                      "PASSED")
endif(FFT_LIB MATCHES FFTW)

add_executable(test_poisson_2d_triangular test_poisson_2d_tri.F90)
target_link_libraries(test_poisson_2d_triangular sll_poisson_solvers
                      sll_mesh_calculus sll_file_io)
add_test(NAME poisson_2d_triangular COMMAND test_poisson_2d_triangular)
set_tests_properties(poisson_2d_triangular PROPERTIES PASS_REGULAR_EXPRESSION
                                                      "PASSED")

if(BUILD_MUDPACK)

  add_executable(test_poisson_2d_mudpack_cartesian test_mudpack_cartesian.F90)
  target_link_libraries(test_poisson_2d_mudpack_cartesian sll_poisson_solvers)
  add_test(NAME poisson_2d_mudpack_cartesian
           COMMAND test_poisson_2d_mudpack_cartesian)
  set_tests_properties(poisson_2d_mudpack_cartesian
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  add_executable(test_poisson_2d_mudpack_polar test_mudpack_polar.F90)
  target_link_libraries(test_poisson_2d_mudpack_polar sll_poisson_solvers)
  add_test(NAME poisson_2d_mudpack_polar COMMAND test_poisson_2d_mudpack_polar)
  set_tests_properties(poisson_2d_mudpack_polar
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

  add_executable(test_poisson_2d_mudpack_curvilinear
                 test_poisson_2d_mudpack_curvilinear.F90)
  target_link_libraries(test_poisson_2d_mudpack_curvilinear sll_poisson_solvers)
  add_test(NAME poisson_2d_mudpack_curvilinear
           COMMAND test_poisson_2d_mudpack_curvilinear)
  set_tests_properties(poisson_2d_mudpack_curvilinear
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

endif(BUILD_MUDPACK)

add_executable(test_poisson_2d_fem test_poisson_2d_fem.F90)
target_link_libraries(test_poisson_2d_fem sll_poisson_solvers sll_memory
                      sll_utilities ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
add_test(NAME poisson_2d_fem COMMAND test_poisson_2d_fem)
set_tests_properties(poisson_2d_fem PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

add_executable(test_hex_poisson test_hex_poisson.F90)
target_link_libraries(test_hex_poisson sll_poisson_solvers)
add_test(NAME hex_poisson COMMAND test_hex_poisson)
set_tests_properties(hex_poisson PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

if(BUILD_MUDPACK)
  add_executable(test_poisson_2d_mudpack_curvilinear_old
                 test_poisson_2d_mudpack_curvilinear_old.F90)
  target_link_libraries(test_poisson_2d_mudpack_curvilinear_old
                        sll_poisson_solvers)
  add_test(NAME poisson_2d_mudpack_curvilinear_old
           COMMAND test_poisson_2d_mudpack_curvilinear_old)
  set_tests_properties(poisson_2d_mudpack_curvilinear_old
                       PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
endif(BUILD_MUDPACK)

# Test FEM Poisson solver based on smooth polar splines (dense version)
add_executable(test_poisson_2d_fem_sps_dense test_poisson_2d_fem_sps_dense.F90)
target_link_libraries(test_poisson_2d_fem_sps_dense sll_poisson_solvers)

# -------------------------------------------------------------------------------
# Test FEM Poisson solver based on smooth polar splines (stencil version)
# -------------------------------------------------------------------------------

# Copy reference output files
file(GLOB OUTPUT_DATA "*.h5")
file(COPY ${OUTPUT_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(TEST_NAME test_poisson_2d_fem_sps_stencil)

# Copy reference input file
file(GLOB INPUT_FILE ${TEST_NAME}.nml)
file(COPY ${INPUT_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(${TEST_NAME} ${TEST_NAME}.F90)
target_link_libraries(${TEST_NAME} sll_poisson_solvers)
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
set_tests_properties(${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "PASSED"
                                             TIMEOUT 60)

# compile executable to produce grid for manufactured solution
add_executable(${TEST_NAME}_grid ${TEST_NAME}_grid.F90)
target_link_libraries(${TEST_NAME}_grid sll_poisson_solvers)
