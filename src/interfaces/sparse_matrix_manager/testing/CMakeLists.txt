# PN sll_collective is necessary when PasTiX is used

add_executable(test_sparse_matrix_cg test_sparse_matrix_cg.F90)
target_link_libraries(test_sparse_matrix_cg sll_sparse_matrix)
add_test(NAME sparse_matrix_cg COMMAND test_sparse_matrix_cg)
set_tests_properties(sparse_matrix_cg PROPERTIES PASS_REGULAR_EXPRESSION
                                                 "PASSED" TIMEOUT 20)

add_executable(test_lapack test_lapack.F90)
target_link_libraries(test_lapack sll_memory ${LAPACK_LIBRARIES})

if(SUITESPARSE_FOUND)

  add_executable(test_sparse_matrix_umfpack test_sparse_matrix_umfpack.F90)
  target_link_libraries(test_sparse_matrix_umfpack sll_sparse_matrix)
  add_test(NAME sparse_matrix_umfpack COMMAND test_sparse_matrix_umfpack)
  set_tests_properties(sparse_matrix_umfpack PROPERTIES PASS_REGULAR_EXPRESSION
                                                        "PASSED" TIMEOUT 20)

endif(SUITESPARSE_FOUND)

if(PASTIX_FOUND)

  add_executable(test_sparse_matrix_pastix test_sparse_matrix_pastix.F90)
  target_link_libraries(test_sparse_matrix_pastix sll_sparse_matrix)
  add_test(NAME sparse_matrix_pastix COMMAND test_sparse_matrix_pastix)
  set_tests_properties(sparse_matrix_pastix PROPERTIES PASS_REGULAR_EXPRESSION
                                                       "PASSED" TIMEOUT 20)

endif(PASTIX_FOUND)

if(MUMPS_FOUND)

  add_executable(test_sparse_matrix_mumps test_sparse_matrix_mumps.F90)
  target_link_libraries(test_sparse_matrix_mumps sll_sparse_matrix)
  add_test(NAME sparse_matrix_mumps COMMAND test_sparse_matrix_mumps)
  set_tests_properties(sparse_matrix_mumps PROPERTIES PASS_REGULAR_EXPRESSION
                                                      "PASSED" TIMEOUT 20)

endif(MUMPS_FOUND)

add_library(sll_csr_one_constraint STATIC sll_csr_one_constraint_module.F90)
target_link_libraries(sll_csr_one_constraint sll_utilities sll_memory)
add_executable(test_csr_one_constraint test_csr_one_constraint.F90)
target_link_libraries(test_csr_one_constraint sll_csr_one_constraint)

if(Fortran_COMPILER MATCHES IBM)
  set_source_files_properties(dPackgmres.f PROPERTIES COMPILE_FLAGS -qfixed)
endif(Fortran_COMPILER MATCHES IBM)

add_executable(test_gmres test_gmres_driver.F90 dPackgmres.f)
target_link_libraries(test_gmres ${BLAS_LIBRARIES} sll_constants)
add_test(NAME gmres COMMAND test_gmres)
set_tests_properties(gmres PROPERTIES TIMEOUT 2)

if(BUILD_PARALUTION)

  link_directories(
    ${CMAKE_BINARY_DIR}/paralution-prefix/src/paralution-build/lib/)
  add_executable(test_paralution test_paralution.F90)
  target_link_libraries(test_paralution sll_sparse_matrix sll_paralution
                        paralution)
  add_test(NAME paralution-${PARALUTION_VERSION} COMMAND test_paralution)
  set_tests_properties(paralution-${PARALUTION_VERSION} PROPERTIES TIMEOUT 10)

endif(BUILD_PARALUTION)

if(PASTIX_FOUND) # Simple test from the library

  add_executable(test_pastix test_pastix.F90)
  target_link_libraries(
    test_pastix
    sll_sparse_matrix
    sll_collective
    sll_memory
    ${PASTIX_LIBRARIES}
    ${PTSCOTCH_LIBRARIES}
    ${SCOTCH_LIBRARIES}
    ${BLAS_LIBRARIES})

endif(PASTIX_FOUND)

if(MUMPS_FOUND) # Simple test from the library

  add_executable(test_mumps test_mumps.F90)
  target_link_libraries(test_mumps sll_sparse_matrix sll_collective
                        ${MUMPS_LIBRARIES})
  add_test(NAME mumps COMMAND test_mumps)
  set_tests_properties(mumps PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")

endif(MUMPS_FOUND)
