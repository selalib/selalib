#PN sll_collective is necessary when PasTiX is used

IF(PASTIX_FOUND OR MUMPS_FOUND)
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
ENDIF(PASTIX_FOUND OR MUMPS_FOUND)

ADD_EXECUTABLE( test_sparse_matrix_cg test_sparse_matrix_cg.F90 )
TARGET_LINK_LIBRARIES( test_sparse_matrix_cg sll_sparse_matrix )
ADD_TEST(NAME sparse_matrix_cg COMMAND test_sparse_matrix_cg )
SET_TESTS_PROPERTIES( sparse_matrix_cg PROPERTIES 
  PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20 )

ADD_EXECUTABLE(test_lapack test_lapack.F90)
TARGET_LINK_LIBRARIES( test_lapack ${LAPACK_LIBRARIES} )

IF(SUITESPARSE_FOUND)

  ADD_EXECUTABLE( test_sparse_matrix_umfpack test_sparse_matrix_umfpack.F90 )
  TARGET_LINK_LIBRARIES( test_sparse_matrix_umfpack sll_sparse_matrix )
  ADD_TEST(NAME sparse_matrix_umfpack COMMAND test_sparse_matrix_umfpack )
  SET_TESTS_PROPERTIES( sparse_matrix_umfpack PROPERTIES 
    PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20 )

ENDIF(SUITESPARSE_FOUND)

IF(PASTIX_FOUND)

  ADD_EXECUTABLE( test_sparse_matrix_pastix test_sparse_matrix_pastix.F90 )
  TARGET_LINK_LIBRARIES( test_sparse_matrix_pastix sll_sparse_matrix )
  ADD_TEST(NAME sparse_matrix_pastix COMMAND test_sparse_matrix_pastix )
  SET_TESTS_PROPERTIES( sparse_matrix_pastix PROPERTIES 
    PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20 )

ENDIF(PASTIX_FOUND)

IF(MUMPS_FOUND)

  ADD_EXECUTABLE( test_sparse_matrix_mumps test_sparse_matrix_mumps.F90 )
  TARGET_LINK_LIBRARIES( test_sparse_matrix_mumps sll_sparse_matrix )
  ADD_TEST(NAME sparse_matrix_mumps COMMAND test_sparse_matrix_mumps )
  SET_TESTS_PROPERTIES( sparse_matrix_mumps PROPERTIES 
    PASS_REGULAR_EXPRESSION "PASSED" TIMEOUT 20 )

ENDIF(MUMPS_FOUND)


ADD_LIBRARY(sll_csr_one_constraint STATIC sll_csr_one_constraint_module.F90)
TARGET_LINK_LIBRARIES(sll_csr_one_constraint sll_utilities sll_memory)
ADD_EXECUTABLE( test_csr_one_constraint test_csr_one_constraint.F90 )
TARGET_LINK_LIBRARIES( test_csr_one_constraint sll_csr_one_constraint )

IF(Fortran_COMPILER MATCHES IBM)
 SET_SOURCE_FILES_PROPERTIES(dPackgmres.f PROPERTIES COMPILE_FLAGS -qfixed)
ENDIF(Fortran_COMPILER MATCHES IBM)

ADD_EXECUTABLE( test_gmres test_gmres_driver.F90 dPackgmres.f )
TARGET_LINK_LIBRARIES( test_gmres ${BLAS_LIBRARIES} )
ADD_TEST( NAME gmres COMMAND test_gmres )
SET_TESTS_PROPERTIES( gmres PROPERTIES TIMEOUT 2 )

IF(BUILD_PARALUTION)

  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/paralution-prefix/src/paralution-build/lib/)
  ADD_EXECUTABLE(test_paralution test_paralution.F90) 
  TARGET_LINK_LIBRARIES(test_paralution sll_sparse_matrix sll_paralution paralution)
  ADD_TEST( NAME paralution-${PARALUTION_VERSION} COMMAND test_paralution )
  SET_TESTS_PROPERTIES( paralution-${PARALUTION_VERSION} PROPERTIES TIMEOUT 10 )

ENDIF(BUILD_PARALUTION)

IF(PASTIX_FOUND) #Simple test from the library 
  
  ADD_EXECUTABLE(test_pastix test_pastix.F90)
  TARGET_LINK_LIBRARIES(test_pastix sll_sparse_matrix
                                    sll_collective
                                    sll_memory
                                    ${PASTIX_LIBRARIES} 
                                    ${PTSCOTCH_LIBRARIES}
                                    ${SCOTCH_LIBRARIES}
                                    ${BLAS_LIBRARIES})

ENDIF(PASTIX_FOUND)

IF(MUMPS_FOUND) #Simple test from the library 

  ADD_EXECUTABLE( test_mumps test_mumps.F90 )
  TARGET_LINK_LIBRARIES( test_mumps sll_sparse_matrix 
                           sll_collective ${MUMPS_LIBRARIES} )
  ADD_TEST( NAME mumps COMMAND test_mumps )
  SET_TESTS_PROPERTIES( mumps PROPERTIES PASS_REGULAR_EXPRESSION "PASSED" )

ENDIF(MUMPS_FOUND)
