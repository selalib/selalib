IF(PASTIX_FOUND)

  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  INCLUDE_DIRECTORIES(${PASTIX_INCLUDE_DIRS})
  ADD_LIBRARY(sll_sparse_matrix STATIC 
              sll_pastix.F90
  	        sll_sparse_matrix_pastix_module.F90
	        sll_sparse_matrix_mp_module.F90)

  TARGET_LINK_LIBRARIES( sll_sparse_matrix sll_collective
                                           sll_working_precision 
                                           sll_memory
                                           sll_assert
                                           ${PASTIX_LIBRARIES} 
                                           ${PTSCOTCH_LIBRARIES}
                                           ${SCOTCH_LIBRARIES}
                                           ${BLAS_LIBRARIES})
ELSEIF(SUITESPARSE_FOUND)

  INCLUDE_DIRECTORIES(${SUITESPARSE_INCLUDE_DIRS})

  ADD_LIBRARY(sll_sparse_matrix STATIC 
  	        sll_sparse_matrix_umfpack_module.F90
	        sll_sparse_matrix_mp_module.F90)

  ADD_LIBRARY(sll_umfpack_wrapper STATIC mod_umfpack.F90 umf4_f77wrapper.c)

  TARGET_LINK_LIBRARIES( sll_sparse_matrix sll_working_precision
                                           sll_memory
                                           sll_assert
                                           sll_umfpack_wrapper 
                                           ${SUITESPARSE_LIBRARIES}
                                           ${BLAS_LIBRARIES})


ELSE()

  ADD_LIBRARY(sll_sparse_matrix STATIC 
  	        sll_sparse_matrix_module.F90
	        sll_sparse_matrix_mp_module.F90 )

ENDIF()


#PN Do not know why sll_collective is necessary here...
IF(BUILD_TESTING)

  IF(MPI_ENABLED)
    ADD_EXECUTABLE( test_sparse_matrix unit_test_sparse_matrix.F90 )
    TARGET_LINK_LIBRARIES( test_sparse_matrix sll_sparse_matrix sll_collective )
  ENDIF(MPI_ENABLED)

  ADD_LIBRARY(sll_csr_one_constraint STATIC sll_csr_one_constraint_module.F90)
  TARGET_LINK_LIBRARIES(sll_csr_one_constraint sll_utilities sll_memory)
  ADD_EXECUTABLE( test_csr_one_constraint unit_test_csr_one_constraint.F90 )
  TARGET_LINK_LIBRARIES( test_csr_one_constraint sll_csr_one_constraint )

  IF(Fortran_COMPILER MATCHES IBM)
   SET_SOURCE_FILES_PROPERTIES(dPackgmres.f PROPERTIES COMPILE_FLAGS -qfixed)
  ENDIF(Fortran_COMPILER MATCHES IBM)
  ADD_EXECUTABLE( test_gmres test_gmres_driver.F90 dPackgmres.f )
  TARGET_LINK_LIBRARIES( test_gmres ${BLAS_LIBRARIES} )

  ADD_EXECUTABLE(test_csr_matrix unit_test_csr_matrix.F90)
  TARGET_LINK_LIBRARIES( test_csr_matrix sll_sparse_matrix dagmg_seq)

ENDIF(BUILD_TESTING)

#IF(BUILD_PARALUTION)
#
#  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/paralution/src)
#  ADD_EXECUTABLE(paralution_solver 
#                  ${CMAKE_BINARY_DIR}/paralution/src/plug-ins/paralution_fortran.cpp
#                  sll_paralution.f90)
#
#ENDIF(BUILD_PARALUTION)
