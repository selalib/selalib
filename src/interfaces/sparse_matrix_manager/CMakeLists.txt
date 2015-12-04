IF(PASTIX_FOUND)

  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  INCLUDE_DIRECTORIES(${PASTIX_INCLUDE_DIRS})
  ADD_LIBRARY(sll_sparse_matrix STATIC 
              sll_m_pastix.F90
              sll_m_qsort_partition.F90
  	        sll_m_sparse_matrix_pastix.F90
	        sll_m_sparse_matrix_mp.F90)

  TARGET_LINK_LIBRARIES( sll_sparse_matrix sll_collective
                                           sll_memory
                                           sll_assert
                                           ${PASTIX_LIBRARIES} 
                                           ${PTSCOTCH_LIBRARIES}
                                           ${SCOTCH_LIBRARIES}
                                           ${BLAS_LIBRARIES})
ELSEIF(SUITESPARSE_FOUND)

  INCLUDE_DIRECTORIES(${SUITESPARSE_INCLUDE_DIRS})

  ADD_LIBRARY(sll_sparse_matrix STATIC 
              sll_m_qsort_partition.F90
  	        sll_m_sparse_matrix_umfpack.F90
	        sll_m_sparse_matrix_mp.F90
              sll_m_umfpack.F90
              umf4_f77wrapper.c)

  TARGET_LINK_LIBRARIES( sll_sparse_matrix sll_memory
                                           sll_assert
                                           ${SUITESPARSE_LIBRARIES}
                                           ${BLAS_LIBRARIES})


ELSE()

  ADD_LIBRARY(sll_sparse_matrix STATIC 
              sll_m_qsort_partition.F90
 	        sll_m_sparse_matrix.F90
	        sll_m_sparse_matrix_mp.F90 )

  TARGET_LINK_LIBRARIES( sll_sparse_matrix sll_memory
                                           sll_assert)


ENDIF()


IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
