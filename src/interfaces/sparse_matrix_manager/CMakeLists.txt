SET(SPM_SOURCES sll_m_qsort_partition.F90 
                sll_m_sparse_matrix.F90
                sll_m_sparse_matrix_mp.F90)

IF(SUITESPARSE_FOUND)

  INCLUDE_DIRECTORIES(${SUITESPARSE_INCLUDE_DIRS})
  ADD_DEFINITIONS(-DUMFPACK)

  SET(SPM_SOURCES ${SPM_SOURCES}
              sll_m_umfpack.F90
              umf4_f77wrapper.c)

ENDIF(SUITESPARSE_FOUND)

#IF(PASTIX_FOUND)
#
#  #PasTiX is parallel
#  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
#  INCLUDE_DIRECTORIES(${PASTIX_INCLUDE_DIRS})
#
#  ADD_LIBRARY(sll_sparse_matrix STATIC 
#              sll_m_pastix.F90
#              sll_m_sparse_matrix_pastix.F90
#              ${SPM_SOURCES})
#
#  TARGET_LINK_LIBRARIES( sll_sparse_matrix sll_collective
#                                           sll_memory
#                                           sll_assert
#                                           ${PASTIX_LIBRARIES} 
#                                           ${PTSCOTCH_LIBRARIES}
#                                           ${SCOTCH_LIBRARIES}
#                                           ${BLAS_LIBRARIES})

ADD_LIBRARY(sll_sparse_matrix STATIC ${SPM_SOURCES})

TARGET_LINK_LIBRARIES( sll_sparse_matrix sll_memory sll_assert)

IF(SUITESPARSE_FOUND)

  TARGET_LINK_LIBRARIES( sll_sparse_matrix ${SUITESPARSE_LIBRARIES}
                                           ${BLAS_LIBRARIES})


ENDIF()

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
