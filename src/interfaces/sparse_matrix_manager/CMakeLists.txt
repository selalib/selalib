IF(MUMPS_ENABLED)
  FIND_PACKAGE(MUMPS)
ENDIF(MUMPS_ENABLED)

IF(PASTIX_ENABLED)
  FIND_PACKAGE(SCOTCH)
  FIND_PACKAGE(PASTIX)
ENDIF(PASTIX_ENABLED)

IF(SUITESPARSE_ENABLED)
  FIND_PACKAGE(SUITESPARSE)
ENDIF(SUITESPARSE_ENABLED)

SET(SPM_SOURCES sll_m_qsort_partition.F90 
                sll_m_sparse_matrix.F90
                sll_m_sparse_matrix_mp.F90)

IF(SUITESPARSE_FOUND)

  INCLUDE_DIRECTORIES(${SUITESPARSE_INCLUDE_DIRS})
  SET(SPM_SOURCES ${SPM_SOURCES} sll_m_umfpack.F90 umf4_f77wrapper.c)

ENDIF(SUITESPARSE_FOUND)

IF(PASTIX_FOUND)

  #PasTiX is parallel
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  INCLUDE_DIRECTORIES(${PASTIX_INCLUDE_DIRS})
  SET(SPM_SOURCES ${SPM_SOURCES} sll_m_pastix.F90)

ENDIF(PASTIX_FOUND)

IF(MUMPS_FOUND)

  #MUMPS is parallel
  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  INCLUDE_DIRECTORIES(${MUMPS_INCLUDE_DIRS})
  SET(SPM_SOURCES ${SPM_SOURCES} sll_m_mumps.F90)

ENDIF(MUMPS_FOUND)

ADD_LIBRARY(sll_sparse_matrix STATIC ${SPM_SOURCES})

TARGET_LINK_LIBRARIES( sll_sparse_matrix sll_memory sll_assert)

IF(SUITESPARSE_FOUND)
  TARGET_LINK_LIBRARIES( sll_sparse_matrix ${SUITESPARSE_LIBRARIES}
                                           ${BLAS_LIBRARIES})
ENDIF()

IF(PASTIX_FOUND)
  TARGET_LINK_LIBRARIES( sll_sparse_matrix sll_collective
                                           ${PASTIX_LIBRARIES} 
                                           ${PTSCOTCH_LIBRARIES}
                                           ${SCOTCH_LIBRARIES}
                                           ${BLAS_LIBRARIES})
ENDIF(PASTIX_FOUND)

IF(MUMPS_FOUND)
  TARGET_LINK_LIBRARIES( sll_sparse_matrix sll_collective ${MUMPS_LIBRARIES} )
ENDIF(MUMPS_FOUND)

IF(BUILD_PARALUTION)
  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/paralution-prefix/src/paralution/src/)
  ADD_LIBRARY(sll_paralution paralution_fortran.cpp sll_m_paralution.F90)
ENDIF()


IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
