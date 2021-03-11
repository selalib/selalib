if(MUMPS_ENABLED)
  find_package(MUMPS)
endif(MUMPS_ENABLED)

if(PASTIX_ENABLED)
  find_package(SCOTCH)
  find_package(PASTIX)
endif(PASTIX_ENABLED)

if(SUITESPARSE_ENABLED)
  find_package(SUITESPARSE)
endif(SUITESPARSE_ENABLED)

set(SPM_SOURCES sll_m_qsort_partition.F90 sll_m_sparse_matrix.F90
                sll_m_sparse_matrix_mp.F90)

if(SUITESPARSE_FOUND)

  set(INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS})
  set(SPM_SOURCES ${SPM_SOURCES} sll_m_umfpack.F90 umf4_f77wrapper.c)

endif(SUITESPARSE_FOUND)

if(PASTIX_FOUND)

  # PasTiX is parallel
  set(INCLUDE_DIRS ${PASTIX_INCLUDE_DIRS})
  set(SPM_SOURCES ${SPM_SOURCES} sll_m_pastix.F90)

endif(PASTIX_FOUND)

if(MUMPS_FOUND)

  # MUMPS is parallel
  set(INCLUDE_DIRS ${MUMPS_INCLUDE_DIRS})
  set(SPM_SOURCES ${SPM_SOURCES} sll_m_mumps.F90)

endif(MUMPS_FOUND)

add_library(sll_sparse_matrix STATIC ${SPM_SOURCES})

target_include_directories(sll_sparse_matrix PRIVATE ${INCLUDE_DIRS})
target_link_libraries(sll_sparse_matrix sll_memory sll_assert)

if(SUITESPARSE_FOUND)
  target_link_libraries(sll_sparse_matrix ${SUITESPARSE_LIBRARIES}
                        ${BLAS_LIBRARIES})
endif()

if(PASTIX_FOUND)
  target_link_libraries(
    sll_sparse_matrix sll_collective ${PASTIX_LIBRARIES} ${PTSCOTCH_LIBRARIES}
    ${SCOTCH_LIBRARIES} ${BLAS_LIBRARIES})
endif(PASTIX_FOUND)

if(MUMPS_FOUND)
  target_link_libraries(sll_sparse_matrix sll_collective ${MUMPS_LIBRARIES})
endif(MUMPS_FOUND)

if(BUILD_PARALUTION)
  include_directories(${CMAKE_BINARY_DIR}/paralution-prefix/src/paralution/src/)
  add_library(sll_paralution paralution_fortran.cpp sll_m_paralution.F90)
endif()

if(BUILD_TESTING)
  add_subdirectory(testing)
endif(BUILD_TESTING)
