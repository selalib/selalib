SET(SPM_FILES getdia.F90 
              getu.F90 
              getl.F90 
              lsol.F90 
              remove_zeros.F90  
              basic_cg.F90 
              gauss_seidel.F90 
              spm_def.F90 
              spm_tools.F90) 

IF(SUITESPARSE_FOUND)
   LIST(APPEND SPM_FILES spm_umfpack.F90)
ELSE()
   LIST(APPEND SPM_FILES spm_basic.F90)
ENDIF()

ADD_LIBRARY(spm STATIC ${SPM_FILES})

IF(PASTIX_FOUND)

  SET(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  INCLUDE_DIRECTORIES(${PASTIX_INCLUDE_DIRS})
  ADD_LIBRARY(sll_sparse_matrix STATIC 
              sll_sparse.F90
              sll_pastix.F90
  	        sll_sparse_matrix_pastix_module.F90
	        sll_sparse_matrix_mp_module.F90)


ELSEIF(SUITESPARSE_FOUND)

  INCLUDE_DIRECTORIES(${SUITESPARSE_INCLUDE_DIRS})
  ADD_LIBRARY(sll_sparse_matrix STATIC 
  	        sll_sparse_matrix_umfpack_module.F90
	        sll_sparse_matrix_mp_module.F90)

ElSE()

  ADD_LIBRARY(sll_sparse_matrix STATIC 
  	        sll_sparse_matrix_module.F90
	        sll_sparse_matrix_mp_module.F90 )

ENDIF()

TARGET_LINK_LIBRARIES( sll_sparse_matrix sll_working_precision 
                                         sll_memory
                                         sll_assert)

IF(PASTIX_FOUND)

  TARGET_LINK_LIBRARIES( sll_sparse_matrix ${PASTIX_LIBRARIES} 
                                           ${PTSCOTCH_LIBRARIES}
                                           ${SCOTCH_LIBRARIES}
                                           ${BLAS_LIBRARIES}
                                           sll_collective)

ELSEIF(SUITESPARSE_FOUND)

  ADD_LIBRARY(sll_umfpack_wrapper STATIC mod_umfpack.F90 umf4_f77wrapper.c)

  TARGET_LINK_LIBRARIES( sll_sparse_matrix sll_umfpack_wrapper 
                                           ${SUITESPARSE_LIBRARIES}
                                           ${BLAS_LIBRARIES})

ENDIF()

ADD_LIBRARY(sll_csr_one_constraint STATIC sll_csr_one_constraint_module.F90)
TARGET_LINK_LIBRARIES(sll_csr_one_constraint sll_utilities sll_memory)


# Tests
#PN Do not know why sll_collective is necessary here...
IF(BUILD_TESTING AND MPI_ENABLED)

  ADD_EXECUTABLE( test_sparse_matrix unit_test_sparse_matrix.F90 )
  TARGET_LINK_LIBRARIES( test_sparse_matrix sll_sparse_matrix sll_collective )

  ADD_EXECUTABLE( test_csr_one_constraint unit_test_csr_one_constraint.F90 )
  TARGET_LINK_LIBRARIES( test_csr_one_constraint sll_csr_one_constraint )

  IF(Fortran_COMPILER MATCHES IBM)
   SET_SOURCE_FILES_PROPERTIES(dPackgmres.f PROPERTIES COMPILE_FLAGS -qfixed)
  ENDIF(Fortran_COMPILER MATCHES IBM)
  ADD_EXECUTABLE( test_gmres test_gmres_driver.F90 dPackgmres.f )
  TARGET_LINK_LIBRARIES( test_gmres ${BLAS_LIBRARIES} )

ENDIF(BUILD_TESTING AND MPI_ENABLED)
