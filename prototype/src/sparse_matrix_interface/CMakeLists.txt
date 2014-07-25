SET(SPM_FILES getdia.F90 
              getu.F90 
              getl.F90 
              lsol.F90 
              remove_zeros.F90  
              basic_cg.F90 
              gauss_seidel.F90 
              spm_def.F90 
              spm_tools.F90
   ) 


IF(PASTIX_FOUND)
#   LIST(APPEND SPM_FILES spm_murge.F90)
ELSEIF(SUITESPARSE_FOUND)
   LIST(APPEND SPM_FILES spm_umfpack.F90)
ElSE()
   LIST(APPEND SPM_FILES spm_basic.F90)
ENDIF()

ADD_LIBRARY(spm STATIC ${SPM_FILES})

#set(spm_mods basic_cg_mod.mod gauss_seidel_mod.mod spm_def.mod spm.mod)
#install(TARGETS spm DESTINATION lib)
#install(FILES ${spm_mods} DESTINATION include)

IF(PASTIX_FOUND)
  ADD_LIBRARY(sll_sparse_matrix STATIC 
  	      sll_sparse_matrix_pastix_module.F90
	      sll_sparse_matrix_mp_module.F90)
ELSEIF(SUITESPARSE_FOUND)
  ADD_LIBRARY(sll_sparse_matrix STATIC 
  	      sll_sparse_matrix_umfpack_module.F90
	      sll_sparse_matrix_mp_module.F90)
ElSE()
  ADD_LIBRARY(sll_sparse_matrix STATIC 
  	      sll_sparse_matrix_module.F90
	      sll_sparse_matrix_mp_module.F90 )
ENDIF()

TARGET_LINK_LIBRARIES( sll_sparse_matrix 
  sll_working_precision 
  sll_memory
  sll_assertion)

IF(PASTIX_FOUND)

  TARGET_LINK_LIBRARIES( sll_sparse_matrix 
    sll_pastix 
    ${BLAS_LIBRARIES})
    
ELSEIF(SUITESPARSE_FOUND)

  ADD_LIBRARY(sll_umfpack_wrapper STATIC mod_umfpack.F90 umf4_f77wrapper.c)

  TARGET_LINK_LIBRARIES( sll_sparse_matrix 
    sll_umfpack_wrapper 
    ${SUITESPARSE_LIBRARIES}
    ${BLAS_LIBRARIES})

ENDIF()




ADD_EXECUTABLE( test_sparse_matrix unit_test_sparse_matrix.F90)

#to be changed o that it compiles also without mpi
TARGET_LINK_LIBRARIES( test_sparse_matrix sll_sparse_matrix sll_collective)


ADD_LIBRARY(sll_csr_one_constraint STATIC sll_csr_one_constraint_module.F90)
TARGET_LINK_LIBRARIES(sll_csr_one_constraint sll_utilities sll_memory)
ADD_EXECUTABLE(test_csr_one_constraint unit_test_csr_one_constraint.F90)
TARGET_LINK_LIBRARIES(test_csr_one_constraint sll_csr_one_constraint)

#=======
#set(spm_mods basic_cg_mod.mod gauss_seidel_mod.mod spm_def.mod spm.mod)
#install(TARGETS spm DESTINATION lib)
#install(FILES ${spm_mods} DESTINATION include)

