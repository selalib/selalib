SET(SPM_FILES getdia.F90 
              getu.F90 
              getl.F90 
              lsol.F90 
              remove_zeros.F90  
              basic_cg.F90 
              gauss_seidel.F90 
              spm_def.F90 
              spm_tools.F90) 


MESSAGE(STATUS "SUITESPARSE_FOUND:${SUITESPARSE_FOUND}")

IF(PASTIX_FOUND)
   LIST(APPEND SPM_FILES spm_murge.F90)
ELSEIF(SUITESPARSE_FOUND)
   LIST(APPEND SPM_FILES spm_umfpack.F90)
ElSE()
   LIST(APPEND SPM_FILES spm_basic.F90)
ENDIF()

ADD_LIBRARY(spm STATIC ${SPM_FILES})


set(spm_mods basic_cg_mod.mod gauss_seidel_mod.mod spm_def.mod spm.mod)
install(TARGETS spm DESTINATION lib)
install(FILES ${spm_mods} DESTINATION include)

ADD_LIBRARY(sll_sparse_matrix STATIC sll_sparse_matrix_module.F90)

TARGET_LINK_LIBRARIES( sll_sparse_matrix 
  sll_working_precision 
  sll_memory
  sll_assertion)

ADD_EXECUTABLE( test_sparse_matrix unit_test_sparse_matrix.F90)

TARGET_LINK_LIBRARIES( test_sparse_matrix sll_sparse_matrix)


