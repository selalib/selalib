
IF(FULL_FORTRAN2003)
  ADD_LIBRARY( sll_linear_solvers STATIC 
    sll_m_tridiagonal.F90
    sll_m_penta_diagonal.F90
    sll_m_choleski.F90
    sol.f muls.f 
    sll_m_pivotbande.F90
    sll_m_linear_operator_base.F90
    sll_m_linear_operator_matrix_dense.F90
    sll_m_conjugate_gradient.F90
    )
ELSE()
  ADD_LIBRARY( sll_linear_solvers STATIC 
    sll_m_tridiagonal.F90
    sll_m_penta_diagonal.F90
    sll_m_choleski.F90
    sol.f muls.f 
    sll_m_pivotbande.F90
    )
ENDIF()

IF(Fortran_COMPILER MATCHES IBM)
   SET_SOURCE_FILES_PROPERTIES( sol.f muls.f PROPERTIES COMPILE_FLAGS '-qfixed' )
ENDIF()
IF(FULL_FORTRAN2003)
  TARGET_LINK_LIBRARIES( sll_linear_solvers sll_vector_space sll_errors sll_memory sll_assert )
ELSE()
  TARGET_LINK_LIBRARIES( sll_linear_solvers sll_errors sll_memory sll_assert )
ENDIF()

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
