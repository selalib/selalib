ADD_LIBRARY( sll_linear_solvers STATIC 
  sll_m_tridiagonal.F90
  sll_m_penta_diagonal.F90
  sll_m_choleski.F90
  sol.f muls.f 
  sll_m_pivotbande.F90
)

IF(Fortran_COMPILER MATCHES IBM)
   SET_SOURCE_FILES_PROPERTIES( sol.f muls.f PROPERTIES COMPILE_FLAGS '-qfixed' )
ENDIF()

TARGET_LINK_LIBRARIES( sll_linear_solvers sll_errors sll_memory sll_assert )

IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
