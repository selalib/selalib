SET( SRC_UTILITIES
  sll_utilities.F90
  tridiagonal.F90
  sll_toeplitz_penta_diagonal.F90
  sll_primitives.F90
  choleski.F90
  fornberg.F90
  kiss.F90
  sol.f muls.f 
  banfac.f90 
  banslv.f90  
  pivotbande.F90
)

IF(Fortran_COMPILER MATCHES IBM)
   SET_SOURCE_FILES_PROPERTIES( sol.f muls.f 
     PROPERTIES COMPILE_FLAGS '-qfixed' )
ENDIF()

ADD_LIBRARY( sll_utilities STATIC ${SRC_UTILITIES} )
TARGET_LINK_LIBRARIES( sll_utilities
  sll_errors sll_memory sll_assert sll_working_precision )


IF(BUILD_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(BUILD_TESTING)
