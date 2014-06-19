IF( NOT STDF95)

  ADD_LIBRARY (sll_fcisl STATIC
		fcisl.F90)

  TARGET_LINK_LIBRARIES(sll_fcisl sll_utilities sll_advection)

  ADD_EXECUTABLE(test_fcisl unit_test.F90)

  TARGET_LINK_LIBRARIES(test_fcisl 
    sll_fcisl
    sll_logical_meshes 
    sll_advection
    sll_constants)

ENDIF(NOT STDF95)
