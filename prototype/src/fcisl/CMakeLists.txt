  ADD_LIBRARY (sll_fcisl STATIC
		fcisl.F90)

  TARGET_LINK_LIBRARIES(sll_fcisl sll_advection sll_utilities)

  ADD_EXECUTABLE(test_fcisl unit_test.F90)

  TARGET_LINK_LIBRARIES(test_fcisl 
    sll_fcisl
    sll_logical_meshes 
    sll_advection
    sll_constants)
