IF( NOT STDF95)

  SET( ADVECTION_SRC_FILES 
	sll_advection_1d_periodic.F90
	sll_advection_1d_non_uniform_cubic_splines.F90
	sll_advection_1d_BSL.F90
	sll_advection_1d_CSL.F90
	sll_advection_1d_PSM.F90
	sll_advection_1d_base.F90
	sll_advection_2d_BSL.F90
	sll_advection_2d_tensor_product.F90
	sll_advection_2d_CSL.F90
	sll_advection_2d_base.F90)


  ADD_LIBRARY(
    sll_advection
    STATIC
    ${ADVECTION_SRC_FILES})

  TARGET_LINK_LIBRARIES( 
    sll_advection
    sll_characteristics
    sll_interpolators
    sll_periodic_interp
    sll_splines
    sll_utilities
  )

  ADD_EXECUTABLE(
    test_advection_1d_periodic
    unit_test_advection_1d_periodic.F90)

  TARGET_LINK_LIBRARIES(
    test_advection_1d_periodic
    sll_advection)

  ADD_EXECUTABLE(
    test_advection_1d_non_uniform_cubic_splines
    unit_test_advection_1d_non_uniform_cubic_splines.F90)

  TARGET_LINK_LIBRARIES(
    test_advection_1d_non_uniform_cubic_splines
    sll_advection)

  ADD_EXECUTABLE(
    test_advection_1d_BSL
    unit_test_advection_1d_BSL.F90)

  TARGET_LINK_LIBRARIES(
    test_advection_1d_BSL
    sll_advection)

  ADD_EXECUTABLE(
    test_advection_1d_CSL
    unit_test_advection_1d_CSL.F90)

  TARGET_LINK_LIBRARIES(
    test_advection_1d_CSL
    sll_advection)

  ADD_EXECUTABLE(
    test_advection_1d_PSM
    unit_test_advection_1d_PSM.F90)

  TARGET_LINK_LIBRARIES(
    test_advection_1d_PSM
    sll_advection)




  ADD_EXECUTABLE(
    test_advection_2d_BSL
    unit_test_advection_2d_BSL.F90)

  TARGET_LINK_LIBRARIES(
    test_advection_2d_BSL
    sll_advection)

  ADD_EXECUTABLE(
    test_advection_2d_tensor_product
    unit_test_advection_2d_tensor_product.F90)

  TARGET_LINK_LIBRARIES(
    test_advection_2d_tensor_product
    sll_advection)


  ADD_EXECUTABLE(
    test_advection_2d_CSL
    unit_test_advection_2d_CSL.F90)

  TARGET_LINK_LIBRARIES(
    test_advection_2d_CSL
    sll_advection)

    
  ADD_EXECUTABLE(
    test_aligned_translation_2d
    aligned_translation_2d.F90)

  TARGET_LINK_LIBRARIES(
    test_aligned_translation_2d
    sll_advection
    sll_utilities)

  ADD_EXECUTABLE(
    test_rotation_2d
    rotation_2d.F90)

  TARGET_LINK_LIBRARIES(
    test_rotation_2d
    sll_advection)




ENDIF(NOT STDF95)
