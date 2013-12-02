IF( NOT STDF95)

  SET( ADVECTION_SRC_FILES 
	sll_advection_1d_periodic.F90
	sll_advection_1d_base.F90
	sll_advection_2d_BSL.F90
	sll_advection_2d_base.F90)


  ADD_LIBRARY(
    sll_advection
    STATIC
    ${ADVECTION_SRC_FILES})

  TARGET_LINK_LIBRARIES( 
    sll_advection
    sll_periodic_interp
    sll_characteristics
    sll_interpolators
  )

  ADD_EXECUTABLE(
    test_advection_1d_periodic
    unit_test_advection_1d_periodic.F90)

  TARGET_LINK_LIBRARIES(
    test_advection_1d_periodic
    sll_advection)

  ADD_EXECUTABLE(
    test_advection_2d_BSL
    unit_test_advection_2d_BSL.F90)

  TARGET_LINK_LIBRARIES(
    test_advection_2d_BSL
    sll_advection)
    
    
    

  ADD_EXECUTABLE(
    test_aligned_translation_2d
    aligned_translation_2d.F90)

  TARGET_LINK_LIBRARIES(
    test_aligned_translation_2d
    sll_utilities
    sll_advection)



ENDIF(NOT STDF95)
